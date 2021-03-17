#include "drawCFM.hpp"

ImageGen::ImageGen(const Planet3D& p, const double& L, const int& N) : p(p), L(L), N(N)
{
	grid = new QGridLayout(this);
	plt = new QCustomPlot(this);

	drawCFM();

	plt->setFixedSize(900, 700);
	grid->addWidget(plt, 0, 0);
	this->setWindowTitle("Central Flash Map");
}

void ImageGen::drawCFM()
{
	int n_phi = 10*N;
	int n_tet = 10*n_phi, ray_conter = 0;
	std::array<double, 2> angle, ex;
	std::array<int, 2> ij;
	QCPColorMap *colorMap = new QCPColorMap(plt->xAxis, plt->yAxis);
	QCPColorScale *colorScale = new QCPColorScale(plt);
	QCPMarginGroup *marginGroup = new QCPMarginGroup(plt);
	double **data = createMatrix<double>(N+1);
	double max = acos(sqrt(1 - sq(p.r_max/L)));
	double min = acos(sqrt(1 - sq(p.R/L))), opt = 0;
	double dt = (max - min)/n_tet, h = 2*OP*max/N;
	double dphi = 1/((double)n_phi-1);
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> t;
	bool in_picture;

	plt->axisRect()->setupFullAxesBox(true);
	plt->xAxis->setLabel("x - km");
	plt->yAxis->setLabel("y - km");

	colorMap->data()->setSize(N+1, N+1);
	colorMap->data()->setRange(QCPRange(-OP*L*max, OP*L*max), QCPRange(-OP*L*max, OP*L*max));

	/*colorMap->data()->cellToCoord(N/2, N/2, &x, &y);
	Print(1.5*max);
	Print("x = " << x);
	Print("y = " << y);*/

	start = std::chrono::system_clock::now();
	for(int i = 0; i < n_phi; ++i){
		angle[PHI] = 2*M_PI*i*dphi - M_PI/2;
		angle[TET] = (opt)?opt:max;
		in_picture = false;
		do{
			++ray_conter;
			ex = rayTracing(p, L, angle);
			ij[X] = (int)(N*(ex[X] + OP*max + OP*h)/(2*OP*max));
			ij[Y] = (int)(N*(ex[Y] + OP*max + OP*h)/(2*OP*max));

			if(ij[X] >= 0 && ij[Y] >= 0 && ij[X] < N+1 && ij[Y] < N+1){
				if(!in_picture){
					in_picture = true;
					if(!opt) opt = angle[TET];
				}
				data[ij[X]][ij[Y]] += 1;
			}
			else if(in_picture) break;

			angle[TET] -= dt;
		}while(ex[X] != -100);
		Print(i+1 << " out of " << n_phi);
		Print("Ray count: " << ray_conter << std::endl);
	}
	end = std::chrono::system_clock::now();
	t = end - start;
	Print("Elapsed calculation time: " << t.count());

	for(int i = 0; i < N+1; ++i){
		for(int j = 0; j < N+1; ++j){
			colorMap->data()->setCell(j, i, data[i][j]*dphi*dt/sq(h));
		}
	}

	plt->plotLayout()->addElement(0, 1, colorScale);
	colorScale->setType(QCPAxis::atRight);
	colorMap->setColorScale(colorScale);
	colorScale->axis()->setLabel("Ray density");
	colorMap->setGradient(redGrad());
	colorMap->setInterpolate(true);
	colorMap->rescaleDataRange();
	plt->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
	colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
	plt->rescaleAxes();

	freeMatrix<double>(data, N+1);
}



QCPColorGradient redGrad()
{
	QCPColorGradient grad;
	QMap<double, QColor> c;

	c[0] = Qt::black;
	c[0.5] = Qt::red;
	c[0.8] = Qt::white;
	c[1] = Qt::white;

	grad.setColorStops(c);
	grad.setColorInterpolation(QCPColorGradient::ciRGB);

	return grad;
}
