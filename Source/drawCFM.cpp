#include "drawCFM.hpp"

ImageGen::ImageGen(const Planet3D& p, const double& l, const double& L, const double& S, const int& N) : p(p), L(L), S(S), l(l), N(N)
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
	int n_tet = 10000*N/S, ray_counter = 0;
	std::array<double, 2> angle;
	std::array<int, 2> ij;
	std::array<double, 3> ex;
	QCPColorMap *colorMap = new QCPColorMap(plt->xAxis, plt->yAxis);
	QCPColorScale *colorScale = new QCPColorScale(plt);
	QCPMarginGroup *marginGroup = new QCPMarginGroup(plt);
	double **data = createMatrix<double>(N+1);
	double max = acos(sqrt(1 - sq(p.r_max/L)));
	double min = acos(sqrt(1 - sq(p.R/L))), opt = 0;
	double dt = (max - min)/n_tet, scale = S/L;
	double dphi = 1/((double)n_phi-1), hh = 0, h = scale/N;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> t;
	bool in_picture;

	plt->axisRect()->setupFullAxesBox(true);
	plt->xAxis->setLabel("y - km");
	plt->yAxis->setLabel("x - km");
	plt->xAxis->setScaleRatio(plt->yAxis, 1);

	colorMap->data()->setSize(N+1, N+1);
	colorMap->data()->setRange(QCPRange(-S, S), QCPRange(-S, S));

	/*colorMap->data()->cellToCoord(N/2, N/2, &x, &y);
	Print(1.5*max);
	Print("x = " << x);
	Print("y = " << y);*/

	p.updateSigma(l);
	Print(p.o);
	LOG
	getchar();
	start = std::chrono::system_clock::now();
	for(int i = 0; i < n_phi; ++i){
		angle[PHI] = 2*M_PI*i*dphi - M_PI/2;
		if(p.dev) angle[TET] = ((opt)?(0.97*opt+0.03*max):max) - hh*cos(angle[PHI])/L;
		else angle[TET] = (opt)?opt:max - hh*cos(angle[PHI])/L;
		in_picture = false;
		do{
			++ray_counter;
			ex = rayTracing(p, hh, L, angle);

			ij[X] = (int)(N*(ex[X] + scale + h + hh/L)/(2*scale));
			ij[Y] = (int)(N*(ex[Y] + scale + h)/(2*scale));

			if(ij[X] >= 0 && ij[Y] >= 0 && ij[X] < N+1 && ij[Y] < N+1){
				if(!in_picture){
					in_picture = true;
					if(!(opt))/* || p.dev))*/ {opt = angle[TET];}
				}
				data[ij[X]][ij[Y]] += ex[BRT];
			}
			else if(in_picture) break;

			angle[TET] -= dt;
		}while(ex[X] != -100);
		Print(i+1 << " out of " << n_phi);
		Print("Ray count: " << ray_counter << std::endl);
	}
	end = std::chrono::system_clock::now();
	t = end - start;
	Print("Elapsed calculation time: " << t.count());

	for(int i = 0; i < N+1; ++i){
		for(int j = 0; j < N+1; ++j){
			colorMap->data()->setCell(i, j, data[j][i]*sq(L*h)/(sq(p.r_max)*dphi*dt));
		}
	}

	plt->plotLayout()->addElement(0, 1, colorScale);
	colorScale->setType(QCPAxis::atRight);
	colorMap->setColorScale(colorScale);
	colorScale->axis()->setLabel("Amplification Factor");
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
