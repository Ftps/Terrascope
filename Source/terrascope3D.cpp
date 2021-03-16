#include "terrascope3D.hpp"

#define OP 0.05

Planet3D::Planet3D(const double &R, const double& r_max, const double& obf, const std::function<dddd>& n) : R(R), r_max(r_max), obf(obf), n(n) {}

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
		angle[PHI] = 2*M_PI*i/(n_phi-1) - M_PI/2;
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
			colorMap->data()->setCell(j, i, data[i][j]);
		}
	}

	plt->plotLayout()->addElement(0, 1, colorScale);
	colorScale->setType(QCPAxis::atRight);
	colorMap->setColorScale(colorScale);
	colorScale->axis()->setLabel("Number of rays");
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




double gx3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x+h, y, z) - f(x-h, y, z))/(2*h);
}

double gy3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x, y+h, z) - f(x, y-h, z))/(2*h);
}

double gz3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x, y, z+h) - f(x, y, z-h))/(2*h);
}


std::array<double,2> rayTracing(const Planet3D& p, const double& L, const std::array<double,2>& a, const double& dz)
{
	double vx, vy, vz, x, y, z, n2, gx, gy, dvx, dvy, r;
	double beta = 1/(1 - p.obf);
	double R = p.R*p.R, r_max = p.r_max*p.r_max;
	double ct = cos(a[TET]), st = sin(a[TET]), cp = cos(a[PHI]);
	double k = (sq(beta)-1)*sq(st*cp) + 1;
	double det = sq(L*ct) + k*(r_max - L*L);

	if(det <= 0) return a;

	vx = st*cp;
	vy = st*sin(a[PHI]);
	vz = ct;

	k = L*ct - sqrt(det)/k;
	x = k*vx;
	y = k*vy;
	z = k*vz - L;

	do{
		n2 = 1/(vz*sq(p.n(x, y, z)));
		gx = gx3D(p.n, x, y, z);
		gy = gy3D(p.n, x, y, z);

		k = gx*vx + gy*vy + gz3D(p.n, x, y, z)*vz;
		dvx = n2*(gx - k*vx);
		dvy = n2*(gy - k*vy);

		vx += dvx*dz;
		vy += dvy*dz;
		vz = sqrt(1 - vx*vx - vy*vy);

		x += vx*dz/vz;
		y += vy*dz/vz;
		z += dz;

		r = sq(beta*x) + y*y + z*z;

		if(r < R) return {-100, -100};

	}while(r < r_max);

	return {atan(vx/vz), atan(vy/vz)};
}
