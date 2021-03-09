#include "drawRay.hpp"

DrawRay::DrawRay(const Planet2D& p, const double& L, const double& a_entry, const int& type)
{
	int n_ray = 250;
	double y, h = (p.r_max - p.R)/n_ray;
	grid = new QGridLayout(this);
	plt = new QCustomPlot(this);

	if(type == 1){
		drawRayDensity(plt, p, Qt::red);

		plt->xAxis->setLabel("x / km");
		plt->yAxis->setLabel("log10(ray density)");
	}
	else if(type == 2){
		drawCircle(plt, p.r_max, p.obf, Qt::white, QColor(0, 0, 255, 40));
		drawCircle(plt, p.R, p.obf, QColorConstants::Svg::brown, QColorConstants::Svg::brown);

		for(int i = 0; i < n_ray; ++i){
			y = p.r_max - i*h;
			drawRayHor(plt, p, y, Qt::black);
		}

		plt->xAxis->setRange(-3.5*p.r_max, 3.5*p.r_max);
		plt->yAxis->setRange(-2*p.r_max, 2*p.r_max);
		plt->xAxis->setLabel("x / km");
		plt->yAxis->setLabel("y / km");
		plt->xAxis2->setVisible(true);
		plt->xAxis2->setTickLabels(false);
		plt->yAxis2->setVisible(true);
		plt->yAxis2->setTickLabels(false);
	}
	else if(type == 3){
		drawCircle(plt, p.r_max, p.obf, Qt::white, QColor(0, 0, 255, 40));
		drawCircle(plt, p.R, p.obf, QColorConstants::Svg::brown, QColorConstants::Svg::brown);
		drawRay(plt, p, L, a_entry, Qt::black);

		plt->xAxis->setRange(-3.5*p.r_max, 3.5*p.r_max);
		plt->yAxis->setRange(-2*p.r_max, 2*p.r_max);
		plt->xAxis->setLabel("x / km");
		plt->yAxis->setLabel("y / km");
		plt->xAxis2->setVisible(true);
		plt->xAxis2->setTickLabels(false);
		plt->yAxis2->setVisible(true);
		plt->yAxis2->setTickLabels(false);
	}

	plt->setFixedSize(1280, 720);
	grid->addWidget(plt, 0, 0);
	this->setWindowTitle("Draw Ray");
	//Print("Number of rays: " << n_ray);
}

void drawCircle(QCustomPlot* plt, const double& r, const double& obf, const QColor& line_c, const QColor& fill_c)
{
	QCPGraph* aux = plt->addGraph();
	QVector<double> x(CIRC_NUM), y1(CIRC_NUM), y2(CIRC_NUM), yl(CIRC_NUM);
	double rr = 1 - obf;
	double h = 2*r/(rr*(CIRC_NUM-1)), r2 = r*r;

	x[0] = -r/rr;
	y1[0] = 0;
	y2[0] = 0;

	for(int i = 1; i < CIRC_NUM-1; ++i){
		x[i] = x[i-1] + h;
		y1[i] = qSqrt(r2 - pwr(rr*x[i], 2));
		y2[i] = -y1[i];
	}

	x[CIRC_NUM-1] = r/rr;
	y1[CIRC_NUM-1] = 0;
	y2[CIRC_NUM-1] = 0;

	aux->setPen(line_c);
	aux->setBrush(fill_c);
	aux->setData(x, y1);

	aux = plt->addGraph();
	aux->setPen(line_c);
	aux->setBrush(fill_c);
	aux->setData(x, y2);

	aux = plt->addGraph();
	aux->setPen(fill_c);
	aux->setData(x, yl);
}

int drawRay(QCustomPlot* plt, const Planet2D& p, const double& L, const double& a_entry, const QColor& line_c)
{
	QCPGraph* aux = plt->addGraph();
	Ray2D r(p);
	int ex;

	ex = r.ray_tracer(L, a_entry);

	aux->setPen(line_c);
	aux->setData(QVector<double>(r.x.begin(), r.x.end()), QVector<double>(r.y.begin(), r.y.end()));

	return ex;
}

int drawRayHor(QCustomPlot* plt, const Planet2D& p, const double& Y, const QColor& line_c)
{
	QCPGraph* aux = plt->addGraph();
	Ray2D r(p);
	int ex;

	ex = r.ray_tracer_hor(Y);

	aux->setPen(line_c);
	aux->setData(QVector<double>(r.x.begin(), r.x.end()), QVector<double>(r.y.begin(), r.y.end()));

	return ex;
}

void drawRayDensity(QCustomPlot *plt, const Planet2D& p, const QColor& line_c)
{
	QCPGraph* aux = plt->addGraph();
	double h_min = p.R, dh = 0.01, a_min, h;
	QVector<double> x, dx, hh;
	int n;

	do{
		h_min += dh;
		a_min = ray_tracer2D_hor(p.n, p.R, p.r_max, h_min, 1);
	}while(a_min == -2);

	n = floor((p.r_max - h_min)/dh);

	Print("Number of rays: " << n);

	x.resize(n);
	dx.resize(n-1);
	hh.resize(n-1);

	x[0] = a_min;
	h = h_min;

	for(int i = 1; i < n; ++i){
		h += dh;
		x[i] = ray_tracer2D_hor(p.n, p.R, p.r_max, p.obf, h, 5);
		dx[i-1] = log10(dh/((x[i] - x[i-1])));
		hh[i-1] = 0.5*(x[i] - x[i-1]);
	}

	Print(a_min);

	aux->setPen(line_c);
	aux->setData(hh, dx);
	aux->rescaleAxes();
}
