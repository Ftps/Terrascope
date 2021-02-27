#include "drawRay.hpp"

DrawRay::DrawRay(const Planet2D& p, const double& L, const double& a_entry)
{
	grid = new QGridLayout(this);
	plt = new QCustomPlot(this);

	drawCircle(plt, p.r_max, Qt::white, QColor(0, 0, 255, 40));
	drawCircle(plt, p.R, QColorConstants::Svg::brown, QColorConstants::Svg::brown);
	drawRay(plt, p, L, a_entry, Qt::black);

	plt->xAxis->setRange(-3.5*p.r_max, 3.5*p.r_max);
	plt->yAxis->setRange(-2*p.r_max, 2*p.r_max);
	plt->xAxis->setLabel("x / km");
	plt->yAxis->setLabel("y / km");
	plt->xAxis2->setVisible(true);
	plt->xAxis2->setTickLabels(false);
	plt->yAxis2->setVisible(true);
	plt->yAxis2->setTickLabels(false);

	grid->addWidget(plt, 0, 0);
	resize(QDesktopWidget().availableGeometry(this).size() * 0.7);
	this->setWindowTitle("Draw Ray");
}

void drawCircle(QCustomPlot* plt, const double& r, const QColor& line_c, const QColor& fill_c)
{
	QCPGraph* aux = plt->addGraph();
	QVector<double> x(CIRC_NUM), y1(CIRC_NUM), y2(CIRC_NUM), yl(CIRC_NUM);
	double h = 2*r/(CIRC_NUM-1), r2 = r*r;

	x[0] = -r;
	y1[0] = 0;
	y2[0] = 0;

	for(int i = 1; i < CIRC_NUM-1; ++i){
		x[i] = x[i-1] + h;
		y1[i] = qSqrt(r2 - x[i]*x[i]);
		y2[i] = -y1[i];
	}

	x[CIRC_NUM-1] = r;
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

void drawRay(QCustomPlot* plt, const Planet2D& p, const double& L, const double& a_entry, const QColor& line_c)
{
	QCPGraph* aux = plt->addGraph();
	Ray2D r(p);

	r.ray_tracer(L, a_entry);

	aux->setPen(line_c);
	aux->setData(QVector<double>(r.x.begin(), r.x.end()), QVector<double>(r.y.begin(), r.y.end()));
}
