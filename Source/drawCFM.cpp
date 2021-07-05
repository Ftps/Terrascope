#include "drawCFM.hpp"

ImageGen::ImageGen(Planet3D& p, const double& l, const double& L, const double& S, const int& N, const double& h) : L(L), S(S), l(l), h(h), N(N)
{
	grid = new QGridLayout(this);
	plt = new QCustomPlot(this);

	drawCFM(p);

	plt->setFixedSize(900, 700);
	grid->addWidget(plt, 0, 0);
	this->setWindowTitle("Central Flash Map");
}

void ImageGen::drawCFM(Planet3D& p)
{
	QCPColorMap *colorMap = new QCPColorMap(plt->xAxis, plt->yAxis);
	QCPColorScale *colorScale = new QCPColorScale(plt);
	QCPMarginGroup *marginGroup = new QCPMarginGroup(plt);
	FlashMap *map;

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
	//map = mapGen2(p, N, S, L, h);
	Print("h = " << h);
	map = mapThread(p, N, S, L, h, 4);

	for(int i = 0; i < N+1; ++i){
		for(int j = 0; j < N+1; ++j){
			colorMap->data()->setCell(i, j, map->map[j][i]/map->Int2);
		}
	}

	plt->plotLayout()->addElement(0, 1, colorScale);
	colorScale->setType(QCPAxis::atRight);
	colorMap->setColorScale(colorScale);
	colorScale->axis()->setLabel("Amplification Factor");
	colorMap->setGradient(redGrad());
	colorMap->setInterpolate(false);
	colorMap->rescaleDataRange();
	plt->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
	colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
	plt->rescaleAxes();

	delete map;
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
