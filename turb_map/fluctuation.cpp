#include "fluctuation.hpp"

Turbulence_Map::Turbulence_Map(const MapOpt& m) : MapOpt(m), Ysh(l_max), min(300), max(-300), Yml(l_max), E(l_max)
{
	//double max_t = 0, max_p = 0, h = 1e-3, gtet, gphi;
	double max_t = -300, min_t = 300;
	double tet, phi;
	double dtet = 2*max_tet/w, dphi = 2*M_PI/l;
	double aux_coef;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<double> dist;
	image = new QImage(l, w, QImage::Format_RGB32);

	for(int i = 0; i < l_max; ++i){
		Yml[i].resize(2*(i+1)+1);
		dist = std::normal_distribution<double>(0, sqrt(n_ref*CLM/(2*(i+1)+1)));
		aux_coef = std::pow(i+1, -alpha/2);
		E[i] = 0;
		for(int j = 0; j < 2*(i+1)+1; ++j){
			Yml[i][j] = dist(gen)*aux_coef;
			E[i] += sq(Yml[i][j]);
		}
	}

	map = new double*[w];
	grad_tet = new double*[w];
	grad_phi = new double*[w];
	for(int i = 0; i < w; ++i){
		map[i] = new double[l];
		grad_tet[i] = new double[l];
		grad_phi[i] = new double[l];
		std::fill_n(map[i], l, n_ref);
		std::fill_n(grad_tet[i], l, 0);
		std::fill_n(grad_phi[i], l, 0);

		tet = i*dtet - max_tet + M_PI/2;
		std::cout << "In i = " << i+1 << " out of " << w << ", theta = " << tet << std::endl;
		for(int j = 0; j < l; ++j){
			phi = j*dphi;
			for(int ll = 1; ll <= l_max; ++ll){
				for(int mm = -ll; mm <= ll; ++mm){
					map[i][j] += Yml[ll-1][mm+ll]*Ysh(tet, phi, mm, ll);
					grad_tet[i][j] += Yml[ll-1][mm+ll]*Ysh.gtet(tet, phi, mm, ll);
					grad_phi[i][j] += Yml[ll-1][mm+ll]*Ysh.gphi(tet, phi, mm, ll);
				}
			}
			if(map[i][j] > max) max = map[i][j];
			else if(map[i][j] < min) min = map[i][j];

			if(grad_tet[i][j] > max_t) max_t = grad_tet[i][j];
			else if(grad_tet[i][j] < min_t) min_t = grad_tet[i][j];
		}
	}

	std::cout << "Max = " << max_t << ", Min = " << min_t << std::endl;
}

Turbulence_Map::~Turbulence_Map()
{
	for(int i = 0; i < w; ++i){
		delete map[i];
	}

	delete map;
	delete image;
}

void Turbulence_Map::draw(const std::string& name) const
{
	QPainter paint;
	QPen pen;
	int rgb;

	paint.begin(image);
	paint.setRenderHint(QPainter::Antialiasing, true);
    pen.setWidth(1);

	for(int i = 0; i < w; ++i){
		for(int j = 0; j < l; ++j){
            rgb = 255.0*(map[i][j] - min)/(max - min);
            pen.setColor(QColor(rgb, rgb, rgb));
            paint.setPen(pen);
            paint.drawPoint(j, i);
        }
    }
    paint.end();
	image->save(name.c_str());
}

void Turbulence_Map::write(const std::string& name) const
{
	std::ofstream fp(name + ".dat", std::ofstream::out | std::ofstream::trunc);
	std::ofstream fp_t(name + "_gtet.dat", std::ofstream::out | std::ofstream::trunc);
	std::ofstream fp_p(name + "_gphi.dat", std::ofstream::out | std::ofstream::trunc);

	fp << l << " " << w << " " << max_tet << std::endl;
	fp_t << l << " " << w << " " << max_tet << std::endl;
	fp_p << l << " " << w << " " << max_tet << std::endl;

	for(int i = 0; i < w; ++i){
		for(int j = 0; j < l; ++j){
			fp << map[i][j] << " ";
			fp_t << grad_tet[i][j] << " ";
			fp_p << grad_phi[i][j] << " ";
		}
		fp << std::endl;
		fp_t << std::endl;
		fp_p << std::endl;
	}
	fp.close();
	fp_t.close();
	fp_p.close();
}





TurbPlot::TurbPlot(const MapOpt& m) : Turbulence_Map(m)
{
	grid = new QGridLayout(this);
	plt = new QCustomPlot*[3];
	plt2 = new QCustomPlot(this);

	draw_plot2();

	for(int i = 0; i < 3; ++i){
		draw_plot(i);
		plt[i]->setFixedSize(1800, std::ceil(1375*max_tet/M_PI + 112.5));
		grid->addWidget(plt[i], i, 0);
	}
	plt2->setFixedSize(550, 300);
	grid->addWidget(plt2, 3, 0);
	this->setWindowTitle("Turbulence Map");
}

TurbPlot::~TurbPlot()
{
	delete grid;
	for(int i = 0; i < 3; ++i) delete plt[i];
	delete plt;
	delete plt2;
}

void TurbPlot::draw_plot(const int& P)
{
	plt[P] = new QCustomPlot(this);
	QCPColorMap *colorMap = new QCPColorMap(plt[P]->xAxis, plt[P]->yAxis);
	QCPColorScale *colorScale = new QCPColorScale(plt[P]);
	QCPMarginGroup *marginGroup = new QCPMarginGroup(plt[P]);


	plt[P]->axisRect()->setupFullAxesBox(true);
	plt[P]->xAxis->setLabel("phi");
	plt[P]->yAxis->setLabel("theta");
	plt[P]->xAxis->setScaleRatio(plt[P]->yAxis, 1);


	colorMap->data()->setSize(l, w);
	colorMap->data()->setRange(QCPRange(0, 2*M_PI), QCPRange(M_PI/2 - max_tet, M_PI/2 + max_tet));

	for(int i = 0; i < l; ++i){
		for(int j = 0; j < w; ++j){
			if(P == 0) colorMap->data()->setCell(i, j, 100*(map[j][i]/n_ref - 1));
			else if(P == 1) colorMap->data()->setCell(i, j, 100*grad_tet[j][i]/n_ref);
			else colorMap->data()->setCell(i, j, 100*grad_phi[j][i]/n_ref);
		}
	}

	plt[P]->plotLayout()->addElement(0, 1, colorScale);
	colorScale->setType(QCPAxis::atRight);
	colorMap->setColorScale(colorScale);
	if(!P) colorScale->axis()->setLabel("Relative refractivity in %");
	else colorScale->axis()->setLabel("Gradient value");
	colorMap->setGradient(blackGrad());
	colorMap->setInterpolate(true);
	colorMap->rescaleDataRange();
	plt[P]->axisRect()->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
	colorScale->setMarginGroup(QCP::msBottom|QCP::msTop, marginGroup);
	plt[P]->rescaleAxes();

}

void TurbPlot::draw_plot2()
{
	QVector<QCPGraphData> data(l_max), loglaw(l_max);
	QPen pen;
	QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);

	plt2->addGraph();
	pen.setColor(QColor(255,170,100));
	pen.setWidth(2);
	pen.setStyle(Qt::DotLine);
	plt2->graph(0)->setPen(pen);
	plt2->graph(0)->setName("Energy equivalent");

	plt2->addGraph();
	pen.setColor(QColor(0, 0, 0, 100));
	pen.setWidth(2);
	pen.setStyle(Qt::SolidLine);
	plt2->graph(1)->setPen(pen);
	plt2->graph(1)->setName("Equivalent power law");

	for(int i = 0; i < l_max; ++i){
		data[i].key = i+1;
		data[i].value = E[i];
		loglaw[i].key = i+1;
		loglaw[i].value = CLM*n_ref*std::pow(i+1, -alpha);
	}

	plt2->graph(0)->data()->set(data);
	plt2->graph(1)->data()->set(loglaw);
	plt2->xAxis->setLabel("l-value");
	plt2->yAxis->setLabel("Refractivity fluctuation");
	plt2->xAxis->setScaleType(QCPAxis::stLogarithmic);
	plt2->yAxis->setScaleType(QCPAxis::stLogarithmic);
	plt2->yAxis->setTicker(logTicker);
	plt2->xAxis->setTicker(logTicker);
	plt2->yAxis->setNumberFormat("eb"); // e = exponential, b = beautiful decimal powers
	plt2->yAxis->setNumberPrecision(0);
	plt2->yAxis->grid()->setSubGridVisible(true);
	plt2->xAxis->grid()->setSubGridVisible(true);
	plt2->xAxis->setRange(1, l_max);
	plt2->yAxis->setRange(1e-14, 1e-6);

	plt2->legend->setVisible(true);
	plt2->legend->setBrush(QBrush(QColor(255,255,255,150)));
	plt2->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignBottom);
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

QCPColorGradient blackGrad()
{
	QCPColorGradient grad;
	QMap<double, QColor> c;

	c[0] = Qt::black;
	c[1] = Qt::white;

	grad.setColorStops(c);
	grad.setColorInterpolation(QCPColorGradient::ciRGB);

	return grad;
}
