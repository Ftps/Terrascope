#include "terrascope2D.hpp"
/*#include <QApplication>
#include <QWidget>
#include <QPushButton>
#include <QLabel>
#include <QGridLayout>
#include <QSignalMapper>
#include <QLineEdit>
#include <QComboBox>
#include <QPixmap>
#include <QStringList>
#include <QProgressBar>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QColorDialog>
#include <QHeaderView>
#include <QBrush>

//#include "qcustomplot.hpp"

class PlotGraph : public QWidget {
	public:
		PlotGraph(QWidget *parent = 0) : QWidget(parent)
		{
			grid = new QGridLayout(this);
			customPlot = new QCustomPlot(this);

			// add two new graphs and set their look:
			customPlot->addGraph();
			customPlot->graph(0)->setPen(QPen(Qt::blue)); // line color blue for first graph
			customPlot->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20))); // first graph will be filled with translucent blue
			customPlot->addGraph();
			customPlot->graph(1)->setPen(QPen(Qt::red)); // line color red for second graph
			customPlot->graph(1)->setBrush(QBrush(QColor(255, 0, 0, 20)));
			// generate some points of data (y0 for first, y1 for second graph):
			QVector<double> x(251), y0(251), y1(251);
			for (int i=0; i<251; ++i)
			{
			  x[i] = i;
			  y0[i] = qExp(-i/150.0)*qCos(i/10.0); // exponentially decaying cosine
			  y1[i] = qExp(-i/150.0);              // exponential envelope
			}
			// configure right and top axis to show ticks but no labels:
			// (see QCPAxisRect::setupFullAxesBox for a quicker method to do this)
			customPlot->xAxis2->setVisible(true);
			customPlot->xAxis2->setTickLabels(false);
			customPlot->yAxis2->setVisible(true);
			customPlot->yAxis2->setTickLabels(false);
			// make left and bottom axes always transfer their ranges to right and top axes:
			connect(customPlot->xAxis, SIGNAL(rangeChanged(QCPRange)), customPlot->xAxis2, SLOT(setRange(QCPRange)));
			connect(customPlot->yAxis, SIGNAL(rangeChanged(QCPRange)), customPlot->yAxis2, SLOT(setRange(QCPRange)));
			// pass data points to graphs:
			customPlot->graph(0)->setData(x, y0);
			customPlot->graph(1)->setData(x, y1);
			// let the ranges scale themselves so graph 0 fits perfectly in the visible area:
			customPlot->graph(0)->rescaleAxes();
			// same thing for graph 1, but only enlarge ranges (in case graph 1 is smaller than graph 0):
			customPlot->graph(1)->rescaleAxes(true);
			// Note: we could have also just called customPlot->rescaleAxes(); instead
			// Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
			customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

			grid->addWidget(customPlot, 0, 0, 10, 5);
			resize(QDesktopWidget().availableGeometry(this).size() * 0.7);
		}
	private:
		QGridLayout *grid;
		QCustomPlot *customPlot;
};

int main(int argc, char* argv[])
{
	QApplication a(argc, argv);
	PlotGraph w;
	w.show();

	return a.exec();
}*/


int main(int argc, char* argv[])
{
	double R = 2574;			// planet radius
	double H = 25;				// atmospheric scale height
	double r_max = R + 10*H;	// top layer of the atmosphere
	//double N = 0.000293;		// surface refractivity
	double L = 500000;
	double a, a_init = 0.5*(R + r_max)/L;
	int N = 350000, iter, tot = 0;
	double h = (1.05-0.96)*a_init/(double)N;

	std::function<ddd> n = [](double x, double y){ return 1 + 0.1*exp(-(sqrt(x*x + y*y)-2574)/25); };

	for(int i = 0; i <= N; ++i){
		a = 1.05*a_init - i*h;
		if(!(i % 1000)) Print(i);
		ray_tracer2D(n, R, r_max, L, a, iter, 10);
		tot += iter;
	}


	Print("Total rays: " << N);
	Print("Total Iterations: " << tot);
	Print("Average Iterations per ray: " << tot/(double)N);

	Print("\nSingle ray test:");
	std::cout << "Entry ray angle: " << 180*a_init/3.14159 << " º" << std::endl;
	std::cout << "Exit ray angle: " << 180*ray_tracer2D(n, R, r_max, L, a_init, iter)/3.14159 << " º" << std::endl;

	return 0;
}
