#include <QApplication>
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

#include <chrono>

#include "qcustomplot.hpp"
#include "terrascope3D.hpp"

#define ARCSEC 206264.806247
#define OP 0.05

class ImageGen : public QWidget {
	public:
		ImageGen(const Planet3D& p, const double& L, const int& N = 200);

	private:
		const Planet3D p;
		const double L;
		const int N;

		QGridLayout *grid;
		QCustomPlot *plt;

		void drawCFM(); // CFM - Central Flash Map
		//std::array<double,2> angleToCoord(const std::array<double,2>& a, const double& max);
};

QCPColorGradient redGrad();