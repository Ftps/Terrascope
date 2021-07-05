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

#include "qcustomplot.hpp"
#include "terrascope3D.hpp"

#define ARCSEC 206264.806247

class ImageGen : public QWidget {
	public:
		ImageGen(Planet3D& p, const double& l, const double& L, const double& S, const int& N = 200, const double& h = 0);

	private:
		const double L, S, l, h;
		const int N;

		QGridLayout *grid;
		QCustomPlot *plt;

		void drawCFM(Planet3D& p); // CFM - Central Flash Map
		//std::array<double,2> angleToCoord(const std::array<double,2>& a, const double& max);
};

QCPColorGradient redGrad();
