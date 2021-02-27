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
#include "terrascope2D.hpp"

class DrawRay : public QWidget {
	public:
		DrawRay(const Planet2D& p, const double& L, const double& a_entry);
	private:
		QGridLayout *grid;
		QCustomPlot *customPlot;
}
