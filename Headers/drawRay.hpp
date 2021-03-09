#ifndef DRAW_RAY_HPP
#define DRAW_RAY_HPP

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

#include <array>

#include "qcustomplot.hpp"
#include "terrascope2D.hpp"

#define CIRC_NUM 500

class DrawRay : public QWidget {
	public:
		DrawRay(const Planet2D& p, const double& L, const double& a_entry, const int& type);
	private:
		QGridLayout *grid;
		QCustomPlot *plt;
};

void drawCircle(QCustomPlot* plt, const double& r, const QColor& line_c, const QColor& fill_c);
int drawRay(QCustomPlot* plt, const Planet2D& p, const double& L, const double& a_entry, const QColor& line_c);
int drawRayHor(QCustomPlot* plt, const Planet2D& p, const double& Y, const QColor& line_c);
void drawRayDensity(QCustomPlot *plt, const Planet2D& p, const QColor& line_c);

#endif
