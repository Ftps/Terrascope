#ifndef TERRASCOPE_3D_HPP
#define TERRASCOPE_3D_HPP

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

#define dddd double(double, double, double)
#define X 0
#define Y 1
#define TET 0
#define PHI 1

struct Planet3D {
	double R, r_max, obf;
	std::function<dddd> n;

	Planet3D(const double &R, const double& r_max, const double& obf, const std::function<dddd>& n);
};

class ImageGen : Planet3D {
	public:

	private:

};

double gx3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);
double gy3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);
double gz3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);

std::array<double,2> rayTracing(const Planet3D& p, const double& L, const std::array<double,2>& a, const double& dz = 10);

#endif
