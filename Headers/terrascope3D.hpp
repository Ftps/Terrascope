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
#include <chrono>

#include "qcustomplot.hpp"
#include "terrascope2D.hpp"

#define dddd double(double, double, double)
#define X 0
#define Y 1
#define TET 0
#define PHI 1
#define ARCSEC 206264.806247

struct Planet3D {
	double R, r_max, obf;
	std::function<dddd> n;

	Planet3D(const double &R, const double& r_max, const double& obf, const std::function<dddd>& n);
};

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

template<typename T>
T** createMatrix(int N)
{
	T** m = new T*[N];

	for(int i = 0; i < N; ++i){
		m[i] = new T[N]();
	}

	return m;
}

template<typename T>
void freeMatrix(T** m, int N)
{
	for(int i = 0; i < N; ++i){
		delete m[i];
	}

	delete m;
}


double gx3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);
double gy3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);
double gz3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);

std::array<double,2> rayTracing(const Planet3D& p, const double& L, const std::array<double,2>& a, const double& dz = 10);



template<typename T>
inline std::array<T,2> operator*(const double& a, std::array<T,2> v)
{
	v[0] = a*v[0];
	v[1] = a*v[1];

	return v;
}

template<typename T>
inline std::array<T,2> operator+(std::array<T,2> v, const double& a)
{
	v[0] += a;
	v[1] += a;

	return v;
}

template<typename T>
inline std::array<T,2> operator/(std::array<T,2> v, const double& a)
{
	v[0] = v[0]/a;
	v[1] = v[1]/a;

	return v;
}

#endif
