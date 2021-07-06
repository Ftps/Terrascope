#ifndef FLUC_HPP
#define FLUC_HPP

#include "qcustomplot.hpp"
#include "sphere_harm.hpp"

#include <random>
#include <QColor>
#include <QImage>
#include <QPainter>

#define CLM 0.0000002

struct MapOpt {
	int l, w, l_max;
	double alpha, n_ref, max_tet;
};

class Turbulence_Map : protected MapOpt {
	public:
		Turbulence_Map(const MapOpt& m, const int& n_thread);
		~Turbulence_Map();
		void draw(const std::string& name) const;
		void write(const std::string& name) const;
	protected:
		const SphereHarm Ysh;
		double **map, **grad_tet, **grad_phi;
		double min, max;
		QImage *image;
		std::vector<std::vector<double>> Yml;
		std::vector<double> E;

};

class TurbPlot : public Turbulence_Map, public QWidget {
	public:
		TurbPlot(const MapOpt& m, const int& n_thread);
		~TurbPlot();
	private:
		void draw_plot(const int& P);
		void draw_plot2();

		QGridLayout *grid;
		QCustomPlot **plt, *plt2;
};

QCPColorGradient redGrad();
QCPColorGradient blackGrad();

#endif
