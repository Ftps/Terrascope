#ifndef TERRASCOPE_2D_HPP
#define TERRASCOPE_2D_HPP

#include "auxiliar.hpp"

struct Planet2D {
	double R, r_max, obf;
	std::function<ddd> n;

	Planet2D(const double& R, const double& r_max, const double& obf, const std::function<ddd>& n);
};

class Ray2D : Planet2D {
	public:
		std::vector<double> x, y;
		double a_entry, a_exit, L;
		int i;

		Ray2D(const Planet2D& p);
		int ray_tracer(const double& Ll, const double& Aa_entry, const double& dx = 1);
		int ray_tracer_hor(const double& Y, const double& dx = 1);
};




double ray_tracer2D(const std::function<ddd>& n, const double& R, const double& r_max, const double& obf, const double& L, const double& a_init, const double& dz = 10);
double ray_tracer2D_hor(const std::function<ddd>& n, const double& R, const double& r_max, const double& obf, const double& y, const double& dx = 10);

double focalPoint(const Planet2D& p, const int& N);

#endif
