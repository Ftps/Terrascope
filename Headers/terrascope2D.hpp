#ifndef TERRASCOPE_2D_HPP
#define TERRASCOPE_2D_HPP

#include <functional>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#define ddd double(double, double)

#define LOG {std::cout << "IN LINE " << __LINE__  << " OF FILE " << __FILE__ << std::endl; fflush(stdout);}
#define Print(X) std::cout << X << std::endl

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

double fisqrt(double n);
double sq(double a);

double gx2D(const std::function<ddd>& f, const double& x, const double& z, const double& h = 0.001);
double gy2D(const std::function<ddd>& f, const double& x, const double& z, const double& h = 0.001);

double ray_tracer2D(const std::function<ddd>& n, const double& R, const double& r_max, const double& obf, const double& L, const double& a_init, const double& dz = 10);
double ray_tracer2D_hor(const std::function<ddd>& n, const double& R, const double& r_max, const double& obf, const double& y, const double& dx = 10);

double focalPoint(const Planet2D& p, const int& N);

#endif
