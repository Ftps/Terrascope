#ifndef TERRASCOPE_3D_HPP
#define TERRASCOPE_3D_HPP

#include "auxiliar.hpp"

#define X 0
#define Y 1
#define TET 0
#define PHI 1


struct Planet3D {
	double R, r_max, obf;
	std::function<dddd> n;

	Planet3D(const double &R, const double& r_max, const double& obf, const std::function<dddd>& n);
};


std::array<double,2> rayTracing(const Planet3D& p, const double& L, const std::array<double,2>& a, const double& dz = 10);

#endif
