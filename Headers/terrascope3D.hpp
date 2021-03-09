#ifndef TERRASCOPE_3D_HPP
#define TERRASCOPE_3D_HPP

#include "terrascope2D.hpp"

#define dddd double(double, double, double)

struct Planet3D {
	double R, r_max, obf;
	std::function<dddd> n;

	Planet3D(const double &R, const double& r_max, const double& obf, const std::function<dddd>& n);
};

class ImageGen {
	public:

	private:

};

#endif
