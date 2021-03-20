#ifndef TERRASCOPE_3D_HPP
#define TERRASCOPE_3D_HPP

#include "auxiliar.hpp"

#define X 0
#define Y 1
#define Z 2
#define TET 0
#define PHI 1


class Planet3D {
	public:
		const double R, R2, r_max, r2_max, betax, betay;
		double rot[3][3];
		std::function<dddd> n;
		bool dev;
		std::function<std::array<double,3>(std::array<double,3>)> rotate;

		Planet3D(const double &R, const double& r_max, const double& obf_x, const double& obf_y, const std::array<double, 2> r, const std::function<dddd>& n);
		double checkDistance(double x, double y, double z) const;
		double checkRayHit(const double& L, const double& vx, const double& vy, const double& vz) const;
};


std::array<double,2> rayTracing(const Planet3D& p, const double& L, const std::array<double,2>& a, const double& dz = 10);

#endif
