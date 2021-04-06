#ifndef TERRASCOPE_3D_HPP
#define TERRASCOPE_3D_HPP

#include "auxiliar.hpp"

class Planet3D {
	public:
		const double R, R2, r_max, r2_max, betay, betaz;
		FDIV n;
		double rot[3][3], prot[2][2];
		bool dev;
		std::function<std::array<double,3>(std::array<double,3>)> rotate, arotate, protate;

		Planet3D(const double &R, const double& r_max, const double& obf_y, const double& obf_z, const std::array<double, 3>& r, const FDIV& n);
		double checkDistance(const std::array<double, 3>& pos) const;
		double checkRayHit(const std::array<double, 3>& pos, const std::array<double, 3>& vel) const;
};


std::array<double,3> rayTracing(const Planet3D& p, const double& H, const double& L, const std::array<double,2>& a, const double& dt = 10);

#endif
