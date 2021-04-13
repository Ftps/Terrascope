#ifndef TERRASCOPE_3D_HPP
#define TERRASCOPE_3D_HPP

#include "auxiliar.hpp"
#include <fstream>

struct Substance {
	double n;
	std::array<double, 3> F;
	std::vector<double> wave;
	std::vector<std::array<double, 3>> n0;
};

class Atmosphere {
	public:
		Atmosphere(const double& N, const double& H, const std::string& filename = "Config/atmosphere.dat");
		double getSigma(const double& l) const;

		const double N, H;
		std::vector<Substance> S;
};

class Planet3D {
	public:
		const double R, R2, r_max, r2_max, betay, betaz;
		Atmosphere atm;
		double o;
		FDIV n;
		double rot[3][3], prot[2][2];
		bool dev;
		std::function<std::array<double,3>(std::array<double,3>)> rotate, arotate, protate;

		Planet3D(const double &R, const double& H, const double& N, const double& obf_y, const double& obf_z, const std::array<double, 3>& r, const FDIV& n, const std::string& filename = "Config/atmosphere.dat");
		double checkDistance(const std::array<double, 3>& pos) const;
		double checkRayHit(const std::array<double, 3>& pos, const std::array<double, 3>& vel) const;
		void updateSigma(const double& l);
};


std::array<double,3> rayTracing(const Planet3D& p, const double& H, const double& L, const std::array<double,2>& a, const double& dt = 10);

#endif
