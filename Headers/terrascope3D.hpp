#ifndef TERRASCOPE_3D_HPP
#define TERRASCOPE_3D_HPP

#include "fluctuation.hpp"

#include <chrono>
#include <thread>

struct FinalPos {
	std::array<double, 3> vel, pos;
	double intensity;
	FinalPos() : vel({0,0,0}), pos({0,0,0}), intensity(0) { }
	FinalPos(const std::array<double,3>& vel, const std::array<double,3> pos, const double& intensity)
	: vel(vel), pos(pos), intensity(intensity) { }
};

struct Substance {
	double n;
	std::array<double, 3> F;
	std::vector<double> wave;
	std::vector<std::array<double, 3>> n0;
};

struct Refractivity {
	double n;
	std::array<double, 3> gn;
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
		const double R, R2, r_max, r2_max, betay, betaz, n_ref;
		Atmosphere atm;
		double o;
		bool dev, turb;

		Planet3D(const double &R, const double& H, const double& N, const double& n_ref, const double& obf_y, const double& obf_z, const std::array<double, 3>& r, const std::string& filemap = "empty", const std::string& filename = "Config/atmosphere.dat");
		double checkDistance(const std::array<double, 3>& pos) const;
		double checkRayHit(const std::array<double, 3>& pos, const std::array<double, 3>& vel) const;
		void updateSigma(const double& l);
		void updateRefrac(const std::array<double, 3>& pos, Refractivity& n) const;

		double rot[3][3], prot[2][2];
		Turbulence_Map m;
		std::function<std::array<double,3>(std::array<double,3>)> rotate, arotate, protate, g_refrac;
		std::function<double(std::array<double,3>)> refrac;

		std::array<double, 3> cartToSphere(const std::array<double, 3>& pos, double& r_s) const;
		std::array<double, 3> sphereGrad(const std::array<double, 3>& pos_s, const Refractivity& n) const;
};

struct FlashMap {
	const int N, S;
	const double L, hh, alpha;
	int ray_hit, ray_counter;
	double Int1, Int2;
	double** map;

	FlashMap(const int& N, const int& S, const double& L, const double& hh, const double& alpha) : N(N), S(S), L(L), hh(hh), alpha(alpha), ray_hit(0), ray_counter(0), Int1(0), Int2(0), map(createMatrix(N+1)) {}
	~FlashMap() { freeMatrix(map, N+1); }
};

std::array<double,3> rayTracing(const Planet3D& p, const double& H, const double& L, const std::array<double,2>& a, double& min, const double& dt = 10);
FinalPos rayTracing2(const Planet3D& p, const double& H, const double& L, const std::array<double,2>& a, double& min, const double& dt = 10);
double rayTrace(const Planet3D& p, std::array<double,3>& pos, std::array<double,3>& vel, const double& dt = 10);

FlashMap* mapGen(const Planet3D& p, const int& N, const double& S, const double& L, const double& hh);
FlashMap* mapGen2(const Planet3D& p, const int& N, const double& S, const double& L, const double& hh);
FlashMap* mapThread(const Planet3D& p, const int& N, const double& S, const double& L, const double& hh, const int& thread_n);
void map_threaded(FlashMap* m, const Planet3D& p, const double& opt, const double& dr, const double& phi1, const double& phi2, const double& dphi);

#endif
