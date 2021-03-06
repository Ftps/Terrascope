#ifndef AUXILIAR_HPP
#define AUXILIAR_HPP

#include <functional>
#include <iostream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <array>

#define ddd double(double, double)
#define dddd double(double, double, double)
#define addd std::array<double,3>(double, double, double)

#define da double(std::array<double,3>)
#define aa std::array<double,3>(std::array<double,3>)

#define X 0
#define Y 1
#define Z 2
#define BRT 2
#define TET 0
#define PHI 1
#define PSI 2

#define LOG {std::cout << "IN LINE " << __LINE__  << " OF FILE " << __FILE__ << std::endl; fflush(stdout);}
#define Print(X) std::cout << X << std::endl
#define ST(X) std::to_string(X)

struct FDIV{
	std::function<da> f;
	std::function<aa> df;
};

double fisqrt(double n);
double sq(double a);

double gx2D(const std::function<ddd>& f, const double& x, const double& z, const double& h = 0.001);
double gy2D(const std::function<ddd>& f, const double& x, const double& z, const double& h = 0.001);

double gx3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);
double gy3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);
double gz3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);

std::function<dddd> generateRefIndex(const double& n_surf, const double& H, const double& R, const double& obf_y, const double& obf_z);
FDIV generateRefFuncs(const double& n_surf, const double& H, const double& R, const double& obf_y, const double& obf_z);

std::array<double, 2> planeIntersect(const std::array<double, 3>& pos, const std::array<double, 3>& vel, const double& L, const double& tet, const double& h);
std::array<double, 2> linReg(const std::vector<std::pair<double, double>>& xy);

template<typename T>
inline std::array<T,3> operator*(const double& a, const std::array<T,3>& v)
{
	std::array<double, 3> res;

	res[0] = a*v[0];
	res[1] = a*v[1];
	res[2] = a*v[2];

	return res;
}

template<typename T>
inline std::array<T,3> operator+(const std::array<T,3> v, const std::array<double,3> w)
{
	std::array<double, 3> res;

	res[0] = v[0] + w[0];
	res[1] = v[1] + w[1];
	res[2] = v[2] + w[2];

	return res;
}

template<typename T>
inline std::array<T,3> operator-(const std::array<T,3> v, const std::array<double,3> w)
{
	std::array<double, 3> res;

	res[0] = v[0] - w[0];
	res[1] = v[1] - w[1];
	res[2] = v[2] - w[2];

	return res;
}

template<typename T>
inline void operator+=(std::array<T,3>& v, const std::array<double,3>& w)
{
	v[0] += w[0];
	v[1] += w[1];
	v[2] += w[2];
}

template<typename T>
inline double operator*(const std::array<T,3>& v, const std::array<double,3>& w)
{
	return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

inline std::array<double,3> operator*(const double A[3][3], const std::array<double,3>& b)
{
	std::array<double,3> res = {0,0,0};

	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			res[i] += A[i][j]*b[j];
		}
	}

	return res;
}

double** createMatrix(const int& N);
void freeMatrix(double const* const* m, const int N);

std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& v);
std::ostream& operator<<(std::ostream& os, const double A[3][3]);

#endif
