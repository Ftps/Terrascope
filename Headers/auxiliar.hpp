#ifndef AUXILIAR_HPP
#define AUXILIAR_HPP

#include <functional>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <array>

#define ddd double(double, double)
#define dddd double(double, double, double)

#define LOG {std::cout << "IN LINE " << __LINE__  << " OF FILE " << __FILE__ << std::endl; fflush(stdout);}
#define Print(X) std::cout << X << std::endl

double fisqrt(double n);
double sq(double a);

double gx2D(const std::function<ddd>& f, const double& x, const double& z, const double& h = 0.001);
double gy2D(const std::function<ddd>& f, const double& x, const double& z, const double& h = 0.001);

double gx3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);
double gy3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);
double gz3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h = 0.001);

template<typename T>
inline std::array<T,2> operator*(const double& a, std::array<T,2> v)
{
	v[0] = a*v[0];
	v[1] = a*v[1];

	return v;
}

template<typename T>
inline std::array<T,2> operator+(std::array<T,2> v, const double& a)
{
	v[0] += a;
	v[1] += a;

	return v;
}

template<typename T>
inline std::array<T,2> operator/(std::array<T,2> v, const double& a)
{
	v[0] = v[0]/a;
	v[1] = v[1]/a;

	return v;
}

template<typename T>
T** createMatrix(int N)
{
	T** m = new T*[N];

	for(int i = 0; i < N; ++i){
		m[i] = new T[N]();
	}

	return m;
}

template<typename T>
void freeMatrix(T** m, int N)
{
	for(int i = 0; i < N; ++i){
		delete m[i];
	}

	delete m;
}

#endif
