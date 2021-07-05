#ifndef SPHERE_HARM_HPP
#define SPEHRE_HARM_HPP

#include "auxiliar.hpp"

#include <tgmath.h>
#include <string>
#include <fstream>

class Poly {
	public:
		Poly(const int& l);
		long double operator()(const int& m, const long double& x) const;
	private:
		int N, l;
		std::vector<std::vector<long double>> poly;
};

class SphereHarm {
	public:
		SphereHarm(const int& max_l);
		long double operator()(const long double& tet, const long double& phi, const int& m, const int& l) const;
		long double gtet(const long double& tet, const long double& phi, const int& m, const int& l) const;
		long double gphi(const long double& tet, const long double& phi, const int& m, const int& l) const;
	private:
		int l_max;
		std::vector<Poly> poly;
		std::vector<std::vector<long double>> Clm;

};

inline long double cotl(const long double& x);

#endif
