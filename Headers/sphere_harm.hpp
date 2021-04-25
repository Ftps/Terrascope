#ifndef SPHERE_HARM_HPP
#define SPEHRE_HARM_HPP

#include <cmath>
#include <tgmath.h>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>

#define LOG {std::cout << "IN LINE " << __LINE__ << " OF FILE " << __FILE__ << std::endl; fflush(stdout);}

class Poly {
	public:
		Poly(const int& l);
		long double operator()(const int& m, const long double& x);
	//private:
		int N, l;
		std::vector<std::vector<long double>> poly;
};

class SphereHarm {
	public:
		SphereHarm(const int& max_l);
		long double operator()(const long double& tet, const long double& phi, const int& m, const int& l);
	//private:
		int l_max;
		std::vector<Poly> poly;
		std::vector<std::vector<long double>> Clm;

};


#endif
