#ifndef SPHERE_HARM_HPP
#define SPEHRE_HARM_HPP

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <gmp.h>

const std::string fileloc = "../Poly/";

long double factorial(const int& a);
long double coeffFact(const int& l, const int& m);
long double binTop(long double a, const int& k);
long double ALpoly(const long double& x, const int& m, const int& l, long double& prev_LP);
long double sphereHarm(const long double& tet, const long double& phi, const int& m, const int& l);
long double sphereHarmIm(const long double& tet, const long double& phi, const int& m, const int& l);

void fillNumbers(const int& l_max, const int& N);
void verifyCoeff();

#endif
