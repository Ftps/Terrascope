#include "sphere_harm.hpp"

const std::string file = "../Poly/", ext = ".dat";

Poly::Poly(const int& l) : l(l)
{
	std::ifstream fp(file + 'P' + std::to_string(l) + ext, std::ifstream::in);

	if(!fp){
		std::cout << "File failed to open" << std::endl;
		exit(-1);
	}

	fp >> N;
	poly.resize(l+1);
	for(int i = 0; i <= l; ++i){
		poly[i].resize(N);
		for(int j = 0; j < N; ++j){
			fp >> poly[i][j];
		}
	}
	--N;
}

long double Poly::operator()(const int& m, const long double& x) const
{
	int i = floor(((x+1)/2.0)*N);

	return (i == N) ? poly[m][i] : poly[m][i] + (poly[m][i+1] - poly[m][i])*(0.5*N*(x+1) - i);
}


SphereHarm::SphereHarm(const int& max_l) : l_max(max_l)
{
	std::ifstream fp(file + "Coeffs" + ext, std::ifstream::in);
	int max;

	if(!fp){
		std::cout << "File failed to open" << std::endl;
		exit(-1);
	}

	fp >> max;
	if(l_max > max){
		std::cout << "Maximum l value larger then available, reducing to l_max = " << max << std::endl;
		l_max = max;
	}

	std::cout << "Loading spherical harmonics . . ." << std::endl;

	Clm.resize(l_max);
	for(int l = 0; l < l_max; ++l){
		Clm[l].resize(l+2);
		for(int m = 0; m < l+2; ++m){
			fp >> Clm[l][m];
		}
		poly.push_back(Poly(l+1));
	}
	fp.close();

	std::cout << "Harmonics loaded." << std::endl;
}

long double SphereHarm::operator()(const long double& tet, const long double& phi, const int& m, const int& l) const
{
	if(m < 0) return poly[l-1](-m, cosl(tet))*sinl(-m*phi);
	else if(m > 0) return poly[l-1](m, cosl(tet))*cosl(m*phi);
	return poly[l-1](0, cosl(tet));
}

long double SphereHarm::gtet(const long double& tet, const long double& phi, const int& m, const int& l) const
{
	if(m < 0) return sinl(-m*phi)*(((l > -m) ? poly[l-1](-m+1, cosl(tet)) : 0) - m*cotl(tet)*poly[l-1](-m, cosl(tet)));
	else if(m > 0) return cosl(m*phi)*(((l > m) ? poly[l-1](m+1, cosl(tet)) : 0) + m*cotl(tet)*poly[l-1](m, cosl(tet)));
	return poly[l-1](1, cosl(tet));
}

long double SphereHarm::gphi(const long double& tet, const long double& phi, const int& m, const int& l) const
{
	if(m < 0) return -m*poly[l-1](-m, cosl(tet))*cosl(-m*phi);
	else if(m > 0) return -m*poly[l-1](m, cosl(tet))*sinl(m*phi);
	return 0;
}




inline long double cotl(const long double& x)
{
	return tanl(M_PI_2 - x);
}
