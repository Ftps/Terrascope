#include "fluctuation.hpp"
#include <cassert>

Turbulence_Map::Turbulence_Map(const std::string& filename)
{
	if(!filename.compare("empty")) return;

	std::ifstream fp(filename + ".dat", std::ifstream::in);
	std::ifstream fp_t(filename + "_gtet.dat", std::ifstream::in);
	std::ifstream fp_p(filename + "_gphi.dat", std::ifstream::in);

	fp >> w;
	fp >> l;
	fp >> max_tet;

	while(fp_t.get() != '\n');
	while(fp_p.get() != '\n');

	map[MAP] = new double*[w];
	map[G_T] = new double*[w];
	map[G_P] = new double*[w];
	for(int i = 0; i < w; ++i){
		map[MAP][i] = new double[l];
		map[G_T][i] = new double[l];
		map[G_P][i] = new double[l];
	}

	for(int i = 0; i < l; ++i){
		for(int j = 0; j < w; ++j){
			fp >> map[MAP][j][i];
			fp_t >> map[G_T][j][i];
			fp_p >> map[G_P][j][i];
		}
	}

	dphi = (2*M_PI)/((double)(w-1));
	dtet = (2*max_tet)/((double)(l-1));
	h2 = 1/(dphi*dtet);
}

void Turbulence_Map::destroy()
{
	for(int i = 0; i < w; ++i){
		delete map[MAP][i];
		delete map[G_T][i];
		delete map[G_P][i];
	}

	delete map[MAP];
	delete map[G_T];
	delete map[G_P];
}

double Turbulence_Map::n(const double& tet, const double& phi) const
{
	return biLin(tet, phi, MAP);
}

double Turbulence_Map::n_tet(const double& tet, const double& phi) const
{
	return biLin(tet, phi, G_T);
}

double Turbulence_Map::n_phi(const double& tet, const double& phi) const
{
	return biLin(tet, phi, G_P);
}




double Turbulence_Map::biLin(const double& tet, const double& phi, const int& m) const
{
	int i, j;
	double x, y;

	x = (phi/(2*M_PI))*(w-2);
	y = ((tet - M_PI_2 + max_tet)/(2*max_tet))*(l-2);
	i = std::floor(x);
	j = std::floor(y);

	//x = phi - i*dphi;
	//y = tet - M_PI_2 + max_tet - j*dtet;
	x = x - i;
	y = y - j;

	//assert(0 <= i && i < w-1);
	//assert(0 <= j && j < l-1);

	return (1-x)*(1-y)*map[m][i][j] + x*(1-y)*map[m][i+1][j] + (1-x)*y*map[m][i][j+1] + x*y*map[m][i+1][j+1];
	//return h2*((dphi-x)*(dtet-y)*map[m][i][j] + x*(dtet-y)*map[m][i+1][j] + (dphi-x)*y*map[m][i][j+1] + x*y*map[m][i+1][j+1]);
}
