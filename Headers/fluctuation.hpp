#ifndef FLUC_HPP
#define FLUC_HPP

#include <fstream>

#include "sphere_harm.hpp"

#define MAP 0
#define G_T 1
#define G_P 2

class Turbulence_Map {
	public:
		Turbulence_Map(const std::string& filename);
		void destroy();
		double n(const double& tet, const double& phi) const;
		double n_tet(const double& tet, const double& phi) const;
		double n_phi(const double& tet, const double& phi) const;
	private:
		int l, w;
		double max_tet, dphi, dtet, h2;
		double **map[3];

		double biLin(const double& tet, const double& phi, const int& m) const;
};



#endif
