#include "terrascope3D.hpp"

Planet3D::Planet3D(const double &R, const double& r_max, const double& obf_y, const double& obf_z, const std::array<double, 3>& r, const FDIV& n_ref) : R(R), R2(sq(R)), r_max(r_max), r2_max(sq(r_max)), betay(1/(1+obf_y)), betaz(1/(1+obf_z)), n(n_ref)
{
	if(r[TET] || r[PHI]){
		dev = true;
		rot[X][X] = cos(r[TET]);
		rot[X][Y] = 0;
		rot[X][Z] = sin(r[TET]);

		rot[Y][X] = sin(r[PHI])*sin(r[TET]);
		rot[Y][Y] = cos(r[PHI]);
		rot[Y][Z] = -sin(r[PHI])*cos(r[TET]);

		rot[Z][X] = -sin(r[TET])*cos(r[PHI]);
		rot[Z][Y] = sin(r[PHI]);
		rot[Z][Z] = cos(r[TET])*cos(r[PHI]);

		rotate = [this](const std::array<double, 3>& pos_init){
			std::array<double, 3> pos = {0, 0, 0};
			for(int i = 0; i < 3; ++i){
				for(int j = 0; j < 3; ++j){
					pos[i] += rot[i][j]*pos_init[j];
				}
			}
			return pos;
		};
		arotate = [this](const std::array<double, 3>& pos_init){
			std::array<double, 3> pos = {0, 0, 0};
			for(int i = 0; i < 3; ++i){
				for(int j = 0; j < 3; ++j){
					pos[i] += rot[j][i]*pos_init[j];
				}
			}
			return pos;
		};
	}
	else{
		dev = false;
		rotate = [](const std::array<double, 3>& pos_init){ return pos_init; };
		arotate = [](const std::array<double, 3>& pos_init){ return pos_init; };
	}

	if(r[PSI]){
		prot[X][X] = cos(r[PSI]);
		prot[X][Y] = -sin(r[PSI]);

		prot[Y][X] = sin(r[PSI]);
		prot[Y][Y] = cos(r[PSI]);

		protate = [this](const std::array<double, 3>& pos_init){
			std::array<double, 3> pos = {0, 0, pos_init[Z]};
			pos[X] = prot[X][X]*pos_init[X] + prot[X][Y]*pos_init[Y];
			pos[Y] = prot[Y][X]*pos_init[X] + prot[Y][Y]*pos_init[Y];

			return pos;
		};
	}
	else{
		protate = [](const std::array<double, 3>& pos_init){ return pos_init; };
	}


}

double Planet3D::checkDistance(const std::array<double,3>& pos) const
{
	double r = sq(pos[X]) + sq(betay*pos[Y]) + sq(betaz*pos[Z]);

	if(r < R2) return 0;
	else return r;
}

double Planet3D::checkRayHit(const std::array<double, 3>& pos, const std::array<double, 3>& vel) const
{
	double a = sq(vel[X]) + sq(betay*vel[Y]) + sq(betaz*vel[Z]);
	double b = pos[X]*vel[X] + sq(betay)*pos[Y]*vel[Y] + sq(betaz)*pos[Z]*vel[Z];
	double c = sq(pos[X]) + sq(betay*pos[Y]) + sq(betaz*pos[Z]) - r2_max;
	double det = sq(b) - a*c;

	if(det < 0) return 0;

	return (-b - sqrt(det))/a;

}


std::array<double,2> rayTracing(const Planet3D& p, const double& H, const double& L, const std::array<double,2>& a, const double& dt)
{
	double n2, n1, r, k, st = sin(a[TET]);
	std::array<double, 3> vel, pos, gn, dv;

	vel = p.rotate({st*cos(a[PHI]), st*sin(a[PHI]), cos(a[TET])});	//
	pos = p.rotate({H, 0, -L});										// rotate velocity and position to planet's reference frame

	if(!(k = p.checkRayHit(pos, vel))){				// check if the ray enters the planet
		vel = p.protate(vel);
		return {atan(vel[X]/vel[Z]), atan(vel[Y]/vel[Z])}; // if not, return angles to the original reference frame
	}

	pos = k*vel + pos;		// update position to the atmosphere's entry point

	do{
		n1 = 1/p.n.f(pos);		//
		n2 = sq(n1);			// necessary constants
		gn = p.n.df(pos);	// gradient of the refraction index

		k = gn*vel;			// inner product between the gradient and the veloicty unit vector
		dv = n2*(gn - k*vel);	// derivative of the velocity unit vector

		vel += dt*dv;		// update velocity
		pos += (n1*dt)*vel;	// update position

		if(!(r = p.checkDistance(pos))) return {-100, -100};	// check if inside the planet, if not return inside planet value {-100, -100}

	}while(r < p.r2_max);	// check if outside the atmosphere

	vel = p.protate(p.arotate(vel));	// rotate velocity back to original reference frame

	return {atan(vel[X]/vel[Z]), atan(vel[Y]/vel[Z])};	// return exit angles
}
