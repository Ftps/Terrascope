#include "terrascope3D.hpp"

Atmosphere::Atmosphere(const double& N, const double& H, const std::string& filename) : N(N), H(H)
{
	std::ifstream fp(filename);
	std::function<double(double)> ff;
	int n_mol, n_l, aux;
	std::vector<int> wave(0);
	std::vector<std::array<double, 3>> factors;
	std::array<double, 3> aux_N;
	Substance s;

	fp >> n_mol;

	for(int i = 0; i < n_mol; ++i){
		fp >> s.n;
		fp >> n_l;

		for(int j = 0; j < n_l; ++j){
			fp >> aux;
			s.wave.push_back(aux*1e-9);

			for(int k = 0; k < 3; ++k) fp >> aux_N[k];
			s.n0.push_back(aux_N);
		}

		for(int k = 0; k < 3; ++k) fp >> aux_N[k];

		s.F = aux_N;

		S.push_back(s);

		s.wave.clear();
		s.n0.clear();
	}

}

double Atmosphere::getSigma(const double& l) const
{
	double o = 2*sq(sq(2*M_PI/l))/(3*M_PI);
	double fN = 0, Fk, n;

	for(Substance s : S){
		Fk = s.F[0] + s.F[1]/sq(l*1e6) + s.F[2]/sq(sq(l*1e6));
		n = 0;
		for(int i = 0; i < (int)s.wave.size(); ++i){
			if(l < s.wave[i] || i == (int)s.wave.size()-1){
				n = s.n0[i][0] + s.n0[i][1]/(s.n0[i][2] - 1/sq(l*1e6));
			}
		}
		Print("n = " << s.n);
		Print("Fk = " << Fk);
		Print("n0 = " << n);
		fN += s.n*sq(n)*Fk;
		//fN += atm_frac[i]*sq(n[i](l))*F[i](l);
	}

	Print("fN = " << fN);
	Print("o = " << o);
	Print("N = " << N);

	return o*fN/N;
}



Planet3D::Planet3D(const double &R, const double& H, const double& N, const double& obf_y, const double& obf_z, const std::array<double, 3>& r, const FDIV& n_ref, const std::string& filename) : R(R), R2(sq(R)), r_max(R + 15*H), r2_max(sq(r_max)), betay(1/(1+obf_y)), betaz(1/(1+obf_z)), atm(N, H, filename), n(n_ref)
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

void Planet3D::updateSigma(const double& l)
{
	o = atm.getSigma(l);
}


std::array<double,3> rayTracing(const Planet3D& p, const double& H, const double& L, const std::array<double,2>& a, const double& dt)
{
	double n2, n1, r, k, tau = 0, st = sin(a[TET]);
	std::array<double, 3> vel, pos, gn, dv;

	vel = p.rotate({st*cos(a[PHI]), st*sin(a[PHI]), cos(a[TET])});	//
	pos = p.rotate({H, 0, -L});										// rotate velocity and position to planet's reference frame

	if(!(k = p.checkRayHit(pos, vel))){				// check if the ray enters the planet
		vel = p.protate(vel);
		return {atan(vel[X]/vel[Z]), atan(vel[Y]/vel[Z]), 1}; // if not, return angles to the original reference frame
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

		if(!(r = p.checkDistance(pos))) return {-100, -100, 0};	// check if inside the planet, if not return inside planet value {-100, -100}

		tau += n1*exp(-(sqrt(r) - p.R)/p.atm.H);

	}while(r < p.r2_max);	// check if outside the atmosphere

	vel = p.protate(p.arotate(vel));	// rotate velocity back to original reference frame

	return {atan(vel[X]/vel[Z]), atan(vel[Y]/vel[Z]), exp(-p.o*tau*dt)};	// return exit angles
}
