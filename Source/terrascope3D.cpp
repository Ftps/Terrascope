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
				break;
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



Planet3D::Planet3D(const double &R, const double& H, const double& N, const double& n_ref, const double& obf_y, const double& obf_z, const std::array<double, 3>& r, const std::string& filemap, const std::string& filename) : R(R), R2(sq(R)), r_max(R + 15*H), r2_max(sq(r_max)), betay(1/(1+obf_y)), betaz(1/(1+obf_z)), atm(N, H, filename), n_ref(n_ref), m(filemap)
{
	int count = 0;

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

	turb = filemap.compare("empty") && !n_ref;
	if(turb){
		for(double phi = 0; phi < 2*M_PI; phi += 2*M_PI/100.0){
			++count;
			this->n_ref += m.n(M_PI/2, phi);
		}
		this->n_ref = this->n_ref/count;
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
	long double a = sq(vel[X]) + sq(betay*vel[Y]) + sq(betaz*vel[Z]);
	long double b = pos[X]*vel[X] + sq(betay)*pos[Y]*vel[Y] + sq(betaz)*pos[Z]*vel[Z];
	long double c = sq(pos[X]) + sq(betay*pos[Y]) + sq(betaz*pos[Z]) - r2_max;
	long double det = sq(b) - a*c;

	if(det < 0) return 0;

	return (-b - sqrt(det))/a;
}

void Planet3D::updateSigma(const double& l)
{
	o = atm.getSigma(l);
}

void Planet3D::updateRefrac(const std::array<double, 3>& pos, Refractivity& n) const
{
	if(turb){
		double r_s;
		std::array<double, 3> pos_s = cartToSphere(pos, r_s);
		n.n = m.n(pos_s[Y], pos_s[Z]);
		std::array<double, 3> grad_s = sphereGrad(pos_s, n);
		double Jac[3][3], ctet = 1.0/(sq(pos_s[X])*r_s), cphi = betay/sq(r_s);

		Jac[X][X] = pos[X]/pos_s[X];
		Jac[Y][X] = sq(betay)*pos[Y]/pos_s[X];
		Jac[Z][X] = sq(betaz)*pos[Z]/pos_s[X];

		Jac[X][Y] = ctet*betaz*pos[X]*pos[Z];
		Jac[Y][Y] = ctet*betay*betaz*pos[Y]*pos[Z];
		Jac[Z][Y] = ctet*betaz*sq(r_s);

		Jac[X][Z] = -cphi*pos[Y];
		Jac[Y][Z] = cphi*pos[X];
		Jac[Z][Z] = 0;

		n.gn = Jac*grad_s;
		n.n = 1 + n.n*exp(-(pos_s[X] - R)/atm.H);

		/*Print("\n\n");
		Print(pos);
		Print(pos_s);
		Print(Jac);
		Print(grad_s);
		Print(n.gn);

		std::array<double, 3> pos_r = {pos[X], sq(betay)*pos[Y], sq(betaz)*pos[Z]};
		double ex = exp(-(pos_s[X] - R)/atm.H);
		Print((-ex*n.n/(atm.H*pos_s[X]))*pos_r);
		LOG
		getchar();*/
	}
	else{
		double r = sqrt(sq(pos[X]) + sq(betay*pos[Y]) + sq(betaz*pos[Z]));
		double ex = exp(-(r - R)/atm.H);
		std::array<double, 3> pos_r = {pos[X], sq(betay)*pos[Y], sq(betaz)*pos[Z]};

		n.n = 1 + n_ref*ex;
		n.gn = (-ex*n_ref/(atm.H*r))*pos_r;
	}
}

std::array<double, 3> Planet3D::cartToSphere(const std::array<double, 3>& pos, double& r_s) const
{
	std::array<double, 3> sphere;

	r_s = sq(pos[X]) + sq(betay*pos[Y]);
	sphere[X] = sqrt(r_s + sq(betaz*pos[Z]));
	r_s = sqrt(r_s);
	sphere[Y] = std::atan2(r_s, pos[Z]);
	sphere[Z] = std::atan2(betay*pos[Y], pos[X]) + M_PI;

	return sphere;
}

std::array<double, 3> Planet3D::sphereGrad(const std::array<double, 3>& pos_s, const Refractivity& n) const
{
	std::array<double, 3> grad_s;

	grad_s[X] = -n.n/atm.H;
	grad_s[Y] = m.n_tet(pos_s[Y], pos_s[Z]);
	grad_s[Z] = m.n_phi(pos_s[Y], pos_s[Z]);

	return exp(-(pos_s[X] - R)/atm.H)*grad_s;
}

void FlashMap::save(const std::string& filename){
	std::ofstream fp(filename, std::ofstream::out);

	fp << N << " " << S << std::endl;
	fp << L << " " << hh << " " << alpha << std::endl;
	fp << ray_hit << " " << ray_counter << std::endl;
	fp << Int1 << " " << Int2 << std::endl;

	for(int i = 0; i < N+1; ++i){
		for(int j = 0; j < N+1; ++j){
			fp << map[i][j] << " ";
		}
		fp << std::endl;
	}
}







std::array<double,3> rayTracing(const Planet3D& p, const double& H, const double& L, const std::array<double,2>& a, double& min, const double& dt)
{
	double n2, n1, r, k, tau = 0, st = sin(a[TET]);
	std::array<double, 3> vel, pos, dv;
	Refractivity n;
	bool inside = false;

	vel = p.rotate({st*cos(a[PHI]), st*sin(a[PHI]), cos(a[TET])});	//
	pos = p.rotate({H, 0, -L});										// rotate velocity and position to planet's reference frame

	if(!(k = p.checkRayHit(pos, vel))){			// check if the ray enters the planet
		vel = p.protate(vel);
		return {atan(vel[X]/vel[Z]), atan(vel[Y]/vel[Z]), 1}; // if not, return angles to the original reference frame
	}

	pos = k*vel + pos;		// update position to the atmosphere's entry point
	min = p.r2_max;

	do{
		p.updateRefrac(pos, n);

		n1 = 1/n.n;		//
		n2 = sq(n1);			// necessary constants

		//gn = p.n.df(pos);	// gradient of the refraction index
		k = n.gn*vel;			// inner product between the gradient and the veloicty unit vector

		dv = n2*(n.gn - k*vel);	// derivative of the velocity unit vector

		vel += dt*dv;		// update velocity
		pos += (n1*dt)*vel;	// update position

		if(!(r = p.checkDistance(pos))){
			min = p.R2;
			return {-100, -100, 0};	// check if inside the planet, if not return inside planet value {-100, -100}
		}
		if(r < p.r2_max) inside = true;

		if(min > r) min = r;

		tau += n1*exp(-(sqrt(r) - p.R)/p.atm.H);

	}while((r < p.r2_max || !inside) && r < 3*p.r2_max);	// check if outside the atmosphere

	vel = p.protate(p.arotate(vel));	// rotate velocity back to original reference frame
	return {atan(vel[X]/vel[Z]), atan(vel[Y]/vel[Z]), exp(-p.o*tau*dt)};	// return exit angles
}

/*FinalPos rayTracing2(const Planet3D& p, const double& H, const double& L, const std::array<double,2>& a, double& min, const double& dt)
{
	double n2, n1, r, k, tau = 0;
	std::array<double, 3> vel, pos, dv;
	Refractivity n;
	bool inside = false;

	vel = p.rotate({0, 0, 1});	//
	pos = p.rotate({a[X], a[Y], -2*p.r_max});										// rotate velocity and position to planet's reference frame

	if(!(k = p.checkRayHit(pos, vel))){			// check if the ray enters the planet
		vel = p.protate(vel);
		pos = p.protate(pos);
		return FinalPos(vel, pos, 1); // if not, return angles to the original reference frame
	}

	pos = k*vel + pos;		// update position to the atmosphere's entry point
	min = p.r2_max;

	do{
		p.updateRefrac(pos, n);

		n1 = 1/n.n;		//
		n2 = sq(n1);			// necessary constants

		//gn = p.n.df(pos);	// gradient of the refraction index
		k = n.gn*vel;			// inner product between the gradient and the veloicty unit vector

		dv = n2*(n.gn - k*vel);	// derivative of the velocity unit vector

		vel += dt*dv;		// update velocity
		pos += (n1*dt)*vel;	// update position

		if(!(r = p.checkDistance(pos))){
			min = p.R2;
			return FinalPos(vel, pos, -1);	// check if inside the planet, if not return inside planet value {-100, -100}
		}
		if(r < p.r2_max) inside = true;

		if(min > r) min = r;

		tau += n1*exp(-(sqrt(r) - p.R)/p.atm.H);

	}while((r < p.r2_max || !inside) && r < 3*p.r2_max);	// check if outside the atmosphere

	vel = p.protate(p.arotate(vel));	// rotate velocity back to original reference frame
	pos = p.protate(p.arotate(pos));	// rotate position back to original reference frame
	return FinalPos(vel, pos, exp(-p.o*tau*dt));
}*/

double rayTrace(const Planet3D& p, std::array<double,3>& pos, std::array<double,3>& vel, const double& dt)
{
	double k, n1, n2, r, tau = 0;
	std::array<double,3> dv;
	Refractivity n;
	bool inside = false;

	vel = p.rotate(vel);	//
	pos = p.rotate(pos);	// rotate velocity and position to planet's reference frame

	if(!(k = p.checkRayHit(pos, vel))){			// check if the ray enters the planet
		vel = p.protate(vel);
		pos = p.protate(pos);
		return 1; // if not, return angles to the original reference frame
	}

	pos = k*vel + pos;		// update position to the atmosphere's entry point

	do{
		p.updateRefrac(pos, n);

		n1 = 1/n.n;		//
		n2 = sq(n1);			// necessary constants

		//gn = p.n.df(pos);	// gradient of the refraction index
		k = n.gn*vel;			// inner product between the gradient and the veloicty unit vector

		dv = n2*(n.gn - k*vel);	// derivative of the velocity unit vector

		vel += dt*dv;		// update velocity
		pos += (n1*dt)*vel;	// update position

		if(!(r = p.checkDistance(pos))){
			return -1;	// check if inside the planet, if not return inside planet value {-100, -100}
		}
		if(r < p.r2_max) inside = true;

		tau += n1*exp(-(sqrt(r) - p.R)/p.atm.H);

	}while((r < p.r2_max || !inside) && r < 3*p.r2_max);	// check if outside the atmosphere

	vel = p.protate(p.arotate(vel));	// rotate velocity back to original reference frame
	pos = p.protate(p.arotate(pos));	// rotate position back to original reference frame
	return exp(-p.o*tau*dt);
}

FlashMap* mapGen(const Planet3D& p, const int& N, const double& S, const double& L, const double& hh)
{
	FlashMap *map = new FlashMap(N, S, L, hh, 0);
	int n_phi = 10*N;
	int n_tet = 10000*N/S, ray_counter = 0, ray_hit = 0;
	std::array<double, 2> angle;
	std::array<int, 2> ij;
	std::array<double, 3> ex, pos, vel;
	double max = acos(sqrt(1 - sq(p.r_max/L)));
	double min = acos(sqrt(1 - sq(p.R/L))), opt = 0;
	double dt = (max - min)/n_tet, scale = S/L;
	double dphi = 1/((double)n_phi-1), h = 2*scale/N;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> t;
	bool in_picture;

	start = std::chrono::system_clock::now();
	for(int i = 0; i < n_phi; ++i){
		angle[PHI] = 2*M_PI*i*dphi - M_PI/2;
		if(p.dev) angle[TET] = ((opt)?(0.97*opt+0.03*max):max) - hh*cos(angle[PHI])/L;
		else angle[TET] = (opt)?opt:max - hh*cos(angle[PHI])/L;
		in_picture = false;
		do{
			++ray_counter;
			//std::cout << "r = " << ray_counter << std::endl;
			//ex = rayTracing(p, hh, L, angle, mIn);
			min = sin(angle[TET]);
			pos = {0, 0, -L};
			vel = {min*cos(angle[PHI]), min*sin(angle[PHI]), cos(angle[TET])};
			min = rayTrace(p, pos, vel);

			ex = {atan(vel[X]/vel[Z]), atan(vel[Y]/vel[Z]), min};

			ij[X] = (int)(N*(ex[X] + scale + h + hh/L)/(2*scale));
			ij[Y] = (int)(N*(ex[Y] + scale + h)/(2*scale));

			if(ij[X] >= 0 && ij[Y] >= 0 && ij[X] < N+1 && ij[Y] < N+1){
				if(!in_picture){
					in_picture = true;
					if(!(opt))/* || p.dev))*/ {opt = angle[TET];}
				}
				map->map[ij[X]][ij[Y]] += ex[BRT];
				map->Int1 += dt*dphi*angle[TET];
				++ray_hit;
			}
			else if(in_picture) break;

			angle[TET] -= dt;
		}while(ex[X] != -100);
		Print(i+1 << " out of " << n_phi);
		Print("Ray count: " << ray_counter << std::endl);
	}
	end = std::chrono::system_clock::now();
	t = end - start;
	Print("Elapsed calculation time: " << t.count());
	map->Int2 = sq(h)*ray_hit/map->Int1;

	return map;
}

FlashMap* mapGen2(const Planet3D& p, const int& N, const double& S, const double& L, const double& hh)
{
	FlashMap *map = new FlashMap(N, S, L, hh, 0);
	FinalPos ex;
	int n_phi = 10*N, n_r = 50000*N/S;
	double intensity, dphi = 2*M_PI/n_phi, dr = (p.r_max - p.R)/n_r, h = 2*S/N, opt = 0;
	std::array<double, 3> pos, vel;
	std::array<double, 2> a;
	std::array<int, 2> ij;
	bool in_picture;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> t_end;

	start = std::chrono::system_clock::now();
	Print("Starting single threaded calculation . . .");
	for(double phi = -M_PI/2; phi < 3*M_PI/2; phi += dphi){
		in_picture = false;
		for(double r = (opt && !hh) ? opt:p.r_max; r > p.R; r -= dr){
			++(map->ray_counter);
			pos = {r*cos(phi), r*sin(phi), -2*p.r_max};
			vel = {0, 0, 1};
			intensity = rayTrace(p, pos, vel);

			//a = {r*cos(phi), r*sin(phi)};
			//ex = rayTracing2(p, 0, L, a, intensity);

			if(intensity == -1) break;
			//if(ex.intensity == -1) break;

			// Project into detector plane
			a = planeIntersect(pos, vel, L, 0, hh);
			//a = planeIntersect(ex.pos, ex.vel, L, 0);

			// Map into detector coordinates
			ij[X] = (int)(N*(a[X] + S + h))/(2*S);
			ij[Y] = (int)(N*(a[Y] + S + h))/(2*S);
			if(ij[X] >= 0 && ij[Y] >= 0 && ij[X] < N+1 && ij[Y] < N+1){
				if(!in_picture){
					in_picture = true;
					if(!opt){
						opt = r;
						Print("Optimal start found . . .");
					}
				}
				map->map[ij[X]][ij[Y]] += intensity;
				map->Int1 += r*dphi*dr;
				++(map->ray_hit);
			}
			else if(in_picture) break;
		}
		//Print(i+1 << " out of " << n_phi);
		//Print("Ray count: " << ray_counter << std::endl);
	}
	end = std::chrono::system_clock::now();
	t_end = end - start;
	Print("Elapsed calculation time: " << t_end.count());
	map->Int2 = sq(h)*map->ray_hit/map->Int1;
	//Print("Int1 = " << map->Int1);
	//Print("Int2 = " << map->Int2);

	return map;
}


FlashMap* mapThread(const Planet3D& p, const int& N, const double& S, const double& L, const double& hh, const int& thread_n)
{
	FlashMap *map = new FlashMap(N, S, L, hh, 0);
	std::vector<FlashMap*> map_t(0);
	std::vector<std::thread> threads(0);
	std::vector<double> thread_angle(thread_n+1);
	int n_phi = ((p.turb) ? 100:10)*N, n_r = 50*p.R*N/S;
	double dphi = 2*M_PI/n_phi, dr = (p.r_max - p.R)/n_r, opt = 0;
	double t_angle = 2*M_PI/thread_n, h = 2*S/N, h0;
	std::array<double,3> pos, vel;
	std::array<int,2> ij;
	std::array<double,2> a;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> t_end;

	start = std::chrono::system_clock::now();

	// allocate thread memory and positions
	Print("Allocating memory for " << thread_n << " threads . . .");
	thread_angle[0] = -M_PI/2;
	for(int i = 0; i < thread_n; ++i){
		map_t.push_back(new FlashMap(N,S,L,hh,0));
		thread_angle[i+1] = thread_angle[i] + t_angle;
	}

	// pre-process to calculate opt
	if(hh == 0){
		Print("Pre-calculating optimal starting distance . . .");
		h0 = 0.5*p.atm.H*(1 - p.betay*p.atm.H/(2*p.R))*log(2*M_PI*p.betay*sq(p.n_ref*L)/(p.R*p.atm.H)) - 9*p.betay*sq(p.atm.H)/(8*p.R);
		opt = h0*(1 + 2*p.betay*p.atm.H/p.R) + p.R/p.betay;
		for(double r = opt+0.5; r > p.R; r -= dr){
			++map->ray_counter;
			pos = {0, r, -2*p.r_max};
			vel = {0, 0, 1};
			rayTrace(p, pos, vel);

			// Project into detector plane
			a = planeIntersect(pos, vel, L, 0, hh);

			// Map into detector coordinates
			ij[X] = (int)(N*(a[X] + S + h))/(2*S);
			ij[Y] = (int)(N*(a[Y] + S + h))/(2*S);
			if(ij[X] >= 0 && ij[Y] >= 0 && ij[X] < N+1 && ij[Y] < N+1){
				opt = r;
				break;
			}
		}
	}

	// process threads
	Print("Starting threads . . .");
	for(int i = 0; i < thread_n; ++i){
		threads.push_back(std::thread(map_threaded, map_t[i], p, opt, dr, thread_angle[i], thread_angle[i+1], dphi));
	}
	// join threads
	for(int i = 0; i < thread_n; ++i){
		threads[i].join();
		Print("Joined thread " << i+1 << " . . .");
	}
	Print("Threading ended . . . ");
	Print("Post-processing results into a single map . . .");
	for(const FlashMap *mt : map_t){
		map->Int1 += mt->Int1;
		map->ray_hit += mt->ray_hit;
		map->ray_counter += mt->ray_counter;
	}
	map->Int2 = sq(h)*map->ray_hit/map->Int1;

	for(int i = 0; i <= N; ++i){
		for(int j = 0; j <= N; ++j){
			for(const FlashMap *mt : map_t){
				map->map[i][j] += mt->map[i][j];
			}
		}
	}
	end = std::chrono::system_clock::now();
	t_end = end - start;
	Print("Elapsed calculation time: " << t_end.count());

	for(const FlashMap *mt : map_t){
		delete mt;
	}
	Print("Hit percentage = " << (100.0*map->ray_hit)/map->ray_counter);
	return map;
}

void map_threaded(FlashMap* m, const Planet3D& p, const double& opt, const double& dr, const double& phi1, const double& phi2, const double& dphi)
{
	bool in_picture;
	std::array<double,3> pos, vel;
	std::array<int,2> ij;
	std::array<double,2> a;
	double intensity, h = 2*m->S/m->N;

	for(double phi = phi1; phi < phi2; phi += dphi){
		in_picture = false;
		for(double r = (opt && !m->hh) ? opt:p.r_max; r > p.R; r -= dr){
			++(m->ray_counter);
			pos = {r*cos(phi), r*sin(phi), -2*p.r_max};
			vel = {0, 0, 1};
			intensity = rayTrace(p, pos, vel);

			//a = {r*cos(phi), r*sin(phi)};
			//ex = rayTracing2(p, 0, L, a, intensity);

			if(intensity == -1) break;
			//if(ex.intensity == -1) break;

			// Project into detector plane
			a = planeIntersect(pos, vel, m->L, m->alpha, m->hh);
			//a = planeIntersect(ex.pos, ex.vel, L, 0);

			// Map into detector coordinates
			ij[X] = (int)(m->N*(a[X] + m->S + h))/(2*m->S);
			ij[Y] = (int)(m->N*(a[Y] + m->S + h))/(2*m->S);
			if(ij[X] >= 0 && ij[Y] >= 0 && ij[X] < m->N+1 && ij[Y] < m->N+1){
				if(!in_picture){
					in_picture = true;
				}
				m->map[ij[X]][ij[Y]] += intensity;
				m->Int1 += r*dphi*dr;
				++(m->ray_hit);
			}
			else if(in_picture) break;
		}
	}
}
