#include "terrascope3D.hpp"

Planet3D::Planet3D(const double &R, const double& r_max, const double& obf_x, const double& obf_y, const std::array<double, 2> r, const std::function<dddd>& n_ref) : R(R), R2(sq(R)), r_max(r_max), r2_max(sq(r_max)), betax(1/(1-obf_x)), betay(1/(1-obf_y))
{
	if(rot[TET] || rot[PHI]){
		dev = true;
		rot[X][X] = cos(r[TET]);
		rot[X][Y] = sin(r[TET])*sin(r[PHI]);
		rot[X][Z] = sin(r[TET])*cos(r[PHI]);

		rot[Y][X] = 0;
		rot[Y][Y] = cos(r[PHI]);
		rot[Y][Z] = -sin(r[PHI]);

		rot[Z][X] = -sin(r[TET]);
		rot[Z][Y] = cos(r[TET])*sin(r[PHI]);
		rot[Z][Z] = cos(r[TET])*cos(r[PHI]);

		rotate = [this](std::array<double, 3> pos_init){
			std::array<double, 3> pos = {0, 0, 0};
			for(int i = 0; i < 3; ++i){
				for(int j = 0; j < 3; ++j){
					pos[i] += rot[i][j]*pos_init[j];
				}
			}
			return pos;
		};

		n = [this, n_ref](double x, double y, double z){
			std::array<double, 3> pos = rotate({x,y,z});
			return n_ref(pos[X], pos[Y], pos[Z]);
		};
	}
	else{
		dev = false;
		n = n_ref;
	}

}

double Planet3D::checkDistance(double x, double y, double z) const
{
	double r;

	if(dev){
		std::array<double,3> pos = rotate({x,y,z});
		x = pos[X];
		y = pos[Y];
		z = pos[Z];
	}

	r = sq(betax*x) + sq(betay*y) + z*z;

	if(r < R2) return 0;
	else return r;
}

double Planet3D::checkRayHit(const double& L, const double& vx, const double& vy, const double& vz) const
{
	std::array<double, 3> vel = rotate({vx, vy, vz}), pos = rotate({0,0,-L});
	double a = sq(betax*vel[X]) + sq(betay*vel[Y]) + sq(vel[Z]);
	double b = sq(betax)*pos[X]*vel[X] + sq(betay)*pos[Y]*vel[Y] + pos[Z]*vel[Z];
	double c = sq(betax*pos[X]) + sq(betay*pos[Y]) + sq(pos[Z]) - r2_max;
	double det = sq(b) - a*c;

	/*Print("rxx = " << rot[X][X]);
	Print("rxy = " << rot[X][Y]);
	Print("rxz = " << rot[X][Z]);
	Print("ryx = " << rot[Y][X]);
	Print("ryy = " << rot[Y][Y]);
	Print("ryz = " << rot[Y][Z]);
	Print("rzx = " << rot[Z][X]);
	Print("rzy = " << rot[Z][Y]);
	Print("rzz = " << rot[Z][Z]);

	Print(L);
	Print("vx = " << vx <<", vy = " << vy << ", vz = " << vz);
	Print("x = " << pos[X] <<", y = " << pos[Y] << ", z = " << pos[Z]);
	Print("vx = " << vel[X] <<", vy = " << vel[Y] << ", vz = " << vel[Z]);
	Print(det);
	getchar();*/

	if(det < 0) return 0;

	return (-b - sqrt(det))/a;

}


std::array<double,2> rayTracing(const Planet3D& p, const double& L, const std::array<double,2>& a, const double& dz)
{
	double vx, vy, vz, x, y, z, n2, gx, gy, dvx, dvy, r;
	double ct = cos(a[TET]), st = sin(a[TET]), cp = cos(a[PHI]);
	double k;

	vx = st*cp;
	vy = st*sin(a[PHI]);
	vz = ct;

	if(!(k = p.checkRayHit(L, vx, vy, vz))) return a;

	x = k*vx;
	y = k*vy;
	z = k*vz - L;

	do{
		n2 = 1/(vz*sq(p.n(x, y, z)));
		gx = gx3D(p.n, x, y, z);
		gy = gy3D(p.n, x, y, z);

		k = gx*vx + gy*vy + gz3D(p.n, x, y, z)*vz;
		dvx = n2*(gx - k*vx);
		dvy = n2*(gy - k*vy);

		vx += dvx*dz;
		vy += dvy*dz;
		vz = sqrt(1 - vx*vx - vy*vy);

		x += vx*dz/vz;
		y += vy*dz/vz;
		z += dz;

		if(!(r = p.checkDistance(x, y, z))) return {-100, -100};

	}while(r < p.r2_max);

	return {atan(vx/vz), atan(vy/vz)};
}
