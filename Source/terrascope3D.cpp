#include "terrascope3D.hpp"

Planet3D::Planet3D(const double &R, const double& r_max, const double& obf, const std::function<dddd>& n) : R(R), r_max(r_max), obf(obf), n(n) {}





double gx3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x+h, y, z) - f(x-h, y, z))/(2*h);
}

double gy3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x, y+h, z) - f(x, y-h, z))/(2*h);
}

double gz3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x, y, z+h) - f(x, y, z-h))/(2*h);
}


std::array<double,2> rayTracing(const Planet3D& p, const double& L, const std::array<double,2>& a, const double& dz)
{
	double vx, vy, vz, x, y, z, n2, gx, gy, dvx, dvy, r;
	double beta = 1/(1 - p.obf);
	double R = p.R*p.R, r_max = p.r_max*p.r_max;
	double ct = cos(a[TET]), st = sin(a[TET]), cp = cos(a[PHI]);
	double k = (sq(beta)-1)*sq(st*cp) + 1;
	double det = sq(L*ct) + k*(L*L - r_max);

	if(det <= 0) return a;

	vx = st*cp;
	vy = st*sin(a[PHI]);
	vz = ct;

	k = - (L*ct + sqrt(det))/k;
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

		r = sq(beta*x) + y*y + z*z;

		if(r < R) return {0, 0};

	}while(r < r_max);

	return {atan(vx/vz), atan(vy/vz)};
}
