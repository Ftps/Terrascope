#include "terrascope2D.hpp"

Planet2D::Planet2D(const double& R, const double& r_max, const std::function<ddd>& n) : R(R), r_max(r_max), n(n) {}
Ray2D::Ray2D(const Planet2D& p) : Planet2D(p), x(0), y(0) {}

int Ray2D::ray_tracer(const double& Ll, const double& Aa_entry, const double& dx)
{
	double vx, vy, r, n2, dv, gv, ggy;
	double ta = tan(Aa_entry);
	double ta2 = ta*ta;
	double det = ta2*ta2*Ll*Ll - (1+ta2)*(ta2*Ll*Ll - r_max*r_max);

	i = 1;
	a_entry = 180*Aa_entry/M_PI;
	L = Ll;

	if(det < 0){
		x.resize(2);
		y.resize(2);
		x[1] = 5*r_max;
		y[1] = ta*(x[1] + L);
		a_exit = a_entry;
		return 0;	// ray does not enter the atmosphere, angle remains constant
	}

	x.resize(2*r_max/dx);
	y.resize(2*r_max/dx);
	x[0] = -L;
	y[0] = 0;

	x[1] = -(ta2*L + sqrt(det))/(1 + ta2); 	//
	y[1] = ta*(x[1] + L);						//
											// initial conditions at the atmosphere's boundary
	vx = cos(Aa_entry);						//
	vy = sin(Aa_entry);						//

	do{
		++i;

		n2 = 1/(vx*pwr(n(x[i-1], y[i-1]), 2));
		ggy = gy2D(n, x[i-1], y[i-1]);
		gv = vx*gx2D(n, x[i-1], y[i-1]) + vy*ggy;
		dv = n2*(ggy - vy*gv);
		vy += dv*dx;
		vx = sqrt(1 - vy*vy);

		x[i] = x[i-1] + dx;
		y[i] = y[i-1] + vy*dx/n(x[i-1],y[i-1]);
		r = x[i]*x[i] + y[i]*y[i];

		if(r < R*R){
			x.resize(i+1);
			y.resize(i+1);
			a_exit = nan("");
			return -1;
		}
	}while(r < r_max*r_max);


	a_exit = 180*atan(vy/vx)/M_PI;
	x.resize(i+2);
	y.resize(i+2);
	x[i+1] = 5*r_max;
	y[i+1] = y[i] + (vy/vx)*(x[i+1] - x[i]);

	return 1;
}

double fisqrt(double n)
{
	long i;
	double x2, y;
	const double threehalfs = 1.5;

	x2 = n*0.5;
	i  = * ( long * ) &n;                       // evil floating point bit level hacking
	i  = 0x5fe6eb50c7b537a9 - ( i >> 1 );               // what the fuck?
	y  = * ( double* ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

	return y;
}

double pwr(double a, int n)
{
	if(n == 1) return a;
	else return pwr(a, n-1);
}



double gx2D(std::function<double(double, double)> f, double x, double y, double h)
{
	return (f(x+h, y) - f(x-h, y))/(2*h);
}

double gy2D(std::function<double(double, double)> f, double x, double y, double h)
{
	return (f(x, y+h) - f(x, y-h))/(2*h);
}



double ray_tracer2D(const std::function<double(double, double)>& n, const double& R, const double& r_max, const double& L, const double& a_init, int& i, const double& dx)
{
	double x, y, vx, vy, r, n2, dv, gv, ggy;
	double ta = tan(a_init);
	double ta2 = ta*ta;
	double det = ta2*ta2*L*L - (1+ta2)*(ta2*L*L - r_max*r_max);

	i = 1;

	if(det < 0) return a_init;	// ray does not enter the atmosphere, angle remains constant

	x = -(ta2*L + sqrt(det))/(1 + ta2); 	//
	y = ta*(x + L);						//
											// initial conditions at the atmosphere's boundary
	vx = cos(a_init);						//
	vy = sin(a_init);						//

	//r = r_max*r_max;

	do{
		++i;
		//LOG
		n2 = 1/(vx*pwr(n(x, y), 2));
		ggy = gy2D(n, x, y);
		gv = vx*gx2D(n, x, y) + vy*ggy;
		dv = n2*(ggy - vy*gv);
		vy += dv*dx;
		vx = sqrt(1 - vy*vy);

		x += dx;
		y += vy*dx/n(x,y);
		r = x*x + y*y;
		//print(n2);
		//print(dv);

		if(r < R*R) return -100;
	}while(r < r_max*r_max);

	return atan(vy/vx);
}
