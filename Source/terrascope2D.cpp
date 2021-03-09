#include "terrascope2D.hpp"

Planet2D::Planet2D(const double& R, const double& r_max, const double& obf, const std::function<ddd>& n) : R(R), r_max(r_max), obf(obf), n(n) {}
Ray2D::Ray2D(const Planet2D& p) : Planet2D(p), x(0), y(0) {}

int Ray2D::ray_tracer(const double& Ll, const double& Aa_entry, const double& dx)
{
	double vx, vy, r, n2, dv, gv, ggy;
	double rr = 1 - obf;
	double ta = tan(Aa_entry);
	double ta2 = ta*ta;
	double det = ta2*ta2*Ll*Ll - (rr*rr+ta2)*(ta2*Ll*Ll - r_max*r_max);

	i = 1;
	a_entry = 180*Aa_entry/M_PI;
	L = Ll;

	if(det < 0){
		x.resize(2);
		y.resize(2);
		x[0] = -L;
		y[0] = 0;
		x[1] = 5*r_max;
		y[1] = ta*(x[1] + L);
		a_exit = a_entry;
		return 0;	// ray does not enter the atmosphere, angle remains constant
	}

	x.resize(2*r_max/dx);
	y.resize(2*r_max/dx);
	x[0] = -L;
	y[0] = 0;

	x[1] = -(ta2*L + sqrt(det))/(rr*rr + ta2); 	//
	y[1] = ta*(x[1] + L);					//
											// initial conditions at the atmosphere's boundary
	vx = cos(Aa_entry);						//
	vy = sin(Aa_entry);						//

	do{
		++i;

		n2 = 1/(vx*pwr(n(x[i-1], y[i-1]), 2));
		ggy = gy2D(n, x[i-1], y[i-1]);
		gv = gx2D(n, x[i-1], y[i-1])*vx + vy*ggy;
		dv = n2*(ggy - vy*gv);
		vy += dv*dx;
		vx = sqrt(1 - vy*vy);

		x[i] = x[i-1] + dx;
		y[i] = y[i-1] + vy*dx/vx;
		r = pwr(rr*x[i], 2) + y[i]*y[i];

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

int Ray2D::ray_tracer_hor(const double& Y, const double& dx)
{
	double vx, vy, r, n2, dv, gv, ggy;
	double rr = 1 - obf;

	a_entry = 0;
	L = INFINITY;
	i = 1;

	if((Y > r_max) || (Y < -r_max)){
		x.resize(2);
		y.resize(2);
		x[0] = -5*r_max;
		y[0] = Y;
		x[1] = 5*r_max;
		y[1] = Y;
		a_exit = a_entry;
		return 0;
	}

	x.resize(2*r_max/dx);
	y.resize(2*r_max/dx);
	x[0] = -5*r_max;
	y[0] = Y;

	x[1] = -sqrt(r_max*r_max - Y*Y)/rr; 		//
	y[1] = Y;								//
											// initial conditions at the atmosphere's boundary
	vx = 1;									//
	vy = 0;									//

	do{
		++i;

		n2 = 1/(vx*pwr(n(x[i-1], y[i-1]), 2));
		ggy = gy2D(n, x[i-1], y[i-1]);
		gv = gx2D(n, x[i-1], y[i-1])*vx + vy*ggy;
		dv = n2*(ggy - vy*gv);
		vy += dv*dx;
		vx = sqrt(1 - vy*vy);

		x[i] = x[i-1] + dx;
		y[i] = y[i-1] + vy*dx/vx;
		r = pwr(rr*x[i], 2) + y[i]*y[i];

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




double fisqrt(double n) // Fast inverse square root algorithim for double
{
	union{
		double y;
		long i;
	} caster;
	double x2;
	const double thf = 1.5;

	x2 = n*0.5;
	caster.y  = n;                       // evil floating point bit level hacking
	caster.i  = 0x5fe6eb50c7b537a9 - ( caster.i >> 1 );       // log magic with exponents

	caster.y  = caster.y * ( thf - ( x2 * caster.y * caster.y ) );   // 1st iteration (accuracy ~1e-3)
	caster.y  = caster.y * ( thf - ( x2 * caster.y * caster.y ) );   // 2nd iteration (accuracy ~1e-6)
	//caster.y  = caster.y * ( thf - ( x2 * caster.y * caster.y ) );   // 3rd iteration, this can be removed for speed (accuracy ~1e-11)

	return caster.y;
}

double pwr(double a, int n)
{
	if(n == 1) return a;
	else return pwr(a, n-1);
}



double gx2D(const std::function<ddd>& f, const double& x, const double& y, const double& h)
{
	return (f(x+h, y) - f(x-h, y))/(2*h);
}

double gy2D(const std::function<ddd>& f, const double& x, const double& y, const double& h)
{
	return (f(x, y+h) - f(x, y-h))/(2*h);
}



double ray_tracer2D(const std::function<ddd>& n, const double& R, const double& r_max, const double& obf, const double& L, const double& a_init, const double& dx)
{
	double x, y, vx, vy, r, n2, dv, gv, ggy;
	double rr = 1 - obf;
	double ta = tan(a_init);
	double ta2 = ta*ta;
	double det = ta2*ta2*L*L - (rr*rr+ta2)*(ta2*L*L - r_max*r_max);

	if(det < 0) return a_init;	// ray does not enter the atmosphere, angle remains constant

	x = -(ta2*L + sqrt(det))/(rr*rr + ta2); 	//
	y = ta*(x + L);							//
											// initial conditions at the atmosphere's boundary
	vx = cos(a_init);						//
	vy = sin(a_init);						//

	//r = r_max*r_max;

	do{
		n2 = 1/(vx*pwr(n(x, y), 2));
		ggy = gy2D(n, x, y);
		gv = gx2D(n, x, y)*vx + vy*ggy;
		dv = n2*(ggy - vy*gv);
		vy += dv*dx;
		vx = sqrt(1 - vy*vy);

		x += dx;
		y += vy*dx/vx;
		r = pwr(x*rr, 2) + y*y;

		if(r < R*R) return -100;
	}while(r < r_max*r_max);

	return atan(vy/vx);
}

double ray_tracer2D_hor(const std::function<ddd>& n, const double& R, const double& r_max, const double& obf, const double& Y, const double& dx)
{
	double x, y, vx, vy, r, n2, dv, gv, ggy;
	double rr = 1 - obf;

	if((Y > r_max) || (-Y > r_max)) return -1;

	x = -sqrt(r_max*r_max - Y*Y)/rr;
	y = Y;
	vx = 1;
	vy = 0;

	do{
		n2 = 1/(vx*pwr(n(x, y), 2));
		ggy = gy2D(n, x, y);
		gv = gx2D(n, x, y)*vx + vy*ggy;
		dv = n2*(ggy - vy*gv);
		vy += dv*dx;
		vx = sqrt(1 - vy*vy);

		x += dx;
		y += vy/vx*dx;
		r = pwr(x*rr, 2) + y*y;

		if(r < R*R) return -2;
	}while(r < r_max*r_max);

	return x - y*vx/vy;
}

double focalPoint(const Planet2D& p, const int& N)
{
	double h_min = p.R, dh = 0.01, a_min, h, dx_max = 0;
	double *x, *dx;
	int n, i_max = 0;

	do{
		h_min += dh;
		a_min = ray_tracer2D_hor(p.n, p.R, p.r_max, p.obf, h_min);
	}while(a_min == -2);

	n = floor((p.r_max - h_min)/dh);

	x = new double[n];
	dx = new double[n-1];

	x[0] = a_min;
	h = h_min;

	for(int i = 1; i < n; ++i){
		h += dh;
		x[i] = ray_tracer2D_hor(p.n, p.R, p.r_max, p.obf, h);
		dx[i-1] = dh/(x[i] - x[i-1]);
		if(dx_max < dx[i-1]){
			i_max = i-1;
			dx_max = dx[i-1];
		}
	}

	delete x;
	delete dx;

	Print(a_min);
	Print(i_max);
	Print(dx_max);
	Print(0.5*(x[i_max+1] + x[i_max]));

	return 0.5*(x[i_max+1] + x[i_max]);
}
