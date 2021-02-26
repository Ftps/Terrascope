#include <functional>
#include <iostream>
#include <cmath>
#include <cstdlib>

#define LOG {std::cout << "IN LINE " << __LINE__  << " OF FILE " << __FILE__ << std::endl; fflush(stdout);}
#define Print(X) std::cout << X << std::endl

double fexp(double n)
{
	return 1 + n*(1 + n*(1/2 + n*(1/6 + n/24)));
}

double fsqrt(double n)
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

	return 1/y;
}

double pwr(double a, int n)
{
	if(n == 1) return a;
	else return pwr(a, n-1);
}

double gx2D(std::function<double(double, double)> f, double x, double y, double h = 0.001)
{
	return (f(x+h, y) - f(x-h, y))/(2*h);
}

double gy2D(std::function<double(double, double)> f, double x, double y, double h = 0.001)
{
	return (f(x, y+h) - f(x, y-h))/(2*h);
}

double ray_tracer2D(std::function<double(double, double)> n, double R, double r_max, double L, double a_init, int& i, double dx = 1)
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
