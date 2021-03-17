#include "auxiliar.hpp"


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

double sq(double a)
{
	return a*a;
}

double gx2D(const std::function<ddd>& f, const double& x, const double& z, const double& h)
{
	return (f(x+h, z) - f(x-h, z))/(2*h);
}

double gy2D(const std::function<ddd>& f, const double& x, const double& z, const double& h)
{
	return (f(x, z+h) - f(x, z-h))/(2*h);
}

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
