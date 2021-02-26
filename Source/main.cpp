#include "terrascope2D.hpp"

int main(int argc, char* argv[])
{
	double R = 2574;			// planet radius
	double H = 25;				// atmospheric scale height
	double r_max = R + 10*H;	// top layer of the atmosphere
	//double N = 0.000293;		// surface refractivity
	double L = 500000;
	double a_init = 0.5*(R + r_max)/L;

	std::function<double(double, double)> n = [](double x, double y){ return 1 + 0.1*exp(-(sqrt(x*x + y*y)-2574)/25); };

	std::cout << "Entry ray angle: " << 180*a_init/3.14159 << " º" << std::endl;
	std::cout << "Exit ray angle: " << 180*ray_tracer2D(n, R, r_max, L, a_init)/3.14159 << " º" << std::endl;

	return 0;
}
