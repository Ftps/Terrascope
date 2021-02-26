#include "terrascope2D.hpp"
#include "qcustomplot.hpp"

int main(int argc, char* argv[])
{
	double R = 2574;			// planet radius
	double H = 25;				// atmospheric scale height
	double r_max = R + 10*H;	// top layer of the atmosphere
	//double N = 0.000293;		// surface refractivity
	double L = 500000;
	double a, a_init = 0.5*(R + r_max)/L;
	int N = 350000, iter, tot = 0;
	double h = (1.05-0.96)*a_init/(double)N;

	std::function<double(double, double)> n = [](double x, double y){ return 1 + 0.1*exp(-(sqrt(x*x + y*y)-2574)/25); };

	for(int i = 0; i <= N; ++i){
		a = 1.05*a_init - i*h;
		if(!(i % 1000)) Print(i);
		ray_tracer2D(n, R, r_max, L, a, iter, 10);
		tot += iter;
	}


	Print("Total rays: " << N);
	Print("Total Iterations: " << tot);
	Print("Average Iterations per ray: " << tot/(double)N);

	Print("\nSingle ray test:");
	std::cout << "Entry ray angle: " << 180*a_init/3.14159 << " º" << std::endl;
	std::cout << "Exit ray angle: " << 180*ray_tracer2D(n, R, r_max, L, a_init, iter)/3.14159 << " º" << std::endl;

	return 0;
}
