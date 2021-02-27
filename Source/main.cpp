#include "drawRay.hpp"

int main(int argc, char* argv[])
{
	double R = 2574;			// planet radius
	double H = 25;				// atmospheric scale height
	double r_max = R + 10*H;	// top layer of the atmosphere
	double L = 500000;
	double a_init = 0.5*(R + r_max)/L;

	if(argc == 1){
		std::function<ddd> n = [](double x, double y){ return 1 + 0.5*exp(-(sqrt(x*x + y*y)-2574)/25); };
		Planet2D p(R, r_max, n);

		QApplication a(argc, argv);

		DrawRay w(p, L, a_init);
		w.show();

		Print("Surface level refractivity of 0.5");
		a.exec();
	}
	else{
		std::function<ddd> n = [](double x, double y){ return 1 + 0.15*exp(-(sqrt(x*x + y*y)-2574)/25); };
		Planet2D p(R, r_max, n);
		int N = 10000000, iter, tot = 0;
		double a, h = (1-0.98)*a_init/(double)N;

		Ray2D ray(p);

		ray.ray_tracer(L, a_init);

		for(int i = 0; i <= N; ++i){
			a = a_init - i*h;
			if(!(i % 100000)) Print(i);
			ray_tracer2D(n, R, r_max, L, a, iter, 10);
			tot += iter;
		}


		Print("Total rays: " << N);
		Print("Total Iterations: " << tot);
		Print("Average Iterations per ray: " << tot/(double)N);

		Print("\nSingle ray test:");
		Print("Entry ray angle: " << ray.a_entry << " º");
		Print("Exit ray angle: " << ray.a_exit << " º");
		Print("Surface level refractivity of 0.15");
	}


	return 0;
}

/*
int main(int argc, char* argv[])
{
	double R = 2574;			// planet radius
	double H = 25;				// atmospheric scale height
	double r_max = R + 10*H;	// top layer of the atmosphere
	//double N = 0.000293;		// surface refractivity
	double L = 500000;
	double a_init = 0.5*(R + r_max)/L;
	//int N = 350000, iter, tot = 0;
	//double a, h = (1.05-0.96)*a_init/(double)N;
	int i;

	std::function<ddd> n = [](double x, double y){ return 1 + 0.1*exp(-(sqrt(x*x + y*y)-2574)/25); };

	Planet2D p(R, r_max, n);
	Ray2D ray(p);

	ray.ray_tracer(L, a_init);

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
	Print("Entry ray angle: " << ray.a_entry << " º");
	Print("Exit ray angle: " << ray.a_exit << " º");

	return 0;
}
*/
