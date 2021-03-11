#include "drawRay.hpp"
#include <chrono>

#define N_REF 0.000293
#define R_REF 6376
#define H_REF 8.5
#define OBF 0.003292568
#define RR (1 - OBF)

int main(int argc, char* argv[])
{
	double R = R_REF;			// planet radius
	double H = H_REF;			// atmospheric scale height
	double r_max = R + 15*H;	// top layer of the atmosphere
	double obf = OBF;
	double L = 500000;
	double a_init = 0.5*(R + r_max)/L;
	std::function<ddd> n = [](double x, double y){ return 1 + N_REF*exp(-(sqrt(sq(RR*x) + y*y)-R_REF)/H_REF); };
	Planet2D p(R, r_max, obf, n);

	if(argc == 1){
		int type;
		char b[2];

		QApplication a(argc, argv);

		try {
			Print("Input graph number (1 - ray density, 2 - multiple horizontal rays, 3 - single ray)");
			fgets(b, 2, stdin);
			type = std::atoi(b);
		}
		catch(...) {
			type = 1;
		}

		DrawRay w(p, L, a_init, type);
		w.show();

		a.exec();
	}
	else if(!strcmp(argv[1], "1")){
		int N = 10000000;
		double a, h = (1-0.98)*a_init/(double)N;

		Ray2D ray(p);

		ray.ray_tracer(L, a_init);

		for(int i = 0; i <= N; ++i){
			a = a_init - i*h;
			if(!(i % 100000)) Print(i);
			ray_tracer2D(n, R, r_max, obf, L, a, 10);
		}

		Print("\nSingle ray test:");
		Print("Entry ray angle: " << ray.a_entry << " º");
		Print("Exit ray angle: " << ray.a_exit << " º");
		Print("Surface level refractivity of " << N_REF);
	}
	else if(!strcmp(argv[1], "2")){
		int N = 1000;

		Print(focalPoint(p, N));
	}
	else if(!strcmp(argv[1], "3")){ // Testing purposes
		double h = R + 0.01;

		Print(ray_tracer2D_hor(n, R, r_max, obf, h, 1));
	}

	return 0;
}
