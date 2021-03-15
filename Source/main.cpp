#include "drawRay.hpp"
#include "terrascope3D.hpp"
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
	std::function<ddd> n = [](double x, double y){ return 1 + N_REF*exp(-(sqrt(sq((1/RR)*x) + y*y)-R_REF)/H_REF); };
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
			//if(!(i % 100000)) Print(i);
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
		std::function<dddd> n2 = [](double x, double y, double z){ return 1 + N_REF*exp(-(sqrt(sq((1/RR)*x) + y*y + z*z)-R_REF)/H_REF); };
		Planet3D p2(R, r_max, obf, n2);
		int N = 10000000;
		double a, h = (1-0.98)*a_init/(double)N;
		std::array<double, 2> init = {a_init, 0};
		std::array<double, 2> ex = rayTracing(p2, L, init);

		for(int i = 0; i < N; ++i){
			a = a_init - i*h;
			if(!(i % 100000)) Print(i);
			rayTracing(p2, L, {a,0});
		}

		Print("X angle: " << 180*init[X]/M_PI << "; Y angle: " << init[Y]);
		Print("X angle: " << 180*ex[X]/M_PI << "; Y angle: " << ex[Y]);
	}
	else if(!strcmp(argv[1], "4")){
		std::function<dddd> n2 = [](double x, double y, double z){ return 1 + N_REF*exp(-(sqrt(sq((1/RR)*x) + y*y + z*z)-R_REF)/H_REF); };
		Planet3D p2(R, r_max, obf, n2);
		std::array<double, 2> init = {a_init, 0}, ex;
		double e;

		Print("2D");
		e = ray_tracer2D(n, R, r_max, obf, L, a_init);

		Print("3D");
		ex = rayTracing(p2, L, init);

		Print("2D angle:" << 180*e/M_PI);
		Print("X angle: " << 180*ex[X]/M_PI << "; Y angle: " << ex[Y]);
	}

	return 0;
}
