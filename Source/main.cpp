#include "tests.hpp"

int main(int argc, char* argv[])
{
	//double R = R_REF;			// planet radius
	double R = R_REF;
	double H = H_REF;			// atmospheric scale height
	//double r_max = R + 15*H;	// top layer of the atmosphere
	double L = 1e9;
	double N = 2e25;
	std::array<double, 3> r = {0, 0, -M_PI/2};
	FinalPos ex;
	double size = 50, l = 1500e-9; //l = (86.50415919381337933e-9)/1.3;
	int res = 100;
	int n_thread = 4;

	if(argc > 1){
		n_thread = std::atoi(argv[1]);
	}

	//Planet3D p(R, H, N, 0, OBF, OBF, r, "Config/map_0,002000_1,666667_0,100000");
	Planet3D p(R, H, N, N_REF, OBF, OBF, r, "Config/map");
	QApplication a(argc, argv);
	ImageGen w(p, l, L, size, res, 0, 4);

	w.show();
	a.exec();

	p.m.destroy();

	return 0;
}

/*int main(int argc, char* argv[])
{
	double R = R_REF;			// planet radius
	double H = H_REF;			// atmospheric scale height
	//double r_max = R + 15*H;	// top layer of the atmosphere
	double L = 1000000000;
	double N = 2e25;
	std::vector<double> n = {0.00001, 0.0002, 0.001, N_REF, 0.005};
	//Planet3D p(R, H, N, N_REF, 0, 0, r, "Config/map");
	int n_thread = 4;

	if(argc > 1){
		n_thread = std::atoi(argv[1]);
	}
	//bendHeight(n, R, H, L, N);
	//detecHeight(n, R, H, N);
	//horPeaks(N_REF, R, H, N, n_thread);
	//peakDist(N_REF, R, H, N, n_thread);
	//crossSec(N_REF, R, H, N);
	//ampWave(N_REF, R, H, N, n_thread);
	//resAmp(N_REF, R, H, L, N, n_thread);
	//ampHeight(N_REF, L, N, n_thread);
	//optim(N_REF, L, N);
	//atmosphericDensity(N_REF, R, H, L, n_thread);
	diamondSize(N_REF, 5, L, N, n_thread);
	//turbvs(R, H, L, N, n_thread);
	//turbsize(R, H, L, N, n_thread);

	//p.m.destroy();

	return 0;
}*/
