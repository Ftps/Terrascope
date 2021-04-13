#include "drawCFM.hpp"

#define N_REF 0.000293
#define R_REF 2500
#define H_REF 25
#define OBF 0.003292568
#define OBFY 0.003292568
#define OBFZ 0.003292568

int main(int argc, char* argv[])
{
	double R = R_REF;			// planet radius
	double H = H_REF;			// atmospheric scale height
	//double r_max = R + 15*H;	// top layer of the atmosphere
	double L = 1500000;
	double N = 2e25;
	FDIV n = generateRefFuncs(N_REF, H, R, OBFY, OBFZ);
	std::array<double, 3> r = {0, 0, 0};
	double size = 50, l = 200e-9;
	int res = 200;

	Planet3D p(R, H, N, OBFY, OBFZ, r, n);
	QApplication a(argc, argv);
	ImageGen w(p, l, L, size, res);

	w.show();
	a.exec();

	return 0;
}
