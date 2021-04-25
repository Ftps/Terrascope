/*#include "drawCFM.hpp"

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
*/

#include "sphere_harm.hpp"
#include "gnuplot-iostream.hpp"

int main()
{

	/*Gnuplot gp;
	auto plots = gp.plotGroup();
	int N = 1200;
	std::vector<std::pair<double, double>> xy;
	std::vector<Poly> poly;
	double h = 2.0/(N-1);
	int m = 32;
	SphereHarm s(77);

	for(int l = 74; l < 78; ++l){
		poly.push_back(Poly(l));

		for(int i = 0; i < N; ++i){
			xy.emplace_back(i*h-1, s.Clm[l-1][m]*poly[l-74](m, i*h-1));
		}

		plots.add_plot1d(xy, "with lines title 'P_{" + std::to_string(l) + "}^{" + std::to_string(m) + "}'");
		xy.clear();
	}

	for(int i = 0; i < poly[0].N; ++i){
		std::cout << poly[0].poly[31][i] << " " << poly[0].poly[32][i] << std::endl;
	}

	std::cout << s.Clm[76][31] << " " << s.Clm[76][32] << std::endl;

	gp << "set xrange [-1:1]\n";
	gp << plots;*/

	Gnuplot gp;
	auto plots = gp.splotGroup();
	double r, tet, phi, r_max = 0;
	int N = 1000, m = -55, l = 77;
	std::vector<std::vector<double>> xx(N), yy(N), zz(N);
	double dtet = M_PI/(N-1), dphi = 2*M_PI/(N-1);
	SphereHarm s(l);
	std::string ran;

	gp << "set xlabel 'x'\nset ylabel 'y'\nset zlabel 'z'\n";
	gp << "set size 0.5,1.0\n";
	//gp << "set style fill transparent solid 0.5\n";
	gp << "set style fill\n";
	gp << "set view 100, 315, 1, 1\n";
	//gp << "set view 0, 0, 1, 1\n";

	for(int i = 0; i < N; ++i){
		xx[i].resize(N);
		yy[i].resize(N);
		zz[i].resize(N);
		tet = i*dtet;
		//std::cout << "i = " << i << std::endl;
		for(int j = 0; j < N; ++j){
			phi = j*dphi;
			r = std::abs(s(tet, phi, m, l));
			xx[i][j] = r*sin(tet)*cos(phi);
			yy[i][j] = r*sin(tet)*sin(phi);
			zz[i][j] = r*cos(tet);

			if(std::abs(r) > r_max) r_max = std::abs(r);
		}
	}

	ran = "[-" + std::to_string(r_max) + ":" +  std::to_string(r_max) + "]";
	std::cout << r_max << std::endl;
	std::cout << ran << std::endl;
	gp << "set xrange " + ran + "\nset yrange "+ ran +"\nset zrange " + ran + "\n";


	plots.add_plot2d(std::make_tuple(xx, yy, zz), "title 'Im(Y_{" + std::to_string(l) + "}^{" + std::to_string(m) + "})' w pm3d");

	gp << plots;

	return 0;
}
