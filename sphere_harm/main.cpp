#include "sphere_harm.hpp"
#include "../Headers/gnuplot-iostream.hpp"

int main(int argc, char* argv[])
{

	int l = 500;
	int N_min = 1000;
	int N_max = 3000;

	fillNumbers(l, N_min, N_max);

	verifyCoeff();

	/*Gnuplot gp;
	auto plots = gp.plotGroup();
	int l = 500, m = 0;
	int N = 3000;
	double *y = new double[N];
	long double dx = M_PI/N, aux, coeff;
	std::vector<std::pair<long double, long double>> xy;
	std::vector<std::complex<double>> dft;

	aux = sqrtl(coeffFact(l, m));
	coeff = sqrtl(((2.0*l + 1.0)/(4.0*M_PI)))*aux;

	for(int i = 0; i < N; ++i){
		y[i] = coeff*ALpoly(cos(i*dx-M_PI/2), m, l, aux);
	}

	dft = ditfft(y, N);

	for(int i = 0; i < N; ++i){
		//xy.emplace_back(i/M_PI, 2*std::abs(dft[i])/N);
		//std::cout << dft[i] << std::endl;
		xy.emplace_back(i*dx - M_PI/2, coeff*ALpoly(cos(i*dx-M_PI/2), m, l, aux));
	}

	plots.add_plot1d(xy, "with lines title 'F(P_{" + std::to_string(l) + "}^{" + std::to_string(m) + "})'");
	//gp << "set xrange [0:" + std::to_string(N/M_PI) + "]\n";
	gp << "set xrange [-1.58:1.58]\n";
	gp << plots;

	std::cout << "Max frequency = " << (aux = maxfreq(m, l)) << " rad^-1" << std::endl;
	std::cout << "Minimum number of points = " << 10*M_PI*aux << std::endl;
	std::cout << "Smallest wavelength = " << 2200/aux << " km" << std::endl;

	delete y;*/

	return 0;
}

/*
Gnuplot gp;
auto plots = gp.plotGroup();
int N = 1000;
std::vector<std::pair<double, double>> xy;
double h = 2.0/(N-1), aux = 0;
int m = 1;


gp << "set xrange [-1:1]\nset yrange [-4:3]\n";

for(int l = m; l < m+5; ++l){
	for(int i = 0; i < N; ++i){
		xy.emplace_back(i*h-1, ALpoly2(i*h-1, m, l, aux));
	}

	plots.add_plot1d(xy, "with lines title 'P_" + std::to_string(l) + "^" + std::to_string(m) + "'");
	xy.clear();
}

gp << plots;

auto plots = gp.splotGroup();
double r, tet, phi, r_max = 0;
int N = 500, m = 1, l = 500;
std::vector<std::vector<double>> xx(N), yy(N), zz(N);
double dtet = M_PI/(N-1), dphi = 2*M_PI/(N-1), aux;

gp << "set xrange [-1:1]\nset yrange [-1:1]\nset zrange [-1:1]\n";
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
	std::cout << "i = " << i << std::endl;
	for(int j = 0; j < N; ++j){
		phi = j*dphi;
		r = std::abs(sphereHarm(tet, phi, m, l, aux));
		xx[i][j] = r*sin(tet)*cos(phi);
		yy[i][j] = r*sin(tet)*sin(phi);
		zz[i][j] = r*cos(tet);

		if(std::abs(r) > r_max) r_max = std::abs(r);
	}
}

for(int i = 0; i < N; ++i){
	for(int j = 0; j < N; ++j){
		xx[i][j] = xx[i][j]/r_max;
		yy[i][j] = yy[i][j]/r_max;
		zz[i][j] = zz[i][j]/r_max;
	}
}

plots.add_plot2d(std::make_tuple(xx, yy, zz), "title 'Re(Y_" + std::to_string(l) + "^" + std::to_string(m) + ")' w pm3d");
r_max = 0;

for(int i = 0; i < N; ++i){
	tet = i*dtet;
	for(int j = 0; j < N; ++j){
		phi = j*dphi;
		r = sphereHarmIm(tet, phi, m, l);
		xx[i][j] = r*sin(tet)*cos(phi);
		yy[i][j] = r*sin(tet)*sin(phi);
		zz[i][j] = r*cos(tet);

		if(std::abs(r) > r_max) r_max = std::abs(r);
	}
}

for(int i = 0; i < N; ++i){
	for(int j = 0; j < N; ++j){
		xx[i][j] = xx[i][j]/r_max;
		yy[i][j] = yy[i][j]/r_max;
		zz[i][j] = zz[i][j]/r_max;
	}
}

plots.add_plot2d(std::make_tuple(xx, yy, zz), "title 'Im(Y_" + std::to_string(l) + "^" + std::to_string(m) + ")' w pm3d");

gp << plots;


return 0;

 */
