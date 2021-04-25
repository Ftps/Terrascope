#include "sphere_harm.hpp"

int main(int argc, char* argv[])
{
	int l_max = 875;
	int N = 1500;

	fillNumbers(l_max, N);

	verifyCoeff();

	/*int l = 875;
	long double aux;
	for(int m = 0; m < l; ++m){
		std::cout << (aux = coeffFact(l,m)) << std::endl;
		std::cout << (aux = sqrtl(aux)) << std::endl;
		std::cout << sqrtl(((2.0*l + 1.0)/(4.0*M_PI)))*aux << std::endl;
	}*/

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
