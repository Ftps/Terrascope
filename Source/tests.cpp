#include "tests.hpp"

void bendHeight(const std::vector<double>& n, const double& R, const double& H, const double& L, const double& N)
{
	int rr = 10000;
	std::array<double, 3> r = {0, 0, 0};
	double dh = (15.0*H)/rr, last = 0, hh;
	std::array<double, 3> res;

	Gnuplot gp1, gp2;
	auto plot = gp1.plotGroup();
	std::vector<std::pair<double, double>> xy;

	for(int i = 0; i < (int)n.size(); ++i){
		Planet3D p(R, H, N, n[i], 0, 0, r, "Config/map");
		for(double h = R + 14*H; h > R; h -= dh){
			hh = 0;
			res = rayTracing(p, h, L, {0,0}, hh);
			if(res[X] == -100){
				break;
			}
			if(std::abs(res[X])*180/M_PI > last) last = std::abs(res[X]*180/M_PI);
			xy.emplace_back(sqrt(hh) - R, std::abs(res[X])*180/M_PI);
		}
		plot.add_plot1d(xy, "with lines title 'n_0 = " + ST(n[i]) + "'");
		xy.clear();
		p.m.destroy();
	}

	gp1 << "set xrange [0:" << ST(14*H) << "]\n";
	gp1 << "set xlabel \"Min. height - km\"\n";
	gp1 << "set yrange [1e-11:15]\n";
	gp1 << "set ylabel \"Deflection angle - º\"\n";
	gp1 << "set logscale y\n";
	gp1 << plot;

	gp2 << "set xrange [0:" << ST(14*H) << "]\n";
	gp2 << "set xlabel \"Min. height - km\"\n";
	gp2 << "set yrange [0:10]\n";
	gp2 << "set ylabel \"Deflection angle - º\"\n";
	gp2 << plot;
}

void detecHeight(const std::vector<double>& n, const double& R, const double& H, const double& N)
{
	double LL;
	std::array<double, 3> r = {0, 0, 0};
	double hh, longest = 0;
	double tet_init, tet, dtet;
	std::array<double, 3> res;
	std::array<double, 2> angle = {0,0};

	Gnuplot gp;
	auto plot = gp.plotGroup();
	std::vector<std::pair<double, double>> xy;

	for(int i = 0; i < (int)n.size(); ++i){
		LL = 1e4;
		Planet3D p(R, H, N, n[i], 0, 0, r, "Config/map");
		do{
			tet_init = asin((R + 14.8*H)/LL);
			dtet = (tet_init - asin(R/LL))/1000.0;
			for(tet = tet_init; tet > 0; tet -= dtet){
				angle[TET] = tet;
				res = rayTracing(p, 0, LL, angle, hh);
				if(res[X] < 0) break;
			}
			xy.emplace_back(LL, sqrt(hh) - R);
			if(longest < sqrt(hh)-R) longest = sqrt(hh)-R;
		}while((LL = 1.15*LL) < 1e10);
		plot.add_plot1d(xy, "with lines title 'n_0 = " + ST(n[i]) + "'");
		xy.clear();
		p.m.destroy();
	}

	gp << "set xrange [1e4:1e10]\n";
	gp << "set logscale x\n";
	gp << "set xlabel \"Detector distance - km\"\n";
	gp << "set yrange [0:" << ST(1.1*longest) << "]\n";
	gp << "set ylabel \"Atmospheric Height - km\"\n";
	gp << plot;
}

void ampWave(const double& n, const double& R, const double& H, const double& N, const int& n_thread)
{
	Gnuplot gp;
	auto plot = gp.plotGroup();
	std::vector<std::pair<double, double>> xy;
	std::array<double, 3> r = {0, 0, 0};
	FlashMap *map;
	Planet3D p(R, H, N, n, 0, 0, r, "Config/map");
	//std::vector<double> l = {1e-9, 10e-9, 25e-9, 50e-9, 75e-9, 87.5e-9, 100e-9, 110e-9, 122.5e-9, 135e-9, 150e-9, 200e-9, 400e-9};
	double L = 1e4;


	for(int i = 0; i < 5; ++i){
		L *= 10;
		//for(int j = 0; j < (int)l.size(); ++j){
		for(double l = 10e-9; l < 1000e-9; l *= 1.05){
			p.updateSigma(l);
			Print("At l = " << l <<", o = " << p.o);
			map = mapThread(p, 250, 50, L, 0, n_thread);

			xy.emplace_back(l*1e9, findMax(map)/map->Int2);
			delete map;
		}
		plot.add_plot1d(xy, "with lines title 'L = " + ST(L) + "'");
		xy.clear();
	}

	for(int i = 0; i < (int)xy.size(); ++i){
		Print("x = " << xy[i].first << ", y = " << xy[i].second);
	}

	gp << "set xrange [" << std::to_string(10) << ":" << std::to_string(400) << "]\n";
	gp << "set xlabel \"Wave length - nm\"\n";
	gp << "set logscale x\n";
	gp << "set ylabel \"Max amplification\"\n";
	gp << plot;

	p.m.destroy();
}

void horPeaks(const double& n, const double& R, const double& H, const double& N, const int& n_thread)
{
	Gnuplot gp;
	auto plot = gp.plotGroup();
	std::vector<std::pair<double, double>> xy;
	std::array<double, 3> r = {0, 0, 0};
	FlashMap *map;
	//std::vector<double> obf = {0, OBF/100, OBF/10, OBF/2, OBF, 2*OBF};
	std::vector<double> obf = {0, OBF/2, OBF, 2*OBF};
	int RES = 250;
	double L = 1500000, h = 2*50.0/RES;

	for(int i = 0; i < (int)obf.size(); ++i){
		Planet3D p(R, H, N, n, obf[i], obf[i], r, "Config/map");
		map = mapThread(p, RES, 50, L, 0, n_thread);

		for(int i = 0; i <= RES; ++i){
			xy.emplace_back(i*h - 50, map->map[RES/2][i]/map->Int2);
		}

		plot.add_plot1d(xy, "with lines title 'a = " + ST(obf[i]) + "'");

		delete map;
		xy.clear();
		p.m.destroy();
	}


	gp << "set xrange [-50:50]\n";
	gp << "set xlabel \"x - km\"\n";
	gp << "set ylabel \"Amplification\"\n";
	gp << plot;
}

void peakDist(const double& n, const double& R, const double& H, const double& N, const int& n_thread)
{
	Gnuplot gp, gl;
	auto plot = gp.plotGroup(), plotl = gl.plotGroup();
	std::vector<std::pair<double, double>> xy, xyl;
	std::array<double, 3> r = {0, 0, 0};
	FlashMap *map;
	int RES = 250;
	double L = 1e5;
	Planet3D p(R, H, N, n, 0, 0, r, "Config/map");

	do{
		map = mapThread(p, RES, 50, L, 0, n_thread);
		xy.emplace_back(L, map->map[RES/2][RES/2]/map->Int2);
		delete map;
	}while((L = 10*L) < 10e11);

	gp << "set xrange [1e5:1e10]\n";
	gp << "set xlabel \"L - km\"\n";
	gp << "set logscale x\n";
	gp << "set yrange [0:250]\n";
	gp << "set ylabel \"Amplification\"\n";
	gp << "plot '-' with lines title 'Amplification vs Detector distance'\n";
	gp.send1d(xy);

	p.m.destroy();
}

void crossSec(const double& n, const double& R, const double& H, const double& N)
{
	std::array<double, 3> r = {0,0,0};
	Planet3D p(R, H, N, n, 0, 0, r, "Config/map");

	Gnuplot gp;
	std::vector<std::pair<double, double>> xy;

	for(double l = 1e-9; l < 10000e-9; l *= 1.03){
		p.updateSigma(l);
		xy.emplace_back(l*1e9, p.o);
	}

	gp << "set termoption enhanced\n";
	gp << "set xrange [1:10000]\n";
	gp << "set xlabel \"{/Symbol l} - nm\"\n";
	gp << "set logscale x\n";
	gp << "set ylabel \"{/Symbol s} - m^2\"\n";
	gp << "set logscale y\n";
	gp << "plot '-' with lines title 'Cross section'\n";
	gp.send1d(xy);

	p.m.destroy();
}

void resAmp(const double& n, const double& R, const double& H, const double& L, const double& N, const int& n_thread)
{
	std::array<double, 3> r = {0,0,0};
	Planet3D p(R, H, N, n, 0, 0, r, "Config/map");
	//std::vector<int> res = {100, 250, 500};//, 1000, 2500, 5000};
	FlashMap *map;

	Gnuplot gp;
	std::vector<std::pair<double, double>> xy;

	for(int res = 100; res < 2500; res *= 1.2){
		if(res%2) ++res;
		Print("res = " << res);
		map = mapThread(p, res, 15, L, 0, n_thread);

		xy.emplace_back((2000*15.0)/(res), findMax(map)/map->Int2);
		delete map;
	}
	//Print("1m pixel");
	//map = mapThread(p, res, 15, L, 0, n_thread);

	gp << "set termoption enhanced\n";
	//gp << "set xrange [" << ST(res[0]) << ":" << ST(res.back()) << "]\n";
	gp << "set xlabel \"Pixel size - m\"\n";
	//gp << "set logscale x\n";
	gp << "set ylabel \"Amplification\"\n";
	//gp << "set logscale y\n";
	gp << "plot '-' with lines title 'Amplification'\n";
	gp.send1d(xy);

	p.m.destroy();
}

void ampHeight(const double& n, const double& L, const double& N, const int& n_thread)
{
	std::array<double,3> r = {0,0,0};
	std::vector<double> diam = {500, 3200, 20000, 50000, 100000, 150000};
	//std::vector<double> diam = {500, 3200, 6400, 20000};
	std::vector<double> h = {2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
	int res = 500;
	FlashMap *map;

	Gnuplot gp;
	auto plot = gp.plotGroup();
	std::vector<std::pair<double, double>> xy;

	for(double D : diam){
		for(double H : h){
			Planet3D p(D/2, H, N, N_REF, 0, 0, r, "Config/map");
			map = mapThread(p, res, 20, L, 0, n_thread);
			xy.emplace_back(H, findMax(map)/map->Int2);
			delete map;
			p.m.destroy();
		}
		plot.add_plot1d(xy, "with lines title 'R = " + ST(D/2) + "'");
		xy.clear();
	}

	gp << "set xrange [0:50]\n";
	gp << "set xlabel \"H - km\"\n";
	gp << "set ylabel \"Amplification\"\n";
	gp << plot;
}







double findMax(const FlashMap* map)
{
	double max = 0;

	for(int i = 0; i < map->N+1; ++i){
		for(int j = 0; j < map->N+1; ++j){
			if(max < map->map[i][j]) max = map->map[i][j];
		}
	}

	return max;
}
