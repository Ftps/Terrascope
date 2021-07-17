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
	int RES = 500;
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
	//std::vector<double> diam = {50000, 100000, 150000};
	//std::vector<double> diam = {500, 3200, 6400, 20000};
	std::vector<double> h = {2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
	int res = 100;
	FlashMap *map;

	Gnuplot gp;
	auto plot = gp.plotGroup();
	std::vector<std::pair<double, double>> xy;

	for(double D : diam){
		for(double H : h){
			Planet3D p(D/2, H, N, N_REF, 0, 0, r, "Config/map");
			Print("R = " << D/2 << ", H = " << H);
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

void optim(const double& n, const double& L, const double& N)
{
	std::array<double,3> r = {0,0,0};
	std::vector<double> diam = {500, 3200, 20000, 50000, 100000, 150000};
	std::vector<double> h = {2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50};
	std::array<double,3> pos, vel;
	std::array<double,2> a;
	std::array<int,2> ij;
	int res = 100;
	double S = 50, opt_M, opt_m, opt_a, hh = 2*S/res, h0;
	bool in_pic;
	std::ofstream fp("optim2.dat", std::ofstream::out);

	std::vector<std::pair<double, double>> xy1, xy2, xy3;

	fp << "Optimal data:" << std::endl;
	fp.precision(10);

	for(double D : diam){
		Gnuplot gp;
		auto plot = gp.plotGroup();
		for(double H : h){
			in_pic = false;
			opt_M = -1;
			opt_m = -1;
			fp << "D = " << D << ", H = " << H << ": ";
			Planet3D p(D/2, H, N, N_REF, 0, 0, r, "Config/map");
			Print("D = " << D << ", H = " << H << ", r_max = " << p.r_max);
			for(double rr = p.r_max; rr > p.R; rr -= 0.005){
				pos = {0, rr, -2*p.r_max};
				vel = {0, 0, 1};
				rayTrace(p, pos, vel);

				a = planeIntersect(pos, vel, L, 0, hh);
				ij[X] = (int)(res*(a[X] + S + hh))/(2*S);
				ij[Y] = (int)(res*(a[Y] + S + hh))/(2*S);
				if(ij[X] >= 0 && ij[Y] >= 0 && ij[X] < res+1 && ij[Y] < res+1){
					if(!in_pic){
						opt_M = rr;
						in_pic = true;
					}
				}
				else if(in_pic){
					opt_m = rr;
					break;
				}

			}
			h0 = 0.5*p.atm.H*(1 - p.atm.H/(2*p.R))*log(2*M_PI*sq(p.n_ref*L)/(p.R*p.atm.H)) - 9*sq(p.atm.H)/(8*p.R) + (sqrt(2)-1)*sqrt(p.atm.H*p.R/(2*M_PI))*(p.R/L);
			opt_a = h0*(1 + 2*p.atm.H/p.R) + p.R;
			p.m.destroy();
			xy1.emplace_back(H, opt_M);
			xy2.emplace_back(H, opt_m);
			xy3.emplace_back(H, opt_a);
			fp << fixed << opt_M << " " << fixed << opt_m << " " << fixed << opt_a << std::endl;
		}
		fp << std::endl;
		plot.add_plot1d(xy1, "with lines title 'max, R = " + ST(D/2) + "'");
		plot.add_plot1d(xy2, "with lines title 'min, R = " + ST(D/2) + "'");
		plot.add_plot1d(xy3, "with lines title 'calc, R = " + ST(D/2) + "'");
		gp << "set xrange [0:50]\n";
		gp << "set xlabel \"H - km\"\n";
		gp << "set ylabel \"Amplification\"\n";
		gp << plot;
		xy1.clear();
		xy2.clear();
		xy3.clear();
	}

	fp.close();
}

void atmosphericDensity(const double& n, const double& R, const double& H, const double& L, const int& n_thread)
{
	std::array<double,3> r = {0,0,0};
	FlashMap *map;

	Gnuplot gp;
	auto plot = gp.plotGroup();
	std::vector<std::pair<double, double>> xy;

	for(double N = 1e25; N < 6.5e26; N *= 2){
		Planet3D p(R, H, N, n, 0, 0, r, "Config/map");
		p.updateSigma(25e-9);
		map = mapThread(p, 250, 20, L, 0, n_thread);
		xy.emplace_back(N/1e25, findMax(map));
		delete map;
		p.m.destroy();
	}

	plot.add_plot1d(xy, "with lines title 'N vs amplification'");
	gp << "set xrange [0:50]\n";
	gp << "set xlabel \"N - molecular density x10^{25}\"\n";
	gp << "set ylabel \"Amplification\"\n";
	gp << plot;
}

void diamondSize(const double& n, const double& H, const double& L, const double& N, const int& n_thread)
{
	std::array<double,3> r = {0,0,0};
	FlashMap *map;
	//std::vector<double> H = {5, 10, 25};
	std::vector<double> rr = {500, 3250, 10000, 25000};
	std::vector<double> oo = {OBF/2, 0.002, 0.005, 0.01, 0.015, 0.02};

	Gnuplot gp;
	auto plot = gp.plotGroup();
	std::vector<std::pair<double, double>> xy;

	for(double R : rr){
		double S = 25*R/R_REF;
		xy.emplace_back(0,0);
		for(double o : oo){
			Print("R = " << R << ", o = " << o);
			Planet3D p(R, H, N, n, o, o, r, "Config/map");
			map = mapThread(p, 100, S, L, 0, n_thread);
			xy.emplace_back(o*100, 2*findDiamond(map));
			p.m.destroy();
			S *= 1.4;
		}
		plot.add_plot1d(xy, "with lines title 'R = " + ST(R) + "'");
		xy.clear();
	}


	//gp << "set xrange [1e7:1e9]\n";
	gp << "set termoption enhanced\n";
	gp << "set xlabel \"{/Symbol a} - %\"\n";
	gp << "set ylabel \"d - km\"\n";
	gp << plot;
}

void turbvs(const double& R, const double& H, const double& L, const double& N, const int& n_thread)
{
	Gnuplot gp;
	auto plot = gp.plotGroup();
	std::vector<std::pair<double, double>> xy;
	std::array<double, 3> r = {0, 0, 0};
	FlashMap *map;
	//std::vector<double> obf = {0, OBF/100, OBF/10, OBF/2, OBF, 2*OBF};
	std::vector<double> obf = {0, OBF/2, OBF, 2*OBF};
	int RES = 500;
	std::vector<double> n = {0.002, 0};
	std::vector<std::string> s = {"no turbulence", "turbulence"};
	double h = 2*50.0/RES;

	for(int k = 0; k < (int)n.size(); ++k){
		Planet3D p(R, H, N, n[k], 0, 0, r, "Config/map_0,002000_1,666667_0,100000");
		map = mapThread(p, RES, 50, L, 0, n_thread);

		for(int i = 0; i <= RES; ++i){
			xy.emplace_back(i*h - 50, map->map[RES/2][i]/map->Int2);
		}

		plot.add_plot1d(xy, "with lines title '" + s[k] +  "'");

		delete map;
		xy.clear();
		p.m.destroy();
	}

	gp << "set xrange [-50:50]\n";
	gp << "set xlabel \"x - km\"\n";
	gp << "set ylabel \"Amplification\"\n";
	gp << plot;

}

void turbsize(const double& R, const double& H, const double& L, const double& N, const int& n_thread)
{
	std::array<double, 3> r = {0,0,0};
	Planet3D p(R, H, N, 0, 0, 0, r, "Config/map_0,002000_1,666667_0,100000");
	//std::vector<int> res = {100, 250, 500};//, 1000, 2500, 5000};
	FlashMap *map;

	Gnuplot gp;
	std::vector<std::pair<double, double>> xy;

	for(int res = 50; res < 1500; res *= 1.5){
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

double findDiamond(const FlashMap* map)
{
	double max = 0, pos = 0;

	for(int i = 0; i < map->N/2 + 1; ++i){
		if(max < map->map[i][map->N/2+1]){
			max = map->map[i][map->N/2+1];
			pos = map->S*(map->N - 2*i)/map->N;
		}
	}

	return pos;
}
