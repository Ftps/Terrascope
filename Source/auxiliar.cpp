#include "auxiliar.hpp"


double fisqrt(double n) // Fast inverse square root algorithim for double
{
	union{
		double y;
		long i;
	} caster;
	double x2;
	const double thf = 1.5;

	x2 = n*0.5;
	caster.y  = n;                       // evil floating point bit level hacking
	caster.i  = 0x5fe6eb50c7b537a9 - ( caster.i >> 1 );       // log magic with exponents

	caster.y  = caster.y * ( thf - ( x2 * caster.y * caster.y ) );   // 1st iteration (accuracy ~1e-3)
	caster.y  = caster.y * ( thf - ( x2 * caster.y * caster.y ) );   // 2nd iteration (accuracy ~1e-6)
	//caster.y  = caster.y * ( thf - ( x2 * caster.y * caster.y ) );   // 3rd iteration, this can be removed for speed (accuracy ~1e-11)

	return caster.y;
}

double sq(double a)
{
	return a*a;
}

double gx2D(const std::function<ddd>& f, const double& x, const double& z, const double& h)
{
	return (f(x+h, z) - f(x-h, z))/(2*h);
}

double gy2D(const std::function<ddd>& f, const double& x, const double& z, const double& h)
{
	return (f(x, z+h) - f(x, z-h))/(2*h);
}

double gx3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x+h, y, z) - f(x-h, y, z))/(2*h);
}

double gy3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x, y+h, z) - f(x, y-h, z))/(2*h);
}

double gz3D(const std::function<dddd>& f, const double& x, const double& y, const double& z, const double& h)
{
	return (f(x, y, z+h) - f(x, y, z-h))/(2*h);
}

std::function<dddd> generateRefIndex(const double& n_surf, const double& H, const double& R, const double& obf_y, const double& obf_z)
{
	std::function<dddd> n;
	double betay = 1/(1 + obf_y);
	double betaz = 1/(1 + obf_z);

	n = [n_surf, H, R, betay, betaz](double x, double y, double z){
			return 1 + n_surf*exp(-(sqrt(sq(x) + sq(betay*y) + sq(betaz*z)) - R)/H);
	};

	return n;
}

FDIV generateRefFuncs(const double& n_surf, const double& H, const double& R, const double& obf_y, const double& obf_z)
{
	FDIV n;
	double betay = 1/(1 + obf_y);
	double betaz = 1/(1 + obf_z);

	n.f = [n_surf, H, R, betay, betaz](std::array<double,3> r){
			return 1 + n_surf*exp(-(sqrt(sq(r[X]) + sq(betay*r[Y]) + sq(betaz*r[Z])) - R)/H);
	};

	n.df = [n_surf, H, R, betay, betaz](std::array<double,3> r){
		double rr = sqrt(sq(r[X]) + sq(betay*r[Y]) + sq(betaz*r[Z]));
		std::array<double, 3> rl = {r[X], betay*r[Y], betaz*r[Z]};

		return (-n_surf*exp(-(rr-R)/H)/(rr*H))*rl;
	};

	return n;
}

std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& v)
{
	os << "[ ";

	for(int i = 0; i < 3; ++i){
		os << v[i];
		if(i != 2) os << ", ";
	}

	os << " ]";

	return os;
}

std::ostream& operator<<(std::ostream& os, const double A[3][3])
{
	os << "[ ";
	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			os << A[i][j];
			if(j != 2) os << ", ";
		}
		if(i != 2) os << std::endl;
	}

	os << " ]";

	return os;
}

std::array<double, 2> linReg(const std::vector<std::pair<double, double>>& xy)
{
	std::array<double, 2> res = {0,0};
	double y_avg = 0, x_avg = 0;
	double aux = 0;

	for(int i = 0; i < (int)xy.size(); ++i){
		x_avg += std::get<0>(xy[i])/xy.size();
		y_avg += std::get<1>(xy[i])/xy.size();
	}

	for(int i = 0; i < (int)xy.size(); ++i){
		aux += sq(std::get<0>(xy[i]) - x_avg);
		res[0] += (std::get<0>(xy[i]) - x_avg)*(std::get<1>(xy[i]) - y_avg);
	}

	res[0] = res[0]/aux;
	res[1] = y_avg  - res[0]*x_avg;

	return res;
}


double** createMatrix(const int& N)
{
	double** m = new double*[N];

	for(int i = 0; i < N; ++i){
		m[i] = new double[N]();
	}

	return m;
}

void freeMatrix(double const* const* m, const int N)
{
	for(int i = 0; i < N; ++i){
		delete m[i];
	}

	delete m;
}

std::array<double, 2> planeIntersect(const std::array<double, 3>& pos, const std::array<double, 3>& vel, const double& L, const double& tet, const double& h)
{
	double t, t_tet, x;

	if(!tet){
		t = (L - pos[Z])/vel[Z];
		return {pos[X] + vel[X]*t - h, pos[Y] + vel[Y]*t};
	}

	t_tet = tan(tet);
	t = (L + (pos[X] - h)*t_tet - pos[Z])/(vel[Z] - t_tet*vel[X]);
	x = pos[X] + vel[X]*t;
	if(x < 0) return {-sqrt(sq(x) + sq(L - pos[Z] - vel[Z]*t)), pos[Y] + vel[Y]*t};
	return {sqrt(sq(x - h) + sq(L - pos[Z] - vel[Z]*t)), pos[Y] + vel[Y]*t};
}
