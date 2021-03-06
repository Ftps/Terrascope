#include "sphere_harm.hpp"

const std::complex<double> i_c(0.0, 1.0);

long double factorial(const int& a)
{
	return ((a == 1) || (a == 0)) ? 1.0 : a*factorial(a-1);
}

long double coeffFact(const int& l, const int& m)
{
	long double res = 1;
	for(int i = m; i > -m; --i){
		res = res*(l+i);
	}

	return 1/res;
}

long double binTop(long double a, const int& k)
{
	long double l = 1;

	for(int i = 0; i < k; ++i){
		l *= a--;
	}

	return l;
}

long double ALpoly(const long double& x, const int& m, const int& l, long double& prev_LP)
{
	long double res, aux;

	if(abs(m) > l) return 0;
	else if (m == l){
		prev_LP = 0;
		return ((l%2) ? -1 : 1)*powl(2, l)*powl(1 - x*x, l/2.0)*binTop(l-0.5, l);
	}
	else if (m < 0){
		return ((m%2) ? -1 : 1)*factorial(l+m)/factorial(l-m)*ALpoly(x, -m, l, prev_LP);
	}
	else{
		aux = ALpoly(x, m, l-1, prev_LP);
		res = ((2*l-1)*x*aux - (l+m-1)*prev_LP)/(l - m);
		prev_LP = aux;
		return res;
	}
}

long double sphereHarm(const long double& tet, const long double& phi, const int& m, const int& l)
{
	long double prev_LP;
	return cos(m*phi)*ALpoly(cos(tet), m, l, prev_LP);
}

long double sphereHarmIm(const long double& tet, const long double& phi, const int& m, const int& l)
{
	long double prev_LP;
	return sin(m*phi)*ALpoly(cos(tet), m, l, prev_LP);
}

std::vector<std::complex<double>> ditfft(const double* x, const int& N, const int& s)
{
	std::vector<std::complex<double>> res, aux1, aux2;
	std::complex<double> p(x[0],0), q;

	if(N == 1){
		res.push_back(p);
		return res;
	}
	aux1 = ditfft(x, N/2, 2*s);
	aux2 = ditfft(x+s, N/2, 2*s);

	res.resize(N);
	for(int k = 0; k < N/2; ++k){
		p = aux1[k];
		q = aux2[k]*exp(-2*M_PI*i_c*(((double)k)/N));
		res[k] = p + q;
		res[k + N/2] = p - q;
	}

	return res;
}

double maxfreq(const int& m, const int& l)
{
	int N = 2*32768;
	double *y = new double[N];
	long double coeff, aux, dx = 2*M_PI/N;
	double max = 0;
	int i_max = 0;
	std::vector<std::complex<double>> dft;

	aux = sqrtl(coeffFact(l, m));
	coeff = sqrtl(((2.0*l + 1.0)/(4.0*M_PI)))*aux;

	for(int i = 0; i < N; ++i){
		y[i] = coeff*ALpoly(cos(i*dx-M_PI/2), m, l, aux);
	}

	dft = ditfft(y, N);

	for(int i = 0; i < N/2; ++i){
		if(std::abs(dft[i]) > max){
			max = std::abs(dft[i]);
			i_max = i;
		}
	}

	delete y;

	return i_max/(2*M_PI);
}










void fillNumbers(const int& l_max, const int& N_min, const int& N_max)
{
	std::ofstream fp, fp_coeffs;
	int N;
	long double dx, prev_LP, aux;
	const long double alpha = (N_max - N_min)/(long double)(l_max - 1.0);

	fp_coeffs.open(fileloc + "Coeffs.dat", std::ofstream::out | std::ofstream::trunc);
	fp_coeffs << l_max << '\n';

	for(int l = 1; l <= l_max; ++l){
		fp.open(fileloc + "P" + std::to_string(l) + ".dat", std::ofstream::out | std::ofstream::trunc);
		N = (int)(N_min + alpha*(l-1));
		dx = 2.0/(N-1);
		fp << N << '\n';
		std::cout << "l = " << l << std::endl;
		for(int m = 0; m <= l; ++m){
			aux = sqrtl(coeffFact(l, m));
			aux = sqrtl(((2.0*l + 1.0)/(4.0*M_PI)))*aux;
			fp_coeffs << aux << " ";
			for(long double x = -1; x <= 1; x += dx){
				fp << aux*ALpoly(x, m, l, prev_LP) << " ";
			}
			fp << '\n';
		}
		fp.close();
		fp_coeffs << '\n';
	}
	fp_coeffs.close();
}

void verifyCoeff()
{
	std::ifstream fp(fileloc + "Coeffs.dat", std::ifstream::in);
	int l_max;
	long double aux, N;
	bool zero;

	fp >> l_max;

	std::cout << "Checking coeffs . . ." << std::endl;

	for(int l = 1; l <= l_max; ++l){
		for(int m = 0; m <= l; ++m){
			fp >> aux;
			if(aux == 0){
				std::cout << "Error found in l = " << l << " and m = " << m << std::endl;
			}
		}
	}

	fp.close();
	std::cout << "Checking poly values . . ." << std::endl;
	for(int l = 1; l <= l_max; ++l){
		fp.open(fileloc + "P" + std::to_string(l) + ".dat", std::ofstream::in);
		fp >> N;
		if(!(l%50)) std::cout << "l = " << l << std::endl;
		for(int m = 0; m <= l; ++m){
			zero = false;
			for(int i = 0; i < N; ++i){
				fp >> aux;
				if(aux == 0){
					if(zero == true){
						std::cout << "Error found in l = " << l << " and m = " << m << std::endl;
						std::cout << "Two zeros." << std::endl;
					}
					zero = true;
					continue;
				}
				else if(std::isinf(aux) || std::isnan(aux)){
					std::cout << "Error found in l = " << l << " and m = " << m << std::endl;
					std::cout << "Out of bounds value." << std::endl;
				}
				zero = false;
			}

		}
		fp.close();
	}

	std::cout << "Done." << std::endl;
}
