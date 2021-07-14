#include "fluctuation.hpp"

const std::vector<double> alpha = {5.0/3.0, 8.0/3.0};
const std::vector<double> n = {0.000253, 0.002, 0.00005};
const std::vector<double> clm = {CLM, CLM/10, CLM*10};

int main(int argc, char* argv[])
{
	int n_thread = 4;

	if(argc > 1){
		n_thread = std::atoi(argv[1]);
	}

	QApplication a(argc, argv);
	MapOpt m;
	m.l_max = 250;
	//m.alpha = 8.0/3.0;
	//m.n_ref = 0.00253;
	m.max_tet = M_PI/5;
	m.l = 2500;
	m.w = m.l*(m.max_tet/M_PI);
	for(int i = 0; i < (int)clm.size(); ++i){
		m.clm = clm[i];
		for(int j = 0; j < (int)alpha.size(); ++j){
			m.alpha = alpha[j];
			for(int k = 0; k < (int)n.size(); ++k){
				m.n_ref = n[k];
				TurbPlot map(m, n_thread);
				map.write("../Config/map_" + std::to_string(m.n_ref) + "_" + std::to_string(m.alpha) + "_" + std::to_string((5e5)*m.clm));
			}

		}
	}

	//map.show();
	//a.exec();

	return 0;
}
