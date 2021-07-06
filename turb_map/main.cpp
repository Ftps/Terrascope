#include "fluctuation.hpp"

int main(int argc, char* argv[])
{
	int n_thread = 4;

	if(argc > 1){
		n_thread = std::atoi(argv[1]);
	}

	QApplication a(argc, argv);
	MapOpt m;
	m.l_max = 150;
	m.alpha = 8.0/3.0;
	m.n_ref = 0.00253;
	m.max_tet = M_PI/5;
	m.l = 100;
	m.w = m.l*(m.max_tet/M_PI);
	TurbPlot map(m, n_thread);
	map.write("../Config/map_test_single");

	map.show();
	a.exec();

	return 0;
}
