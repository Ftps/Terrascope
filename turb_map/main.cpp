#include "fluctuation.hpp"

int main(int argc, char* argv[])
{
	QApplication a(argc, argv);
	MapOpt m;
	m.l_max = 150;
	m.alpha = 8.0/3.0;
	m.n_ref = 0.00253;
	m.max_tet = M_PI/5;
	m.l = 800;
	m.w = m.l*(m.max_tet/M_PI);
	TurbPlot map(m);
	map.write("../Config/map3");

	map.show();
	a.exec();

	return 0;
}
