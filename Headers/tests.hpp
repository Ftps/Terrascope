#ifndef TESTS_HPP
#define TESTS_HPP

#include "gnuplot-iostream.hpp"
#include "drawCFM.hpp"

#define N_REF 0.00253
#define R_REF 2500
#define H_REF 25
#define OBF 0.003292568
#define OBFY 0.003292568
#define OBFZ 0.003292568

void bendHeight(const std::vector<double>& n, const double& R, const double& H, const double& L, const double& N);
void detecHeight(const std::vector<double>& n, const double& R, const double& H, const double& N);
void ampWave(const double& n, const double& R, const double& H, const double& N, const int& n_thread);
void horPeaks(const double& n, const double& R, const double& H, const double& N, const int& n_thread);
void peakDist(const double& n, const double& R, const double& H, const double& N, const int& n_thread);
void crossSec(const double& n, const double& R, const double& H, const double& N);
void resAmp(const double& n, const double& R, const double& H, const double& L, const double& N, const int& n_thread);
void ampHeight(const double& n, const double& L, const double& N, const int& n_thread);

double findMax(const FlashMap* map);

#endif
