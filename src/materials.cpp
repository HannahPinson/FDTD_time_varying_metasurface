/*
 * materials.cpp
 *
 *  Created on: Mar 27, 2015
 *      Author: Hannah_Pinson
 */


#include "materials.h"


using namespace std;


double Sigma_m_resonances::operator() (double time, double tau){
	double t0 = (*(t0_functor))(time);
	double sum = 0;
	for (int i = m; i < m + 1; i++){ //!!!!!!
		double sin_i = sin( (2*i + 1) * M_PI * tau / t0);
		double cos_i = cos( (2*i + 1) * M_PI * tau / t0);
		double factor_i = log(1/a) / ((2*i + 1)*M_PI);
		sum += cos_i + (factor_i * sin_i);
	}
	double total = (4/t0 * exp(-log(1/a) * tau / t0) * sum *total_ksi_factor);
	return total;
};


/*-------------------- metasurface with dispersion --------------------------- */


Dispersive_Metasurface::Dispersive_Metasurface(int node_, Sigma_m_resonances sigma_functor_e_, Sigma_m_resonances sigma_functor_h_):
				node(node_), sigma_functor_e(sigma_functor_e_), sigma_functor_h(sigma_functor_h_){
	saved_e.push_back(0);
	saved_h.push_back(0);
};



