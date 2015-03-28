/*
 * materials.h
 *
 *  Created on: Sep 14, 2014
 *      Author: Hannah_Pinson
 */

#ifndef MATERIALS_H_
#define MATERIALS_H_

#include <vector>
#include <cmath>

#define _USE_MATH_DEFINES //M_PI = pi

using namespace std;


/* ----------------- functors for time variation of material parameters--------------*/


class Linear_Time_Variation{
public:
	Linear_Time_Variation(double gradient_, double offset_): gradient(gradient_), offset(offset_){};
	double operator() (double time) {return (gradient*time) + offset;};
	double gradient, offset;
};


class Sigma_m_resonances{

public:
	Sigma_m_resonances(Linear_Time_Variation* t0_functor_, double& a_, int m_, double& total_ksi_factor_): t0_functor(t0_functor_), a(a_), m(m_), total_ksi_factor(total_ksi_factor_){};
	double operator() (double time, double tau);
	Linear_Time_Variation* t0_functor;
	double a, total_ksi_factor;
	int m;
};


/*------------------ static material ----------------- */

//no time dependence
//no dispersion
//electric and magnetic conductivity
//-> e.g. PML


class Static_Material{
public:
	//constructor
	Static_Material(double& sigma_e_, double& sigma_h_, double& eps_, double& mu_, int start_node_, int end_node_):
		sigma_e(sigma_e_), sigma_h(sigma_h_),  eps(eps_), mu(mu_), start_node(start_node_), end_node(end_node_)
{thickness = end_node - start_node;};

	double sigma_e, sigma_h; //electric and magnetic conductivity
	double eps, mu; //permittivity, permeability
	int start_node;
	int end_node;
	int thickness;
};



/*-------------------- metasurface with dispersion --------------------------- */

class Dispersive_Metasurface{
public:
	Dispersive_Metasurface(int node_, Sigma_m_resonances sigma_functor_e_, Sigma_m_resonances sigma_functor_h_);

	int node;
	Sigma_m_resonances sigma_functor_e, sigma_functor_h;
	vector<double> saved_e;
	vector<double> saved_h;
};



#endif /* MATERIALS_H_ */
