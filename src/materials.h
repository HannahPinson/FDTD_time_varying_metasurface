/*
 * materials.h
 *
 *  Created on: Sep 14, 2014
 *      Author: Hannah_Pinson
 */

#ifndef MATERIALS_H_
#define MATERIALS_H_

#define _USE_MATH_DEFINES //M_PI = pi

#include <cmath>
#include <vector>


/* ----------------- functors for time variation of material parameters--------------*/


class Linear_Time_Variation{
public:
	Linear_Time_Variation(double gradient_, double offset_): gradient(gradient_), offset(offset_){;};
	double operator() (double time) {
		return (gradient*time) + offset;
	};
	double gradient, offset;
};

/*class Time_Variation_Sigma{
public:
	Time_Variation_Sigma(Linear_Time_Variation* t0_functor_, double& a_, double& total_ksi_factor_): t0_functor(t0_functor_), a(a_), total_ksi_factor(total_ksi_factor_){;};
	double operator() (double time) {
		double t0 = (*(t0_functor))(time);
		double log_factor = log(1/a);
		double sin_factor = sin(M_PI * time / t0);
		double cos_factor = cos(M_PI * time / t0);
		double exp_factor = exp(-log_factor * time / t0);
		double total = (4/t0 * exp_factor * cos_factor + (log_factor * sin_factor / M_PI))*total_ksi_factor;
		if (time>0)
			return total;
		else
			return 0;
	}
	Linear_Time_Variation* t0_functor;
	double a, total_ksi_factor;

};*/

class Sigma_m_resonances{

public:
	Sigma_m_resonances(Linear_Time_Variation* t0_functor_, double& a_, int m_, double& total_ksi_factor_, double& onset_, double& offset_): t0_functor(t0_functor_), a(a_), m(m_), total_ksi_factor(total_ksi_factor_), onset(onset_), offset(offset_){;};

	double operator() (double time, double tau) {

		double shifted_time = time-onset;

		if ((shifted_time)>0 && (time < offset)){ //heaviside theta
			double t0 = (*(t0_functor))(shifted_time);
			double sum = 0;
			/*for (int i = 0; i < m; i++){
				double sin_i = sin( (2*i + 1) * M_PI * tau / t0);
				double cos_i = cos( (2*i + 1) * M_PI * tau / t0);
				double factor_i = log(1/a) / ((2*i + 1)*M_PI);
				sum += cos_i + (factor_i * sin_i);
			}
			double total = (4/t0 * exp(-log(1/a) * tau / t0) * sum *total_ksi_factor);
			return total;*/
			double b = 1 * pow(10,-12);
			for (int i = 0; i < m; i++){
				sum += pow(-a,i) / (M_PI*b) * 1/(1 + ((tau - i*t0)/b)*(tau - i*t0)/b);
			}
			return sum;
		}
		else{
			return 0;
		}

	}
	Linear_Time_Variation* t0_functor;
	double a, total_ksi_factor;
	double onset, offset;
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

/*-------------------- material with time-dependent permittivity & permeability -------------*/

// time-varying permittivity and permeability
// no dispersion
// no conductivity
// -> e.g. 3D time-dependent linear frequency converter

// Only linear time variation (June 2014)
// Can be generalized by passing pointers to general functors

class Time_Varying_Material{
public:
	//constructor
	Time_Varying_Material(Linear_Time_Variation timeVaryingEps_, Linear_Time_Variation timeVaryingMu_, int start_node_, int end_node_):
		start_node(start_node_), end_node(end_node_), epsFunc(timeVaryingEps_),  muFunc(timeVaryingMu_)
{thickness = end_node - start_node;};
	Linear_Time_Variation epsFunc, muFunc; //functional objects -> varying permittivity and permeability
	int start_node;
	int end_node;
	int thickness;
};


/*-------------------- metasurface with dispersion --------------------------- */

class Dispersive_Metasurface{
public:
	Dispersive_Metasurface(int node_, Sigma_m_resonances sigma_functor_e_, Sigma_m_resonances sigma_functor_h_, double& eps_, double& mu_):
		node(node_), sigma_functor_e(sigma_functor_e_), sigma_functor_h(sigma_functor_h_), eps(eps_), mu(mu_){};

	int node;
	int eps;
	int mu;

	Sigma_m_resonances sigma_functor_e, sigma_functor_h;

	vector<double> saved_e;
	vector<double> saved_h;
};



#endif /* MATERIALS_H_ */
