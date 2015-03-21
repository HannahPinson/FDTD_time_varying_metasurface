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
		return offset;
		//return (gradient*time) + offset;
	};
	double gradient, offset;
};

class Sigma_ADE{

public:
	Sigma_ADE(int& m, double& a, double total_ksi_factor, Linear_Time_Variation* t0_functor):
		_m(m), _a(a), _total_ksi_factor(total_ksi_factor), _t0_functor(t0_functor){
		vector<double> Q_initial = vector<double>(_m, 0);
		_Q_future = Q_initial;
		_Q_past =  Q_initial;
		_Q_current =  Q_initial;
	}

	/**
	 * calculate the next values of Q by solving the ADE.
	 * @param field_value : electric of magnetic field value
	 * @param delta_t
	 * @param timestep : integer or half integer value
	 */

	void calculate_next(double& field_value, double& delta_t, double& timestep){
		double t0 = (*(_t0_functor))(delta_t*timestep);
		//cout << "t0: " << t0 << endl;
		double f = 4./t0;
		double gamma = -1./t0 * log(_a);
		for (int i = 0; i < _m ; i++){

			double omega_0 = sqrt(pow( log(_a)/t0 ,2) + pow( (2*i+1)*M_PI/t0 ,2));
			double first_factor =  _total_ksi_factor * f * pow(delta_t,2) / (1 + gamma*delta_t); // f^2 ????
			double second_factor = (2 - pow(delta_t,2) * pow(omega_0,2)) / (1 + gamma*delta_t);
			double third_factor = (gamma*delta_t - 1) / (gamma*delta_t + 1) ;
			_Q_future[i] =  first_factor * field_value + second_factor * _Q_current[i] + third_factor * _Q_past[i];
			/*cout << "--------------------" << endl;
			cout << "first factor: " << first_factor << endl;
			cout << "second factor: " << second_factor << endl;
			cout << "third factor: " << third_factor << endl;
			cout << "_Q_future: " << _Q_future[i] << endl;
			cout << "_Q_current: " << _Q_current[i] << endl;
			cout << "_Q_past: " << _Q_past[i] << endl;*/

		}

	}

	double calculate_J_term(){
		double result = 0;
		for (int i = 0; i < _m ; i++){
			result += _Q_future[i] - _Q_current[i];
		}
		advance_vectors();
		return result;
	}

private:
	int _m; // number of resonances (time domain) / number of Lorentzian poles (frequency domain)
	vector <double> _Q_past, _Q_current, _Q_future;
	double _a, _total_ksi_factor;
	Linear_Time_Variation* _t0_functor;

	/**
	 */
	void advance_vectors(){
		_Q_past = _Q_current;
		_Q_current = _Q_future;
	}

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
	Dispersive_Metasurface(int node_, Sigma_ADE sigma_functor_e_, Sigma_ADE sigma_functor_h_, double& eps_, double& mu_):
		node(node_), sigma_functor_e(sigma_functor_e_), sigma_functor_h(sigma_functor_h_), eps(eps_), mu(mu_){};

	int node;
	int eps;
	int mu;

	Sigma_ADE sigma_functor_e, sigma_functor_h;


	/*vector<double> saved_e;
	vector<double> saved_h;*/
};



#endif /* MATERIALS_H_ */
