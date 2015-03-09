/*
 * initialization_values.h
 *
 *  Created on: Sep 15, 2014
 *      Author: Hannah_Pinson
 */

#ifndef INITIALIZATION_VALUES_H_
#define INITIALIZATION_VALUES_H_

#include "physical_constants.h"


double all_zero (int dim_grid, int x){ //arguments unnecessary here, but needed for 'uniform' passing of function pointers
	return 0;
};

double standing_wave (int dim_grid, int x){
	// standing wave :  sin(kx)cos(wt), at t= 0 : sin(kx)*1
	// kx = 2*pi*m / N_lambda (for vacuum); m = node number, here 'x'
	return sin(2*M_PI*x / (dim_grid/2) ); // N_lambda = dim/2 -> 4th harmonic
};


#endif /* INITIALIZATION_VALUES_H_ */
