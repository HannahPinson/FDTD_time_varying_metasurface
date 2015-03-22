/*
 * physical_constants.h
 *
 *  Created on: Sep 15, 2014
 *      Author: Hannah_Pinson
 */

#ifndef PHYSICAL_CONSTANTS_H_
#define PHYSICAL_CONSTANTS_H_

#include <cmath>
#define _USE_MATH_DEFINES //M_PI = pi

double mu_0 = 4*M_PI*1/pow(10.0,7);
double eps_0 = 8.85 * 1/pow(10.0,12);
double imp0 = sqrt(mu_0/eps_0);
double c = 1/sqrt(mu_0*eps_0);;

#endif /* PHYSICAL_CONSTANTS_H_ */
