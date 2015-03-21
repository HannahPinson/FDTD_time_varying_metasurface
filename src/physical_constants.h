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

double c = 299792458;
double imp0 = 376.73; //impedance of free space = sqrt(mu_0/epsilon_0)
double mu_0 = 4*M_PI*1/pow(10.0,7);
double eps_0 = 8.854 * 1/pow(10.0,12);



#endif /* PHYSICAL_CONSTANTS_H_ */
