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

const double mu0 = 4*M_PI*1/pow(10.0,7);
const double eps0 = 8.85419 * 1/pow(10.0,12);
const double imp0 = sqrt(mu0/eps0);  // +- 377
const double c = 1/sqrt(mu0*eps0);

#endif /* PHYSICAL_CONSTANTS_H_ */
