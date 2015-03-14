/*
 * sources.h
 *
 *  Created on: Sep 15, 2014
 *      Author: Hannah_Pinson
 */

#ifndef SOURCES_H_
#define SOURCES_H_

#include <cmath>

#define _USE_MATH_DEFINES //M_PI = pi
#define SQR(X) ((X)*(X))

#include "physical_constants.h"

using namespace std;


/* ----------- Source functors ----------------- */

// functors (= functional objects, objects that can act as functions by overloading the "( )"-operator) that take timestep as an argument

// used as members of actual source objects (implemented below)

class Source_functor{
public:
  virtual double operator()(int& timestep) = 0; // "blueprint" -> all derived classes should implement ()-operator
};


//derived classes, overloading the ()-operator

class Gaussian_functor : public Source_functor {
public:
    Gaussian_functor(double delay_, double width_, double scale_factor_): delay(delay_), width(width_), scale_factor(scale_factor_){;};
    double operator() (int& timestep) {return  scale_factor * exp(- SQR( (timestep - delay) / width) ); };
    double delay, width, scale_factor;
};

class Sine_functor : public Source_functor{
public:
    Sine_functor(double S_c_, double N_lambda_, double scale_factor_): S_c(S_c_), N_lambda(N_lambda_), scale_factor(scale_factor_) {;};
    double operator() (int& timestep) {return  scale_factor * sin(2 * M_PI *  S_c  / N_lambda * timestep); };
    double S_c, N_lambda, scale_factor;
};

class Gaussian_packet_functor : public Source_functor{
public:
    Gaussian_packet_functor(double S_c_, double N_lambda_, double delay_, double width_, double scale_factor_): S_c(S_c_), N_lambda(N_lambda_), delay(delay_), width(width_), scale_factor(scale_factor_){;};
    double operator() (int& timestep) {return scale_factor * sin(2 * M_PI *  S_c  / N_lambda * (timestep-delay)) * exp(- SQR((timestep - delay) / width) ); };
    double S_c, N_lambda, delay, width, scale_factor;
};


class No_source_functor : public Source_functor{
public:
    //default constructor
    double operator() (int& time){return 0;};
};


/* ----------------- Source  ----------------------*/


class Source{
public:
    Source(int node_, char type_, Source_functor* ptr_source_functor_): node(node_), type_of_source(type_), ptr_source_functor(ptr_source_functor_){;};
    int node;
    char type_of_source; // 'a' = additive or 'h' = hard or 's' = subtractive
    Source_functor* ptr_source_functor; // pointer to source functor, all objects of derived classes (e.g. Gaussian_packet_functor) also allowed
};

//standard no source
No_source_functor no_source_functor;
Source standard_no_source(0, 'a', &no_source_functor);


#endif /* SOURCES_H_ */
