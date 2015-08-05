#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "specrel.h"

#ifndef _REGGAE_H
#define _REGGAE_H


#define SortLength 7
#define DGPI 3.141592653589793238

// GENBOD
// P is total four-momentum
// n is the number of particles 
// *mass is an array to the particle masses
// *op is the array of particle momenta 
// *seed is the pointer to the random generator seed
double Mconserv (vector4 P, int n,double *mass, vector4 *op,long int *seed);

// rescattering routine
// n is the number of particles
// *avec is the array with four-momenta of the particles 
// *seed is the seed for the random number generator
double collision (int n,vector4 *avec, long int *seed);



// here is the random generator from numerical recipes
// first the parameters for the generator

#define KAS_RND_IM1 2147483563
#define KAS_RND_IM2 2147483399
#define KAS_RND_AM (1.0/KAS_RND_IM1)
#define KAS_RND_IMM1 (KAS_RND_IM1-1)
#define KAS_RND_IA1 40014
#define KAS_RND_IA2 40692
#define KAS_RND_IQ1 53668
#define KAS_RND_IQ2 52774
#define KAS_RND_IR1 12211
#define KAS_RND_IR2 3791
#define KAS_RND_NTAB 32
#define KAS_RND_NDIV (1+KAS_RND_IMM1/KAS_RND_NTAB)
#define KAS_RND_EPS 1.2e-7
#define KAS_RND_RNMX (1.0-KAS_RND_EPS)
double KAS_rndm(long * );

//an implementation of the quicksort algorithm
void quicksort(int , double a[]);

#endif
