#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "vector"

#ifdef OMP_PARALLEL
#include <omp.h>
#endif

using namespace std;
class particle
{
   double par_mass;
   double *par_pos;
   double *par_vel;
   double *par_acc;

   public:
       particle( double mass, double *pos, double *vel );
       ~particle();
}; // CLASS : particle


#endif // #ifndef PARTICLE_H
