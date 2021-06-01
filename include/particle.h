#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "vector"
#include "macro.h"
#include "simulation_option.h"

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

       void update_all( const double *vel, const double *acc, const double dt );

       void set_pos( const double *pos );
       void set_vel( const double *vel );
       void set_acc( const double *acc );
       
       void get_pos( double *pos );
       void get_vel( double *vel );
       void get_acc( double *acc );
}; // CLASS : particle


#endif // #ifndef PARTICLE_H
