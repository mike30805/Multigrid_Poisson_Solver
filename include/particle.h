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

       void Par_UpdateAll( const double *vel, const double *acc, const double dt );
       void Par_AddMassToCell( double **source );

       void Par_SetPos( const double *pos );
       void Par_SetVel( const double *vel );
       void Par_SetAcc( const double *acc );
       
       void Par_GetPos( double *pos );
       void Par_GetVel( double *vel );
       void Par_GetAcc( double *acc );

}; // CLASS : particle


#endif // #ifndef PARTICLE_H
