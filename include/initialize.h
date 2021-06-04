#ifndef INITIALIZE_H
#define INITIALIZE_H


#include "macro.h"
#include "simulation_option.h"


bool Init_matrix( matrix &mat, particle *pars );

void Init_SinWave( matrix &mat, particle *pars );

void Init_Particle( matrix &mat, particle * pars );

#endif // #ifndef INITIALIZE_H
