#ifndef SOLVER_H
#define SOLVER_H

matrix V_Cycle( matrix phi, matrix dens );

matrix W_Cycle( matrix phi, matrix dens, const int LR );

matrix FAS_Method( matrix phi, matrix dens );

matrix FMG_Method( matrix phi, matrix dens, int n_fmg);

#endif // #ifndef SOLVER_H
