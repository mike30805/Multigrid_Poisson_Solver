#ifndef CLASSES_H
#define CLASSES_H

#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "vector"

using namespace std;

class matrix
{ 
    int dim;
    double **value;
    double h;
  
    public:
        matrix(int n,double hi);
        ~matrix();
        
        void display();
        void Error(const matrix &b);
  
        void   SOR_smoothing(const matrix& rho,double omega,int steps);
        double averaging(int i,int j);
        matrix Restriction();
        double insertion(int i,int j,int dim_in);
        matrix Interpolation(int);
        matrix Residual(const matrix& rho);
  
        void init_constant_dens();
        void init_sin_dens();
  
        double get_h();
        double get_dim();
        void   input_answer(int i,int j,double ans);

        matrix operator+(const matrix&);

}; // CLASS : matrix


#endif // #ifndef CLASSES_H
