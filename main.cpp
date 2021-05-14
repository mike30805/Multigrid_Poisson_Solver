#include "classes.h"
#define PI 3.14159265359
double f(double x,double y){
  return sin(x)*sin(y);
}
void solved(matrix &m){
  int dim = m.get_dim();
  for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            double h = m.get_h();
            double x = h*i;
            double y = h*j;
            
            m.input_answer(i,j,f(x,y));
        }
    }
}

matrix V_Cycle(matrix phi,matrix dens){

	// Pre-Smoothing
	phi.SOR_smoothing(dens,1.7,3);
	
	// Compute Residual Errors
	matrix r = phi.Residual(dens);
	
	// Restriction
	matrix rhs = r.Restriction();

	matrix eps(rhs.get_dim(),rhs.get_h());
  
	// stop recursion at smallest grid size, otherwise continue recursion
	if (rhs.get_dim()<=5){
    matrix dens2 = dens.Restriction();
    eps.SOR_smoothing(dens2,1.7,100);
  }
        	
	else {
    eps = V_Cycle(eps,rhs);  
  }       

	// Prolongation and Correction
	phi = phi + eps.Interpolation(phi.get_dim());
	
  // Post-Smoothing
	phi.SOR_smoothing(dens,1.7,3); 
  return phi;
}

int main()
{
  int N = 8000;
  double h = PI/(N-1);
  matrix pot(N,h);

  matrix dens(N,h);
  dens.init_sin_dens();
  
  matrix ans(N,h);
  solved(ans);
  
  matrix solution = V_Cycle(pot,dens);
  
  solution.Error(ans);
}
  
