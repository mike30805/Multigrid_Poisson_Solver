#include "classes.h"
#define PI 3.14159265359


double f( double x,double y )
{
  return sin(x)*sin(y);
} // double f(double x,double y)



void solved( matrix &m )
{
  int dim = m.get_dim();
  for(int i=0; i<dim; i++)
  {
    for(int j=0; j<dim; j++)
    {
      double h = m.get_h();
      double x = h*i;
      double y = h*j;
      
      m.input_answer(i,j,f(x,y));
    } // for(int j=0; j<dim; j++)
  } // for(int i=0; i<dim; i++)
} // FUNCTION : solved



matrix V_Cycle(matrix phi, matrix dens)
{
  // Pre-Smoothing
  phi.SOR_smoothing( dens, 1.7, 3 );

  // Compute Residual Errors
  matrix r = phi.Residual( dens );

  // Restriction
  matrix rhs = r.Restriction();
  
  matrix eps( rhs.get_dim(), rhs.get_h() );
  
  // stop recursion at smallest grid size, otherwise continue recursion
  if ( rhs.get_dim() <= 5 )
  {
    matrix dens2 = dens.Restriction();
    eps.SOR_smoothing( dens2, 1.7, 100 );
  } // if ( rhs.get_dim() <= 5 )
  else
  {
    eps = V_Cycle( eps, rhs );  
  } // if ( rhs.get_dim() <= 5 ) ... else ...

  // Prolongation and Correction
  phi = phi + eps.Interpolation( phi.get_dim() );
	
  // Post-Smoothing
  phi.SOR_smoothing( dens, 1.7, 3 );
  
  return phi;
} // FUNCTION : V_Cycle


//--------------------------------------------------------------------------------
// Function    : W_Cycle
// Description : 
// Note        :
// Input       : phi  :
//               dens :
//               LR   : 0 : the left valley of the W cycle
//                      1 : the right valley of the W cycle
// Output      :
//--------------------------------------------------------------------------------
matrix W_Cycle( matrix phi, matrix dens, const int LR )
{
   const int smooth_step  = 3;
   const int exact_step   = 100;
   const double SOR_omega = 1.7;
   const int n = dens.get_dim();
   
   // the smallest grid size, do the exact solver
   if ( n <= 3 )
   {
      // Exact solver
      matrix eps( phi.get_dim(), phi.get_h() );
      matrix dens2 = dens.Restriction();
      eps.SOR_smoothing( dens2, SOR_omega, exact_step );
   }
   // stop recursion at smallest grid size
   else if ( n <= 7 ) // if ( n <= 3 )
   {
      // Pre-Smoothing
      phi.SOR_smoothing( dens, SOR_omega, smooth_step );
      
      // Restriction
      matrix r = phi.Residual( dens );
      matrix rhs = r.Restriction();
      
      // Exact solver
      matrix eps( rhs.get_dim(), rhs.get_h() );
      matrix dens2 = dens.Restriction();
      eps.SOR_smoothing( dens2, SOR_omega, exact_step );
      
      // Prolongation
      phi = phi + eps.Interpolation( phi.get_dim() );
      
      // Don't do the Post-smoothing for the right valley
      if ( LR == 1 )
      {
         // Post-Smoothing
         phi.SOR_smoothing( dens, SOR_omega, smooth_step );
      }
   }
   else // if ( n <= 3 ) ... else if ( n <= 7 ) ...
   {
      // Pre-Smoothing
      phi.SOR_smoothing( dens, SOR_omega, smooth_step );
      
      // Restriction
      matrix r = phi.Residual( dens );
      matrix rhs = r.Restriction();
      matrix eps( rhs.get_dim(), rhs.get_h() );
      
      // Left valley
      eps = W_Cycle( eps, rhs, 0 );

      // Right valley
      eps = W_Cycle( eps, rhs, 1 );
      
      // Prolongation
      phi = phi + eps.Interpolation( phi.get_dim() );
      
      // Don't do the Post-smoothing for the right valley
      if ( LR == 1 )
      {
         // Post-Smoothing
         phi.SOR_smoothing( dens, SOR_omega, smooth_step );
      }
   } // if ( n <= 3 ) ... else if ( n <= 7 ) ... else ...
   
   return phi;
   
} // FUNCTION : W_Cycle


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
  
  matrix solution2 = V_Cycle(pot,dens);
  
  solution.Error(ans);
  solution2.Error(ans);
} // FUNCTION : main 
