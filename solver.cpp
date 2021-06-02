#include "classes.h"

#include "solver.h"



//--------------------------------------------------------------------------------
// Function    : V_Cycle
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix V_Cycle( matrix phi, matrix dens )
{   
    const int smooth_step  = 3;
    const int exact_step   = 100;
    const double SOR_omega = 1.9;
    const int n = dens.get_dim();

    // the smallest grid size, do the exact solver
    if ( n <= 3 )
    {
          // Exact solver
          matrix eps( phi.get_dim(), phi.get_h() );
          matrix dens2 = dens.Restriction();
          eps.SOR_smoothing( dens2, SOR_omega, exact_step );
    }

    else{
        // Pre-Smoothing
        phi.SOR_smoothing( dens, SOR_omega, smooth_step );
        
        // Compute Residual Errors
        matrix r = phi.Residual( dens );
    
        // Restriction
        matrix rhs = r.Restriction();
        
        matrix eps( rhs.get_dim(), rhs.get_h() );
        
        // Go to smaller grid size
        eps = V_Cycle( eps, rhs );  
    
        // Prolongation and Correction
        phi = phi + eps.Interpolation( phi.get_dim() );
        
        // Post-Smoothing
    }
    
    phi.SOR_smoothing( dens, SOR_omega, smooth_step );
    
    return phi;

} // FUNCTION : V_Cycle



//--------------------------------------------------------------------------------
// Function    : W_Cycle
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
//               LR   : 0 : the left valley of the W cycle
//                      1 : the right valley of the W cycle
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix W_Cycle( matrix phi, matrix dens, const int LR )
{
     const int smooth_step  = 3;
     const int exact_step   = 100;
     const double SOR_omega = 1.9;
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
          
          // Don't do the Post-smoothing for the left valley
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
          matrix r   = phi.Residual( dens );
          matrix rhs = r.Restriction();
          matrix eps( rhs.get_dim(), rhs.get_h() );
          
          // Left valley
          eps = W_Cycle( eps, rhs, 0 );
  
          // Right valley
          eps = W_Cycle( eps, rhs, 1 );
          
          // Prolongation
          phi = phi + eps.Interpolation( phi.get_dim() );
          
          // Don't do the Post-smoothing for the left valley
          if ( LR == 1 )
          {
               // Post-Smoothing
               phi.SOR_smoothing( dens, SOR_omega, smooth_step );
          }
     } // if ( n <= 3 ) ... else if ( n <= 7 ) ... else ...
     
     return phi;
   
} // FUNCTION : W_Cycle



//--------------------------------------------------------------------------------
// Function    : FAS_Method
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix FAS_Method( matrix phi, matrix dens )
{   
    const int smooth_step  = 3;
    const int exact_step   = 100;
    const double SOR_omega = 1.9;
    const int n = dens.get_dim();

    // the smallest grid size, do the exact solver
    if ( n <= 3 )
    {
          // Exact solver
          matrix eps( phi.get_dim(), phi.get_h() );
          matrix dens2 = dens.Restriction();
          eps.SOR_smoothing( dens2, SOR_omega, exact_step );
    }

    else{
        // Pre-Smoothing
        phi.SOR_smoothing( dens, SOR_omega, smooth_step );
        
        // Compute Residual Errors
        matrix res = phi.Residual( dens );
    
        // Restriction
        matrix res_coarse = res.Restriction();
        
        matrix phi_coarse = phi.Restriction();

        matrix dens_coarse = phi_coarse.Laplacian() +res_coarse;
        
        // Go to smaller grid size
        phi_coarse = V_Cycle( phi_coarse, dens_coarse );  
    
        // Prolongation and Correction
        matrix corr = phi_coarse -phi.Restriction();
        phi = phi + corr.Interpolation( phi.get_dim() );
        
        // Post-Smoothing
        phi.SOR_smoothing( dens, SOR_omega, smooth_step );
    }
    
    return phi;

} // FUNCTION : FAS_Method



//--------------------------------------------------------------------------------
// Function    : FMG_Method
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
//               n_fmg: On each working level, one applies n_fmg MG cycles
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix FMG_Method( matrix phi, matrix dens, int n_fmg)
{   
    const int smooth_step  = 3;
    const int exact_step   = 100;
    const double SOR_omega = 1.9;
    const int n = dens.get_dim();

    if(n<3){
        phi.SOR_smoothing( dens, SOR_omega, exact_step );
    }
    else{
        matrix phi_coarse = phi.Restriction();
        matrix dens_coarse = dens.Restriction();
        phi_coarse = FMG_Method( phi_coarse, dens_coarse,n_fmg);
        phi = phi_coarse.Interpolation(n);
        for(int i=0;i<n_fmg;i++){
            phi = V_Cycle( phi, dens );
        }
    }
    return phi;

} // FUNCTION : FMG_Method
