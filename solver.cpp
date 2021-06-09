#include "classes.h"

#include "solver.h"


//--------------------------------------------------------------------------------
// Function    : Solver_Potential
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix Solver_Potential( matrix pot, matrix dens )
{
#   if ( POT_SOLVER == SOR ) 
    matrix solution = SOR_Method( pot, dens );
#   elif ( POT_SOLVER == V_CYCLE ) 
    matrix solution = V_Cycle( pot, dens );
#   elif ( POT_SOLVER == W_CYCLE ) 
    matrix solution = W_Cycle( pot, dens, 1 );
#   elif ( POT_SOLVER == FAS ) 
    // empty now
#   elif ( POT_SOLVER == FMG ) 
    matrix solution = FMG_Method( pot, dens, 2 );
#   endif // #if ( POT_SOLVER == ... )

    return solution;

} // FUNCTION : Solver_Potential



//--------------------------------------------------------------------------------
// Function    : Solver_Force
// Description : 
// Note        :
// Input       : pot  : Potential of the test problem
// Output      : Solved force 3D array
//--------------------------------------------------------------------------------
void Solver_Force( matrix pot, double ***force )
{
    const double dx = BOX_DX;
    int di, dj;

    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( d == 0 )
        {
            di = 1;
            dj = 0;
        } else if ( d == 1 )
        {
            di = 0;
            dj = 1;
        } // if ( d == 0 ) ... else if ...

        for ( int i = 0; i < BOX_N; i++ )
        {
            for ( int j = 0; j < BOX_N; j++ )
            {
                if ( (d == 0 && i == 0) || (d == 1 && j == 0) )
                {
                    force[d][i][j] = -( pot.get_value(i+di, j+dj) - pot.get_value(i, j) ) / dx;
                } else if ( (d == 0 && i == BOX_N-1) || ( d == 1 && j == BOX_N-1) )
                {
                    force[d][i][j] = -( pot.get_value(i, j) - pot.get_value(i-di, j-dj) ) / dx;
                } else
                {
                    force[d][i][j] = -0.5 * ( pot.get_value(i+di, j+dj) - pot.get_value(i-di, j-dj) ) / dx;
                } // if ( i == 0 || j == 0 ) ... else if ... else ...

            } // for ( int j = 0; j < BOX_N; j++ )
        } // for ( int i = 0; i < BOX_N; i++ )

    } // for ( int d = 0; d < N_DIMS; d++ )

} // FUNCTION : Solver_Force



//--------------------------------------------------------------------------------
// Function    : SOR_Method
// Description : 
// Note        :
// Input       : phi  : Potential of the test problem
//               dens : Density
// Output      : Solved potential 2D array
//--------------------------------------------------------------------------------
matrix SOR_Method( matrix phi, matrix dens )
{   
    phi.SOR_Exact( dens, SOR_EXACT_STEP );
    
    return phi;

} // FUNCTION : V_Cycle



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
    const int smooth_step  = SOR_SMOOTH_STEP;
    const int exact_step   = SOR_EXACT_STEP;

    const int n = dens.get_dim();

    // the smallest grid size, do the exact solver
    if ( n <= 3 )
    {
          // Exact solver
          matrix eps( phi.get_dim(), phi.get_h() );
          matrix dens2 = dens.Restriction();
          eps.SOR_smoothing( dens2, exact_step );
    }

    else{
        // Pre-Smoothing
        phi.SOR_smoothing( dens, smooth_step );
        
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
    
    phi.SOR_smoothing( dens, smooth_step );
    
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
     const int smooth_step  = SOR_SMOOTH_STEP;
     const int exact_step   = SOR_EXACT_STEP;
     const int n = dens.get_dim();
     
     // the smallest grid size, do the exact solver
     if ( n <= 3 )
     {
          // Exact solver
          matrix eps( phi.get_dim(), phi.get_h() );
          matrix dens2 = dens.Restriction();
          eps.SOR_smoothing( dens2, exact_step );
     }
     // stop recursion at smallest grid size
     else if ( n <= 7 ) // if ( n <= 3 )
     {
          // Pre-Smoothing
          phi.SOR_smoothing( dens, smooth_step );
          
          // Restriction
          matrix r = phi.Residual( dens );
          matrix rhs = r.Restriction();
          
          // Exact solver
          matrix eps( rhs.get_dim(), rhs.get_h() );
          matrix dens2 = dens.Restriction();
          eps.SOR_smoothing( dens2, exact_step );
          
          // Prolongation
          phi = phi + eps.Interpolation( phi.get_dim() );
          
          // Don't do the Post-smoothing for the left valley
          if ( LR == 1 )
          {
               // Post-Smoothing
               phi.SOR_smoothing( dens, smooth_step );
          }
     }
     else // if ( n <= 3 ) ... else if ( n <= 7 ) ...
     {
          // Pre-Smoothing
          phi.SOR_smoothing( dens, smooth_step );
          
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
               phi.SOR_smoothing( dens, smooth_step );
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
    const int smooth_step  = SOR_SMOOTH_STEP;
    const int exact_step   = SOR_EXACT_STEP;
    const int n = dens.get_dim();

    // the smallest grid size, do the exact solver
    if ( n <= 3 )
    {
          // Exact solver
          matrix eps( phi.get_dim(), phi.get_h() );
          matrix dens2 = dens.Restriction();
          eps.SOR_smoothing( dens2, exact_step );
    }

    else{
        // Pre-Smoothing
        phi.SOR_smoothing( dens, smooth_step );
        
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
        phi.SOR_smoothing( dens, smooth_step );
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
    const int smooth_step  = SOR_SMOOTH_STEP;
    const int exact_step   = SOR_EXACT_STEP;
    const int n = dens.get_dim();

    if(n<3){
        phi.SOR_smoothing( dens, exact_step );
    }
    else{
        matrix phi_coarse = phi.Restriction();
        matrix dens_coarse = dens.Restriction();
        phi_coarse = FMG_Method( phi_coarse, dens_coarse,n_fmg);
        phi = phi_coarse.Interpolation(n);
        for(int i=0;i<n_fmg;i++){
            phi = V_Cycle( phi, dens );
            //phi = W_Cycle( phi, dens, 1 );
        }
    }
    return phi;

} // FUNCTION : FMG_Method
