#include "macro.h"
#include "validate.h"
#include "classes.h"
#include "particle.h"
#include "solver.h"
#include "simulation_option.h"

double f( const double x, const double y )
{   
    return sin(x) * sin(y) + background_pot;
} // FUNCTION : f



void solved( matrix &m )
{
    const int dim = m.get_dim();
    const int dim_2 = dim*dim;
    
    for ( int idx = 0; idx < dim_2; idx++ )
    {
        const int i = idx / dim;
        const int j = idx % dim;
        
        const double h = m.get_h();
        const double x = h*i;
        const double y = h*j;
      
        m.input_answer( i, j, f(x, y) );
    } // for ( int idx = 0; idx < dim_2; idx++ )

} // FUNCTION : solved



int main()
{
    if ( not Validate() ) return 0;

    matrix pot( BOX_N, BOX_DX );
    pot.init_potential();
    
    matrix dens( BOX_N, BOX_DX );
    dens.init_density();
    
    matrix ans( BOX_N, BOX_DX );
    solved(ans);
    
    auto start = chrono::steady_clock::now();

    // solve potential
#   if ( POT_SOLVER == SOR ) 
    // empty now
#   elif ( POT_SOLVER == V_CYCLE ) 
    matrix solution = V_Cycle( pot, dens );
#   elif ( POT_SOLVER == W_CYCLE ) 
    matrix solution = W_Cycle( pot, dens, 1 );
#   elif ( POT_SOLVER == FAS ) 
    // empty now
#   elif ( POT_SOLVER == FMG ) 
    matrix solution = FMG_Method( pot, dens, 2 );
#   endif // #if ( POT_SOLVER == ... )

    solution.Error( ans ); // print the error

    auto elapsed = chrono::steady_clock::now() - start;
    auto sec_double = chrono::duration<double>(elapsed);     // double
    cout << sec_double.count() << endl;

} // FUNCTION : main 
