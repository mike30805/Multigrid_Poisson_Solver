#include "classes.h"
#define PI 3.14159265359

#include "solver.h"

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
    int N = 2000;
    double h = PI/(N-1);
    matrix pot( N, h );
    pot.init_potential();
    
    matrix dens( N, h );
    dens.init_density();
    
    matrix ans( N, h );
    solved(ans);
    
    /*matrix solution = V_Cycle( pot, dens );
    pot.init_potential();
    matrix solution2 = W_Cycle( pot, dens, 1 );
    pot.init_potential();
        
    solution.Error( ans );
    solution2.Error( ans );*/

    
    auto start = chrono::steady_clock::now();

    // V,W-Cycle
    matrix solution = V_Cycle( pot, dens );
    pot.init_potential();
    matrix solution2 = W_Cycle( pot, dens, 1 );
    pot.init_potential();
        
    solution.Error( ans );
    solution2.Error( ans );
    
    //FMG Method
    matrix solution3 = FMG_Method( pot, dens,2);
    
    pot.init_potential();

    matrix ans3(solution3.get_dim(),solution3.get_h());
    solved(ans3);
    solution3.Error(ans3);

    auto elapsed = chrono::steady_clock::now() - start;
    auto sec_double = chrono::duration<double>(elapsed);     // double
    cout<<sec_double.count()<<endl;
} // FUNCTION : main 
