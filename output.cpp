#include "classes.h"
#include "particle.h"

#include "output.h"



//--------------------------------------------------------------------------------
// function    : Output_matrix
// description : 
// note        :
// input       : mat      : matrix will be outputed
//               filename : save file name
// output      : 
//--------------------------------------------------------------------------------
void Output_matrix( matrix mat, const char filename[] )
{
    int n = mat.get_dim();
    #if ( N_DIMS == 2 )
    int cells = n*n;
    #elif ( N_DIMS == 3)
    int cells = n*n*n;
    #endif

    FILE *file_mat;
    file_mat = fopen( filename, "a" );
    
    for ( int idx = 0; idx < cells; idx++ )
    {
        const int i = idx%n;
        const int j = ( idx%(n*n) ) / n;
        const int k = idx/(n*n);

        fprintf( file_mat, "%.8e ", mat.get_value( idx ) );

        if ( i == n-1 && j == n-1 && k == n-1 )    fprintf( file_mat, "\n" );
        if ( i == n-1 && j == n-1 )                fprintf( file_mat, "\n" );
        if ( i == n-1 )                            fprintf( file_mat, "\n" );

    } // for ( int idx = 0; idx < cells; idx++ )

    fclose( file_mat );

} // FUNCTION : Output_matrix



//--------------------------------------------------------------------------------
// function    : Output_particles
// description : 
// note        :
// input       : mat      : matrix will be outputed
//               filename : save file name
// output      : 
//--------------------------------------------------------------------------------
void Output_particles( particle *pars, const char filename[] )
{
    double mass;
    double *pos = new double[N_DIMS];
    double *vel = new double[N_DIMS];
    FILE *file_par;

    file_par = fopen( filename, "a" );
    
    // header
    fprintf( file_par, "# idx           mass           " );
    for ( int d = 0; d < N_DIMS; d++ )    fprintf( file_par, "pos%d           ", d );
    for ( int d = 0; d < N_DIMS; d++ )    fprintf( file_par, "vel%d           ", d );
    fprintf( file_par, "\n" );


    for ( int p = 0; p < N_PARS; p++ )
    {
        pars[p].Par_GetMass( mass );
        pars[p].Par_GetPos( pos ); 
        pars[p].Par_GetVel( vel );
        
        fprintf( file_par, "%5d %.8e ", p, mass );
        for ( int d = 0; d < N_DIMS; d++ )    fprintf( file_par, "%.8e ", pos[d] );
        for ( int d = 0; d < N_DIMS; d++ )    fprintf( file_par, "%.8e ", vel[d] );
        fprintf( file_par, "\n" );
    }

    fclose( file_par );
    
    delete[] pos;
    delete[] vel;

} // FUNCTION : Output_particles



//--------------------------------------------------------------------------------
// function    : Output_parameters
// description : 
// note        :
// input       : 
// output      : 
//--------------------------------------------------------------------------------
void Output_parameter()
{
    FILE *file_mat;
    file_mat = fopen( "Simulation_Parameters.txt", "a" );
    
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "Simulation box\n");
    fprintf( file_mat, "N_DIMS            %-7d\n",   N_DIMS );
    fprintf( file_mat, "BOX_N             %-7d\n",   BOX_N );
    fprintf( file_mat, "BOX_L             %-7.5f\n", BOX_L );
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "\n\n" );
    
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "Simulation time\n");
    fprintf( file_mat, "END_TIME          %-7.5f\n", END_TIME );
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "\n\n" );
    
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "Problem number\n");
    fprintf( file_mat, "PROB_NUM          %-3d\n",   PROB_NUM );
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "\n\n" );
    
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "Potential solver\n");
    fprintf( file_mat, "BG_POTENTIAL      %-7.5f\n", BG_POTENTIAL );
    fprintf( file_mat, "POT_SOLVER        %-3d\n",   POT_SOLVER );
    fprintf( file_mat, "SOR_OMEGA         %-7.5f\n", SOR_OMEGA );
    fprintf( file_mat, "SOR_SMOOTH_STEP   %-7d\n",   SOR_SMOOTH_STEP );
    fprintf( file_mat, "SOR_EXACT_STEP    %-7d\n",   SOR_EXACT_STEP );
    fprintf( file_mat, "SOR_ERROR         %-1.7e\n", SOR_ERROR );
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "\n\n" );
    
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "Parallel\n");
    #ifdef OMP_PARALLEL
    fprintf( file_mat, "OMP_PARALLEL      true\n");
    #else
    fprintf( file_mat, "OMP_PARALLEL      false\n");
    #endif
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "\n\n" );
    
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "Particle\n");
    fprintf( file_mat, "N_PARS            %-7d\n",   N_PARS );
    fprintf( file_mat, "MASS_TO_CELL      %-2d\n",   MASS_TO_CELL );
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "\n\n" );
    
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "Time evolution\n");
    fprintf( file_mat, "EVOLVE            %-2d\n",   EVOLVE );
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "\n\n" );
    
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "Output\n");
    fprintf( file_mat, "OUT_DT            %-7.5f\n", OUT_DT );
    fprintf( file_mat, "==========================================\n" );
    fprintf( file_mat, "\n\n" );
    
    fclose( file_mat );

} // FUNCTION : Output_parameter
