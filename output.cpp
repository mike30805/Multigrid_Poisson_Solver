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
    const int N = mat.get_dim();
    

    FILE *file_mat;
    file_mat = fopen( filename, "a" );
    for ( int i = 0; i < N; i++ )
    {
        for ( int j = 0; j < N; j++ )
        {
            fprintf( file_mat, "%.8e ", mat.get_value( i, j ) );
        }
        fprintf( file_mat, "\n" );
    }
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
