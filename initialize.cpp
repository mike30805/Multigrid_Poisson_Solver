#include "classes.h"
#include "particle.h"

#include "initialize.h"



//--------------------------------------------------------------------------------
// Function    : Init_matrix
// Description : Initialize the matrix by different PROB_NUM
// Note        :
// Input       : mat  : matrix
//               pars : particles
// Output      : if the initializaion is done or not.
//--------------------------------------------------------------------------------
bool Init_matrix( matrix &mat, particle *pars )
{
    #if ( PROB_NUM == PROB_SINWAVE )
    Init_SinWave( mat, pars );
    #elif ( PROB_NUM == PROB_PARTICLE )
    Init_Particle( mat, pars );
    #else
    printf("No such a problem number %d.\n", PROB_NUM);
    return false;
    #endif

    return true;

} // FUNCTION : Init_matrix( matrix &mat )



//--------------------------------------------------------------------------------
// Function    : Init_SinWave
// Description : Matrix value set as f(x, y) = -2 * sin(x) * sin(y) 
// Note        :
// Input       : mat : matrix
// Output      : 
//--------------------------------------------------------------------------------
void Init_SinWave( matrix &mat, particle *pars )
{
    for( int i = 0; i < BOX_N; i++ )
    {
        for( int j = 0; j < BOX_N; j++ )
        {
            const double x = BOX_DX*i;
            const double y = BOX_DX*j;

            mat.set_value( i, j, -2.*sin(x)*sin(y) );
        } // for( int j = 0; j < dim; j++ )
    } // for( int i = 0; i < dim; i++ )

} // FUNCTION : Init_SinWave



//--------------------------------------------------------------------------------
// Function    : Init_Particle
// Description : Problem number 2
// Note        :
// Input       : mat : matrix
// Output      : 
//--------------------------------------------------------------------------------
void Init_Particle( matrix &mat, particle *pars )
{
    for( int i = 0; i < BOX_N; i++ )
    {
        for( int j = 0; j < BOX_N; j++ )
        {
            mat.set_value( i, j, 0.0 );
        } // for( int j = 0; j < dim; j++ )
    } // for( int i = 0; i < dim; i++ )

    // particles
    double *pos = new double[N_DIMS];
    double *vel = new double[N_DIMS];
    for ( int p = 0; p < N_PARS; p++ )
    {
        for( int d = 0; d < N_DIMS; d++ )
        {
            pos[d] = p * PI / N_PARS+0.03;
            vel[d] = 0.0;
        }

        pars[p].Par_SetMass( 10.0 );
        pars[p].Par_SetPos( pos );
        pars[p].Par_SetVel( vel );
        pars[p].Par_AddMassToCell( mat );

    } // for ( int p = 0; p < N_PARS; p++ )

    delete[] pos;
    delete[] vel;


} /// FUNCTION : Init_Particle



