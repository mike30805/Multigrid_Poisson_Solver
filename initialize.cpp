#include "classes.h"

#include "initialize.h"



//--------------------------------------------------------------------------------
// Function    : Init_matrix
// Description : Initialize the matrix by different PROB_NUM
// Note        :
// Input       : mat : matrix
// Output      : if the initializaion is done or not.
//--------------------------------------------------------------------------------
bool Init_matrix( matrix &mat )
{
    #if ( PROB_NUM == PROB_SINWAVE )
    Init_SinWave( mat );
    #elif ( PROB_NUM == PROB_PARTICLE )
    Init_Particle( mat );
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
void Init_SinWave( matrix &mat )
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
// Description :
// Note        :
// Input       : mat : matrix
// Output      : 
//--------------------------------------------------------------------------------
void Init_Particle( matrix &mat )
{
} /// FUNCTION : Init_Particle

