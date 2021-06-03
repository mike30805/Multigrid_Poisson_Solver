#include "classes.h"

#include "output.h"



//--------------------------------------------------------------------------------
// Function    : Output_matrix
// Description : 
// Note        :
// Input       : mat      : matrix will be outputed
//               filename : save file name
// Output      : 
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
