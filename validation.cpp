#include "validate.h"


bool Validate()
{

    if ( N_DIMS != 2 )
    {
        printf( "ERROR: Simulation only support 2 dimensions now. N_DIMS = %d\n", N_DIMS );
        return false;
    }


    if ( BOX_L <= 0.0 )
    {
        printf( "ERROR: Simulation box size should be greater than zero. BOX_L = %.5f\n", BOX_L );
        return false;
    }


    if ( BOX_N <= 0 )
    {
        printf( "ERROR: Simulation resolusion should be greater than 1. BOX_N = %d\n", BOX_N );
        return false;
    }


    if ( POT_SOLVER != SOR && POT_SOLVER != V_CYCLE && POT_SOLVER != W_CYCLE && 
         POT_SOLVER != FAS && POT_SOLVER != FMG )
    {
        printf( "ERROR: Potential solver should be one of those: SOR / V_CYCLE / W_CYCLE / FAS / FMG .\n");
        return false;
    } else if ( POT_SOLVER == FAS )
    {
        printf( "ERROR: Potential solver FAS is not support yet.\n");
        return false;
    }


    if ( N_PARS < 0 )
    {
        printf( "ERROR: Particle number should NOT be negative.\n" );
        return false;
    }


    if ( MASS_TO_CELL != NGP && MASS_TO_CELL != CIC && MASS_TO_CELL != TSC )
    {
        printf( "ERROR: Particle mass deposition method should be one of those: NGP / CIC / TSC .\n" );
        return false;
    } else if ( MASS_TO_CELL == TSC )
    {
        printf( "Notice: The edge case will be solved by the CIC method.\n" );
    }


    if ( EVOLVE != EULER && EVOLVE != KDK && EVOLVE != DKD && EVOLVE != RK4 )
    {
        printf( "ERROR: Time evolve method should be one of those: EULER / KDK / DKD / RK4 .\n");
        return false;
    } else if ( EVOLVE == KDK || EVOLVE == RK4 )
    {
        printf( "ERROR: Time evolve method KDK and RK4 are not support yet.\n");
        return false;
    }

    return true;

} // FUNCTION : Validate()
