#include "classes.h"
#include "particle.h"


particle::particle()
{
    par_pos = new double[N_DIMS];
    par_vel = new double[N_DIMS];
    par_acc = new double[N_DIMS];

    par_mass   = 0.0;
    for ( int d = 0; d < N_DIMS; d++ )
    {
        par_pos[d] = 0.0;
        par_vel[d] = 0.0;
    } // for ( int d = 0; d < N_DIMS; d++ )

} // CONSTRUCTER : particle::particle()



particle::particle( double mass, double *pos, double *vel )
{
    par_pos = new double[N_DIMS];
    par_vel = new double[N_DIMS];
    par_acc = new double[N_DIMS];

    par_mass   = mass;
    for ( int d = 0; d < N_DIMS; d++ )
    {
        par_pos[d] = pos[d];
        par_vel[d] = vel[d];
    } // for ( int d = 0; d < N_DIMS; d++ )

} // CONSTRUCTER : particle::particle( double mass, double *pos, double *vel )



particle::~particle()
{
} // DESTRUCTER : particle::~particle()



bool particle::Par_InBox()
{
    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( this->par_pos[d] < 0.0 || BOX_L < this->par_pos[d] )    return false;
    }
    return true;
 
} // FUNCTION : particle::Par_InBox()



void particle::Par_AddMassToCell( matrix &source )
{
    // 0. check if the particle is in the box or not.
    if ( not this->Par_InBox() )    return;
    
    // grid idx   0  |  1  |  2  |  3  |  4  |  5  |
    //            -------------------------------------
    // cell       | L R | L R | L R |     |     |     |
    //            -------------------------------------
    // grid space    |     |     |     |     |     |

    // 1. get particle position left cell index in the grid
    int pos_idx[N_DIMS];
    double dist_to_left[N_DIMS];

    // There is a bug! 1.0/BOX_DX will give wrong result. 
    // printf("_BOX_DX = %.5f", 1.0/BOX_DX);
    const double dx = BOX_DX;
    for ( int d = 0; d < N_DIMS; d++ )
    {
        pos_idx[d]      = this->par_pos[d] / dx;
        dist_to_left[d] = this->par_pos[d] / dx - pos_idx[d]; // in unit of BOX_DX
    } // for ( int d = 0; d < N_DIMS; d++ )

    // 2. calculate particle mass in cell
    #if ( MASS_TO_CELL == NGP )
    this->Par_AddMassToCell_NGP( source, pos_idx, dist_to_left );
    #elif ( MASS_TO_CELL == CIC )
    this->Par_AddMassToCell_CIC( source, pos_idx, dist_to_left );
    #elif ( MASS_TO_CELL == TSC )
    this->Par_AddMassToCell_TSC( source, pos_idx, dist_to_left );
    #endif

} // FUNCTION : Par_AddMassToCell



void particle::Par_AddMassToCell_NGP( matrix &source, const int *pos_idx, const double *dist_to_left )
{
    const double _cell_vol = 1. / pow( BOX_DX, N_DIMS );
    int cell_shift[N_DIMS];     // shift the desposited cell when particle in the right half cell
    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( dist_to_left[d] < 0.5 ) // particle in the left half cell
        {
            cell_shift[d] = 0;
        } else 
        {
            cell_shift[d] = 1;
        }
    } // for ( int d = 0; d < N_DIMS; d++ )

    // deposit the mass to cell
    #if ( N_DIMS == 2 )
    source.add_value( pos_idx[0]+cell_shift[0] + BOX_N*( pos_idx[1]+cell_shift[1] ), 
                          this->par_mass * _cell_vol );
    #elif ( N_DIMS == 3 )
    source.add_value( pos_idx[0]+cell_shift[0] + BOX_N*( pos_idx[1]+cell_shift[1] ) + BOX_N*BOX_N*( pos_idx[2]+cell_shift[2] ), 
                          this->par_mass * _cell_vol );
    #endif

} // FUNCTION : Par_AddMassToCell_NGP



void particle::Par_AddMassToCell_CIC( matrix &source, const int *pos_idx, const double *dist_to_left )
{
    const double _cell_vol = 1. / pow( BOX_DX, N_DIMS );
    int cells = 1;              // Number of cells need to be deposited.
    double frac[N_DIMS][2];     // mass fraction of cell center. [0]:left, [1]:right
    
    for ( int d = 0; d < N_DIMS; d++ )
    {
        frac[d][0] = 1.0 - dist_to_left[d];
        frac[d][1] = dist_to_left[d];

        cells *= 2;
    } // for ( int d = 0; d < N_DIMS; d++ )
    

    // deposit the mass to cell
    for ( int idx = 0; idx < cells; idx++ )
    {
        const int di = idx%2;
        const int dj = (idx%4) / 2;
        const int dk = idx/4;

        #if ( N_DIMS == 2 )
        source.add_value( pos_idx[0]+di + BOX_N*( pos_idx[1]+dj ), 
                            frac[0][di] * frac[1][dj] * this->par_mass * _cell_vol );
        #elif ( N_DIMS == 3 )
        source.add_value( pos_idx[0]+di + BOX_N*( pos_idx[1]+dj ) + BOX_N*BOX_N*( pos_idx[2]+dk ),
                            frac[0][di] * frac[1][dj] * frac[2][dk] * this->par_mass * _cell_vol );
        #endif
    } // for ( int idx = 0; idx < cells; idx++ )

} // FUNCTION : Par_AddMassToCell_CIC



void particle::Par_AddMassToCell_TSC( matrix &source, const int *pos_idx, const double *dist_to_left )
{
    // This is the temporary way to solve the edge case.
    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( pos_idx[d] == 0 || pos_idx[d] == BOX_N-1 )
        {
            this->Par_AddMassToCell_CIC( source, pos_idx, dist_to_left );
            return;
        }
    }

    const double _cell_vol = 1. / pow( BOX_DX, N_DIMS );
    int cell_shift[N_DIMS];     // shift the desposited cell when particle in the right half cell
    int cells = 1;              // Number of cells need to be deposited.
    double frac[N_DIMS][3];     // mass fraction of cell center. [0]:left, [1]:middle, [2]:right
    
    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( dist_to_left[d] < 0.5 ) // particle in the left half cell
        {
            cell_shift[d] = -1;
            frac[d][0] = 0.5  * ( 0.5-dist_to_left[d] ) * ( 0.5-dist_to_left[d] );
            frac[d][1] = 0.75 -       dist_to_left[d]   *       dist_to_left[d];
            frac[d][2] = 0.5  * ( 0.5+dist_to_left[d] ) * ( 0.5+dist_to_left[d] );
        } else 
        {
            cell_shift[d] = 0;
            frac[d][0] = 0.5  * ( 1.5-dist_to_left[d] ) * ( 1.5-dist_to_left[d] );
            frac[d][1] = 0.75 - ( 1.0-dist_to_left[d] ) * ( 1.0-dist_to_left[d] );
            frac[d][2] = 0.5  * ( 0.5-dist_to_left[d] ) * ( 0.5-dist_to_left[d] );
        }

        cells *= 3;

    } // for ( int d = 0; d < N_DIMS; d++ )
    
    // deposit the mass to cell
    for ( int idx = 0; idx < cells; idx++ )
    {
        const int di = idx%3;
        const int dj = (idx%9) / 3;
        const int dk = idx/9;

        #if ( N_DIMS == 2 )
        source.add_value( pos_idx[0]+cell_shift[0]+di + BOX_N*( pos_idx[1]+cell_shift[1]+dj ),
                            frac[0][di] * frac[1][dj] * this->par_mass * _cell_vol );
        #elif ( N_DIMS == 3 )
        source.add_value( pos_idx[0]+cell_shift[0]+di + BOX_N*( pos_idx[1]+cell_shift[1]+dj) + BOX_N*BOX_N*( pos_idx[2]+cell_shift[2]+dk ), 
                            frac[0][di] * frac[1][dj] * frac[2][dk] * this->par_mass * _cell_vol );
        #endif
    } // for ( int idx = 0; idx < cells; idx++ )

} // FUNCTION : Par_AddMassToCell_CIC



void particle::Par_SumAcc( double *acc, double **force )
{
    // 0. check if the particle is in the box or not.
    if ( not this->Par_InBox() )    return;
    
    // grid idx   0  |  1  |  2  |  3  |  4  |  5  |
    //            -------------------------------------
    // cell       | L R | L R | L R |     |     |     |
    //            -------------------------------------
    // grid space    |     |     |     |     |     |

    // 1. get particle position left cell index in the grid
    int pos_idx[N_DIMS];
    double dist_to_left[N_DIMS];

    // There is a bug! 1.0/BOX_DX will give wrong result. 
    // printf("_BOX_DX = %.5f", 1.0/BOX_DX);
    const double dx = BOX_DX;
    for ( int d = 0; d < N_DIMS; d++ )
    {
        pos_idx[d]      = this->par_pos[d] / dx;
        dist_to_left[d] = this->par_pos[d] / dx - pos_idx[d]; // in unit of BOX_DX
    } // for ( int d = 0; d < N_DIMS; d++ )

    // 2. Sum particle acceleration
    #if ( MASS_TO_CELL == NGP )
    this->Par_SumAcc_NGP( acc, force, pos_idx, dist_to_left );
    #elif ( MASS_TO_CELL == CIC )
    this->Par_SumAcc_CIC( acc, force, pos_idx, dist_to_left );
    #elif ( MASS_TO_CELL == TSC )
    this->Par_SumAcc_TSC( acc, force, pos_idx, dist_to_left );
    #endif

} // FUNCTION : particle::Par_SumAcc



void particle::Par_SumAcc_NGP( double *acc, double **force, const int *pos_idx, const double *dist_to_left )
{
    double acc_temp[N_DIMS];
    int cell_shift[N_DIMS];     // shift the desposited cell when particle in the right half cell
    
    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( dist_to_left[d] < 0.5 ) // particle in the left half cell
        {
            cell_shift[d] = 0;
        } else 
        {
            cell_shift[d] = 1;
        }

    } // for ( int d = 0; d < N_DIMS; d++ )

    // sum the acceleration
    for ( int d = 0; d < N_DIMS; d++ )
    {
        const int i = pos_idx[0]+cell_shift[0];
        const int j = pos_idx[1]+cell_shift[1];
        int idx = i + BOX_N * j;
        
        #if ( N_DIMS == 2 )
        acc_temp[d] = force[d][idx];
        #elif ( N_DIMS == 3 )
        const int k = pos_idx[2]+cell_shift[2];
        idx += BOX_N * BOX_N * k;
        acc_temp[d] = force[d][idx];
        #endif
    } // for ( int d = 0; d < N_DIMS; d++ )

    for ( int d = 0; d < N_DIMS; d++ )    acc[d] = acc_temp[d];

} // FUNCTION : particle::Par_SumAcc_NGP



void particle::Par_SumAcc_CIC( double *acc, double **force, const int *pos_idx, const double *dist_to_left )
{
    double acc_temp[N_DIMS];
    int cells = 1;              // Number of cells need to collect.
    double frac[N_DIMS][2];     // mass fraction of cell center. [0]:left, [1]:right
    
    for ( int d = 0; d < N_DIMS; d++ )
    {
        frac[d][0] = 1.0 - dist_to_left[d];
        frac[d][1] = dist_to_left[d];

        cells *= 2;

        acc_temp[d] = 0.0;
    } // for ( int d = 0; d < N_DIMS; d++ )

    // sum the acceleration
    for ( int d = 0; d < N_DIMS; d++ )
    {
        for ( int idx = 0; idx < cells; idx++ )
        {
            const int di = idx%2;
            const int dj = (idx%4) / 2;
            const int dk = idx/4;

            const int i = pos_idx[0] + di;
            const int j = pos_idx[1] + dj;
            int idx_force = i + BOX_N * j;

            #if ( N_DIMS == 2 )
            acc_temp[d] += frac[0][di] * frac[1][dj] * force[d][idx_force];
            #elif ( N_DIMS == 3 )
            const int k = pos_idx[2] + dk;
            idx_force += BOX_N * BOX_N * k;
            acc_temp[d] += frac[0][di] * frac[1][dj] * frac[2][dk] * force[d][idx_force];
            #endif
        } // for ( int idx = 0; idx < cells; idx++ )

    } // for ( int d = 0; d < N_DIMS; d++ )

    for ( int d = 0; d < N_DIMS; d++ )    acc[d] = acc_temp[d];

} // FUNCTION : particle::Par_SumAcc_CIC



void particle::Par_SumAcc_TSC( double *acc, double **force, const int *pos_idx, const double *dist_to_left )
{
    // This is the temporary way to solve the edge case.
    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( pos_idx[d] == 0 || pos_idx[d] == BOX_N-1 )
        {
            this->Par_SumAcc_CIC( acc, force, pos_idx, dist_to_left );
            return;
        }
    }

    double acc_temp[N_DIMS];
    int cell_shift[N_DIMS];     // shift the desposited cell when particle in the right half cell
    int cells = 1;              // Number of cells need to collect.
    double frac[N_DIMS][3];     // mass fraction of cell center. [0]:left, [1]:middle, [2]:right
    
    for ( int d = 0; d < N_DIMS; d++ )
    {
        if ( dist_to_left[d] < 0.5 ) // particle in the left half cell
        {
            cell_shift[d] = -1;
            frac[d][0] = 0.5  * ( 0.5-dist_to_left[d] ) * ( 0.5-dist_to_left[d] );
            frac[d][1] = 0.75 -       dist_to_left[d]   *       dist_to_left[d];
            frac[d][2] = 0.5  * ( 0.5+dist_to_left[d] ) * ( 0.5+dist_to_left[d] );
        } else 
        {
            cell_shift[d] = 0;
            frac[d][0] = 0.5  * ( 1.5-dist_to_left[d] ) * ( 1.5-dist_to_left[d] );
            frac[d][1] = 0.75 - ( 1.0-dist_to_left[d] ) * ( 1.0-dist_to_left[d] );
            frac[d][2] = 0.5  * ( 0.5-dist_to_left[d] ) * ( 0.5-dist_to_left[d] );
        }
        
        cells *= 3;

        acc_temp[d] = 0.0;
    } // for ( int d = 0; d < N_DIMS; d++ )
    
    // sum the acceleration
    for ( int d = 0; d < N_DIMS; d++ )
    {
        for ( int idx = 0; idx < cells; idx++ )
        {
            const int di = idx%3;
            const int dj = (idx%9) / 3;
            const int dk = idx/9;
            
            const int i = pos_idx[0] + di + cell_shift[0];
            const int j = pos_idx[1] + dj + cell_shift[1];
            int idx_force = i + BOX_N * j;

            #if ( N_DIMS == 2 )
            acc_temp[d] += frac[0][di] * frac[1][dj] * force[d][idx_force];
            #elif ( N_DIMS == 3 )
            const int k = pos_idx[2] + dk + cell_shift[2];
            idx_force += BOX_N * BOX_N * k;
        acc_temp[d] += frac[0][di] * frac[1][dj] * frac[2][dk] * force[d][idx_force];
        #endif
        } // for ( int idx = 0; idx < cells; idx++ )

    } // for ( int d = 0; d < N_DIMS; d++ )
    
    for ( int d = 0; d < N_DIMS; d++ )    acc[d] = acc_temp[d];

} // FUNCTION : particle::Par_SumAcc_TSC



void particle::Par_UpdateAll( const double *vel, const double *acc, const double dt )
{
    for ( int d = 0; d < N_DIMS; d++ )
    {
        this->par_pos[d] += vel[d] * dt;
        this->par_vel[d] += acc[d] * dt;
    } // for ( int d = 0; d < N_DIMS; d++ )
    return;

} // FUNCTION : particle::Par_UpdateAll



void particle::Par_SetMass( const double m )
{
    this->par_mass = m;
    return;
} // FUNCTION : particle::Par_SetMass



void particle::Par_SetPos( const double *pos )
{
    for ( int d = 0; d < N_DIMS; d++ )    this->par_pos[d] = pos[d];
    return;
} // FUNCTION : particle::Par_SetPos



void particle::Par_SetVel( const double *vel )
{
    for ( int d = 0; d < N_DIMS; d++ )    this->par_vel[d] = vel[d];
    return;
} // FUNCTION : particle::Par_SetVel



void particle::Par_SetAcc( const double *acc )
{
    for ( int d = 0; d < N_DIMS; d++ )    this->par_acc[d] = acc[d];
    return;
} // FUNCTION : particle::Par_SetAcc



void particle::Par_GetMass( double &m )
{
    m = this->par_mass;
    return;
} // FUNCTION : particle::Par_SetMass



void particle::Par_GetPos( double *pos )
{
    for ( int d = 0; d < N_DIMS; d++ )    pos[d] = this->par_pos[d];
    return;
} // FUNCTION : particle::Par_GetPos



void particle::Par_GetVel( double *vel )
{
    for ( int d = 0; d < N_DIMS; d++ )    vel[d] = this->par_vel[d];
    return;
} // FUNCTION : particle::Par_GetVel



void particle::Par_GetAcc( double *acc )
{
    for ( int d = 0; d < N_DIMS; d++ )    acc[d] = this->par_acc[d];
    return;
} // FUNCTION : particle::Par_GetAcc
