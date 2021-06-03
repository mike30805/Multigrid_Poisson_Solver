#include "particle.h"


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



void particle::Par_AddMassToCell( double **source )
{
    // grid idx   0  |  1  |  2  |  3  |  4  |  5  |
    //            -------------------------------------
    // cell       | L R | L R | L R |     |     |     |
    //            -------------------------------------
    // grid space    |     |     |     |     |     |

    const double _cell_vol = 1. / pow( BOX_DX, N_DIMS );
    
    // 1. get particle position left cell index in the grid
    int pos_idx[N_DIMS];
    double dist_to_left[N_DIMS];

    for ( int d = 0; d < N_DIMS; d++ )
    {
        pos_idx[d]      = this->par_pos[d] / BOX_DX;
        dist_to_left[d] = par_pos[d]/BOX_DX - pos_idx[d]; // in unit of BOX_DX
    } // for ( int d = 0; d < N_DIMS; d++ )

    // 2. calculate particle mass in cell
    #if ( MASS_TO_CELL == NGP )
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

    source[ pos_idx[0]+cell_shift[0] ][ pos_idx[1]+cell_shift[1] ] += this->par_mass * _cell_vol;
    
    #elif ( MASS_TO_CELL == CIC )
    double L_frac[N_DIMS];      // mass fraction of left cell
    for ( int d = 0; d < N_DIMS; d++ )    L_frac[d] = 1.0 - dist_to_left[d];

    // deposit the mass to cell
    source[ pos_idx[0]   ][ pos_idx[1]   ] +=     L_frac[0]  *     L_frac[1]  * this->par_mass * _cell_vol;
    source[ pos_idx[0]   ][ pos_idx[1]+1 ] +=     L_frac[0]  * (1.-L_frac[1]) * this->par_mass * _cell_vol;
    source[ pos_idx[0]+1 ][ pos_idx[1]   ] += (1.-L_frac[0]) *     L_frac[1]  * this->par_mass * _cell_vol;
    source[ pos_idx[0]+1 ][ pos_idx[1]+1 ] += (1.-L_frac[0]) * (1.-L_frac[1]) * this->par_mass * _cell_vol;

    #elif ( MASS_TO_CELL == TSC )
    // empty for now
    #endif


} // FUNCTION : Par_AddMassToCell



void particle::Par_UpdateAll( const double *vel, const double *acc, const double dt )
{
    for ( int d = 0; d < N_DIMS; d++ )
    {
        this->par_pos[d] += vel[d] * dt;
        this->par_vel[d] += acc[d] * dt;
    } // for ( int d = 0; d < N_DIMS; d++ )
    return;

} // FUNCTION : particle::Par_UpdateAll



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
