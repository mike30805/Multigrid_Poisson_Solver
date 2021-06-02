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




void particle::update_all( const double *vel, const double *acc, const double dt )
{
    for ( int d = 0; d < N_DIMS; d++ )
    {
        this->par_pos[d] += vel[d] * dt;
        this->par_vel[d] += acc[d] * dt;
    } // for ( int d = 0; d < N_DIMS; d++ )
    return;

} // FUNCTION : particle::update_all



void particle::set_pos( const double *pos )
{
    for ( int d = 0; d < N_DIMS; d++ )    this->par_pos[d] = pos[d];
    return;
} // FUNCTION : particle::set_pos



void particle::set_vel( const double *vel )
{
    for ( int d = 0; d < N_DIMS; d++ )    this->par_vel[d] = vel[d];
    return;
} // FUNCTION : particle::set_vel



void particle::set_acc( const double *acc )
{
    for ( int d = 0; d < N_DIMS; d++ )    this->par_acc[d] = acc[d];
    return;
} // FUNCTION : particle::set_acc



void particle::get_pos( double *pos )
{
    for ( int d = 0; d < N_DIMS; d++ )    pos[d] = this->par_pos[d];
    return;
} // FUNCTION : particle::get_pos



void particle::get_vel( double *vel )
{
    for ( int d = 0; d < N_DIMS; d++ )    vel[d] = this->par_vel[d];
    return;
} // FUNCTION : particle::get_vel



void particle::get_acc( double *acc )
{
    for ( int d = 0; d < N_DIMS; d++ )    acc[d] = this->par_acc[d];
    return;
} // FUNCTION : particle::get_acc
