#include "macro.h"
#include "validate.h"
#include "classes.h"
#include "particle.h"
#include "initialize.h"
#include "solver.h"
#include "output.h"
#include "evolve.h"
#include "simulation_option.h"
#include "Particle_IC_Constructor.h"

double f( const double x, const double y )
{   
    return sin(x) * sin(y) + BG_POTENTIAL;
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
      
        m.input_answer( idx, f(x, y) );
    } // for ( int idx = 0; idx < dim_2; idx++ )

} // FUNCTION : solved


void test_particle_ic_constructor(){
    double Newton_G =1.0;
    double Rho0          =  1.0;     // peak density [unit: M mass of sun/kpc^3]
    double R0            =  1.0 ;                 // scale radius [unit: kpc]
    double MaxR          =  2.0;                    // maximum radius for particles [unit: kpc]
    int MassProfNBin   = 1000    ;             // number of radial bins in the mass profile table [1000]
    double Alpha         =1     ;               // alpha parameter for Eiasto model [1]
    int r_col         = 0         ;           // number of the column of the radius of density profile, when model is "UNKNOWN"  [0]
    int rho_col        =1          ;          // number of the column of the density of density profile, when model is "UNKNOWN"  [1]
    double truncation    = 0       ;             // whether to turn on a smoothy truncation function of density near MaxR [0]

    Particle_IC_Constructor constructor_Models;
    constructor_Models.init("UNKNOWN",Alpha,Newton_G,Rho0,R0,MassProfNBin,MaxR,truncation,0.7,r_col,rho_col,"halo_profile.txt");

    //Example Code for Constructing Particle ICs (May be modified in your present code)
    const int N = 1000;
    double par_r [N];
    double par_vel[N];
    int count = 0;
    for (int p=0;p<N;p++){
        par_r[p] = constructor_Models.set_radius();
        //cout<<par_r[p]<<endl;
        par_vel[p] = constructor_Models.set_vel(par_r[p]);//par_r[p]/R0
        //cout<<par_vel[p]<<endl;
        //cout<<p<<endl;
        cout<<par_r[p]<<endl;
        if (par_r[p]>MaxR){
            
            count++;
        }
    }
    cout<<count<<endl;
}


int main()
{
    //test_particle_ic_constructor();
    
    if ( not Validate() ) return 0;
    bool init_status;
    
    //particle *pars = NULL;
    particle *pars = new particle[N_PARS];

    matrix pot( BOX_N, BOX_DX );
    pot.init_potential();
    
    matrix dens( BOX_N, BOX_DX );
    init_status = Init_matrix( dens, pars );
    if ( not init_status ) return 0;
    //dens.init_density();
    
    matrix ans( BOX_N, BOX_DX );
    solved(ans);
    
    
    //auto start = chrono::steady_clock::now();

    // solve potential
    matrix solution = Solver_Potential( pot, dens );
    
    //auto elapsed = chrono::steady_clock::now() - start;
    //auto sec_double = chrono::duration<double>(elapsed);     // double
    //cout << "Potential solve time: " << sec_double.count() << "(s)" << endl;
    
    // Output
    int output_counter = 0;
    char density_filename[50], potential_filename[50], particle_filename[50];
    
    sprintf( density_filename, "Density_%d%d.txt", (output_counter%100)/10, output_counter%10 );
    sprintf( potential_filename, "Potential_%d%d.txt", (output_counter%100)/10, output_counter%10 );
    sprintf( particle_filename, "Particle_%d%d.txt", (output_counter%100)/10, output_counter%10 );
    
    Output_matrix( dens, density_filename );
    Output_matrix( solution, potential_filename );
    Output_particles( pars, particle_filename );
    
    solution.Error( ans ); // print the error
    
    /*
    //==============================================================
    int output_counter = 0;
    double time_now = 0.0;
    double dt;
    bool init_status = false, out_stat = false, end_stat = false;
    const double dx = BOX_DX;
    matrix dens( BOX_N, dx );
    particle *pars = new particle[N_PARS];
    char density_filename[50], potential_filename[50], particle_filename[50];
    
    // 1. Validate the simulation options.
    Output_parameter();
    if ( not Validate() ) return 0;
    
    // 2. Initialize the matrix and the particle.
    init_status = Init_matrix( dens, pars );
    if ( not init_status ) return 0;    // End the simulation, if the initialization is failed.
    
    // 2-a. Output the initial condition.
    sprintf( density_filename,  "Density_%d%d.txt",  (output_counter%100)/10, output_counter%10 );
    sprintf( particle_filename, "Particle_%d%d.txt", (output_counter%100)/10, output_counter%10 );
    
    Output_matrix( dens, density_filename );
    Output_particles( pars, particle_filename );
   
    //matrix pot( BOX_N, dx );
    //pot.init_potential();
    //matrix solved_pot = Solver_Potential( pot, dens );
    //sprintf( potential_filename,  "Potential_%d%d.txt",  (output_counter%100)/10, output_counter%10 );
    //Output_matrix( solved_pot, potential_filename );
    output_counter += 1;
    
    // 3. Time evolution.
    while ( not end_stat )
    {
        dt = Evolve_GetDt( pars, time_now, output_counter, out_stat, end_stat );
        dens.reset();
        for ( int p = 0; p < N_PARS; p++ )   pars[p].Par_AddMassToCell( dens );
        Evolve_UpdateParticle( dens, pars, dt );
        
        time_now += dt;
        printf("T_now = %.5f, dt = %.5f\n", time_now, dt);
        
        // Output Data
        if ( out_stat )
        {
            //pot.init_potential();
            //matrix solved_pot = Solver_Potential( pot, dens );
            //sprintf( potential_filename,  "Potential_%d%d.txt",  (output_counter%100)/10, output_counter%10 );
            //Output_matrix( solved_pot, potential_filename );
            
            sprintf( density_filename,  "Density_%d%d.txt",  (output_counter%100)/10, output_counter%10 );
            sprintf( particle_filename, "Particle_%d%d.txt", (output_counter%100)/10, output_counter%10 );
            
            Output_matrix( dens, density_filename );
            Output_particles( pars, particle_filename );
            
            output_counter += 1;
            out_stat = false;
            printf("Data Output at T_now = %.5f\n", time_now);
        } // if ( out_stat )
        
    } // while ( time_now <= END_TIME )

    //dump_data();
    
    //==============================================================

    delete[] pars;
    */
} // FUNCTION : main 
