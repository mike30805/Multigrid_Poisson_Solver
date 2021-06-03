// Simulation option 
#define N_DIMS          2
#define BOX_N           2000
#define BOX_L           PI

// Solver
// SOR / V_CYCLE / W_CYCLE / FAS / FMG
#define SOR_OMEGA       1.9
#define POT_SOLVER      V_CYCLE

// Parallel
#define OMP_PARALLEL

// Particle
// NGP / CIC / TSC
#define MASS_TO_CELL    NGP
