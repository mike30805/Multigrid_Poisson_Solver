// Simulation option 
#define N_DIMS          2
#define BOX_N           200
#define BOX_L           10.0

// Solver
// SOR / V_CYCLE / W_CYCLE / FAS / FMG
#define POT_SOLVER      V_CYCLE

// Parallel
#define OMP_PARALLEL

// Particle
// NGP / CIC / TSC
#define MASS_TO_CELL    NGP
