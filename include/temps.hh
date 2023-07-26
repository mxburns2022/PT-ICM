#include <vector>
#ifndef TEMPS
#define TEMPS
// Temperatures are taken from Malakis and Papakonstaninou 2013 "Comparative Study of Selected Parallel Tempering Methods"

// Temperatures for Constant Acceptance Exchange (CAE) for a 50x50 2D Ising model (ferromagnetic) with M=19
//      Constant exchange probability should be r=0.5
const std::vector<double> cae_ising_2d = {
    1.92,
    1.9669,
    2.0121,
    2.0557,
    2.0975,
    2.1377,
    2.1757,
    2.2115,
    2.2446,
    2.2751,
    2.3050,
    2.3374,
    2.3746,
    2.4167,
    2.4631,
    2.5134,
    2.568,
    2.6268,
    2.6903};
// The correlation time for each temperature in the Canonical ensemble (in MCMC sweeps)
const std::vector<size_t> canonical_correlation_sweeps = {
    3,
    3,
    3,
    4,
    4,
    5,
    7,
    10,
    20,
    33,
    34,
    20,
    10,
    7,
    5,
    5,
    4,
    4,
    3
};
const std::vector<size_t> single_correlation_sweeps = {
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1
};
// Temperatures for a Constant Entropy Increase (CEI) with M=19
const std::vector<double> cei_ising_2d = {
    1.92,
    1.9825,
    2.0375,
    2.0875,
    2.1325,
    2.1725,
    2.2075,
    2.2375,
    2.2625,
    2.2850,
    2.3075,
    2.3325,
    2.3625,
    2.4000,
    2.4425,
    2.4925,
    2.5525,
    2.6200,
    2.6975};
// CAE, M=11 Temperatures for the L=6 cubic EA lattice with bimodal couplings



#endif