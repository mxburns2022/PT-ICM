#include <vector>
#ifndef TEMPS
#define TEMPS
// Temperatures are taken from Malakis and Papakonstaninou 2013 "Comparative Study of Selected Parallel Tempering Methods"

// Temperatures for Constant Acceptance Exchange (CAE) for a 50x50 2D Ising model (ferromagnetic) with M=19
//      Constant exchange probability should be r=0.5
const std::vector<double> cae_ising_2d = {
    1/1.92,
    1/1.9669,
    1/2.0121,
    1/2.0557,
    1/2.0975,
    1/2.1377,
    1/2.1757,
    1/2.2115,
    1/2.2446,
    1/2.2751,
    1/2.3050,
    1/2.3374,
    1/2.3746,
    1/2.4167,
    1/2.4631,
    1/2.5134,
    1/2.568,
    1/2.6268,
    1/2.6903};
const std::vector<double> cae_test = {
    0.37,
    0.383672,
    0.387797,
    0.402753,
    0.41229,
    0.418768,
    0.425246,
    0.431724,
    0.438202,
    0.444679,
    0.451157,
    0.457635,
    0.464113,
    0.470591,
    0.477068,
    0.483546,
    0.490024,
    0.496502,
    0.502979
};
// The correlation time for each temperature in the Canonical ensemble (in MCMC sweeps)
const std::vector<uint64_t> canonical_correlation_sweeps = {
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
const std::vector<uint64_t> single_correlation_sweeps = {
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
    1/1.92,
    1/1.9825,
    1/2.0375,
    1/2.0875,
    1/2.1325,
    1/2.1725,
    1/2.2075,
    1/2.2375,
    1/2.2625,
    1/2.2850,
    1/2.3075,
    1/2.3325,
    1/2.3625,
    1/2.4000,
    1/2.4425,
    1/2.4925,
    1/2.5525,
    1/2.6200,
    1/2.6975};
// CAE, M=11 Temperatures for the L=6 cubic EA lattice with bimodal couplings



#endif