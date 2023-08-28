#include <vector>
#include <cstdint>
#include <random>
#include <string>
#include <list>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>
#include <stack>
#include <boost/dynamic_bitset.hpp>
#include "mpi_defs.hh"
#ifndef PT_ICM
#define PT_ICM
/**
 * Parallel Tempering + Isoenergetic Cluster Moves: Main approach
 *      * The backbone of the algorithm is Parallel Tempering (AKA Replica Exchange Monte Carlo), which provides ergodicity
 *      * Isoenergetic cluster moves are used to "tunnel" through local minima
 * After a Metropolis-Hastings sweep has been performed, each pair of adjacent replicas are considered for an exchange
 * (Depending on the site percolation threshold of the lattice), ICM is then performed
 * For optimized operation, Isakov-style optimizations will be added
*/
namespace pticm
{
    void bfs_icm(Problem& p, int64_t seed_val, bool* state_1, bool* state_2);
    int64_t xor_icm(PTWorker& w, bool* state_1, bool* state_2);
    Replica& get_replica(PTWorker& w, double beta, uint8_t chain);
    void icm(PTWorker& w, size_t target_replica);
    void check_icm(PTWorker& w);

} // namespace pticm

#endif