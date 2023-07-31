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
#include <boost/dynamic_bitset.hpp>
#include "chain_ensemble.hh"
#ifndef PT_ICM
#define PT_ICM
inline
int sign(bool var) {
    return (static_cast<int>(var)*2 - 1);
}
/**
 * Parallel Tempering + Isoenergetic Cluster Moves: Main approach
 *      * The backbone of the algorithm is Parallel Tempering (AKA Replica Exchange Monte Carlo), which provides ergodicity
 *      * Isoenergetic cluster moves are used to "tunnel" through local minima
 * After a Metropolis-Hastings sweep has been performed, each pair of adjacent replicas are considered for an exchange
 * (Depending on the site percolation threshold of the lattice), ICM is then performed
 * For optimized operation, Isakov-style optimizations will be added
*/

std::ostream& operator<<(std::ostream& out, ReplicaChain& rc);

#endif