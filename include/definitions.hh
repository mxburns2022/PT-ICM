#ifndef DEFS
#define DEFS
#include <list>
#include <vector>
#include <iostream>
#include "mpi.h"
struct Parameters {
    uint64_t burn_in;
    uint64_t sweeps;
    uint64_t Nt;
    uint64_t beta_min;
    uint64_t beta_max;
    
};
template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> vec){
    out << "[";
    if (vec.size() > 0) {
        for (uint64_t i = 0; i < vec.size() - 1; i++) {
            out << vec.at(i) << ", ";
        }
        out << vec.back();

    }
    out << "]";
    return out;
};
struct Edge {
    uint64_t dest;
    double weight;
};
typedef std::vector<bool> StateVector; // for convenience, use boolean vectors
std::ostream& operator<<(std::ostream& out, StateVector vec);
typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<std::list<Edge>> SparseMatrix;

#endif