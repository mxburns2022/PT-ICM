#ifndef PROBLEM
#define PROBLEM
#include <random>
#include "definitions.hh"
class Problem {
    private:
        bool sparse; // if the problem is sparse (in which case the matrix is an edgelist)
        Matrix formulation; // use whichever one is determined at runtime
        SparseMatrix formulation_sparse;
        uint64_t N; // number of nodes
        uint64_t M; // number of edges
        double p_c; // Percolation threshold
        double char_ene; // Characteristic energy scale
        double best_ene; // best energy found so far
        StateVector best_solution; // best solution found so far (lowest energy)
        std::mt19937_64 rng_engine; // Mersenne Twister random engine for RNG
        std::uniform_real_distribution<double> unif;
        size_t burn_in; // number of burn-in steps before sampling
    public:
        Problem(std::string path, uint64_t seed = 0, size_t _burn_in = UINT64_MAX);
        double energy (const StateVector& s) const;
        double maxcut (const StateVector& s) const;
        uint64_t vcount () const {return N;};
        double rand_double () {return unif(rng_engine);}
        size_t get_burn_in() const {return burn_in;}
        const std::list<Edge>& get_sparse_neighbors(size_t ind) const {return formulation_sparse[ind];}
        const std::vector<double>& get_dense_neighbors(size_t ind) const {return formulation[ind];}
        bool is_sparse () const {return N;}
        double calculate_dE (std::vector<double>& dE, const StateVector& s) const;
        std::mt19937_64& get_rng() {return rng_engine;};
        // Incremental Isakov step
        void incremental_dE(std::vector<double>& dE, const StateVector& s, uint32_t index) const;
};
#endif