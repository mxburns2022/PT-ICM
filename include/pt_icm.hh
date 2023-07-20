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
template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> vec){
    out << "[";
    if (vec.size() > 0) {
        for (size_t i = 0; i < vec.size() - 1; i++) {
            out << vec.at(i) << ", ";
        }
        out << vec.back();

    }
    out << "]";
    return out;
};
struct Edge {
    size_t dest;
    double weight;
};
typedef std::vector<int8_t> StateVector; // for convenience, use boolean vectors
std::ostream& operator<<(std::ostream& out, StateVector vec);
typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<std::list<Edge>> SparseMatrix;

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
    public:
        Problem(std::string path, uint64_t seed = 0);
        double energy (const StateVector* s) const;
        double energy (StateVector s) const;
        uint64_t vcount () const {return N;};
        double rand_double () {return unif(rng_engine);}
        bool is_sparse () const {return N;}
        double calculate_dE (std::vector<double>* dE, StateVector* s) const;

        // Incremental Isakov step
        void incremental_dE(std::vector<double>* dE, StateVector* s, uint32_t index) const;
};

class Replica {
    private:
        StateVector* state; // Current replica state in the problem state bank
        std::vector<double>* dE; // deltaE pointer, will indicate location in the problem dE bank
        uint32_t flips; // Number of flips (optional counter)
        double energy; // Current configuration energy
        double beta; // Inverse temperature of this replica

        // Compare the boltzmann exponent with a log uniform random number
        bool mh_consider (Problem& p, uint32_t index);
        // Accept a proposed MH move
        void mh_accept (Problem& p, uint32_t index);
    public:
        Replica(): energy(0), beta(0), flips(0) {};
        Replica(Problem& p, double _beta, StateVector* _state, std::vector<double>* _dE);
        inline double get_energy() const {return energy;}
        inline StateVector* get_state() const {return state;}
        inline size_t get_flips() const {return flips;}
        inline void set_state( StateVector* new_state) { state = new_state;}
        inline void set_dE( std::vector<double>* new_dE) { dE = new_dE;}
        inline void set_ene(double new_ene) { energy = new_ene;}
        inline double get_beta() const {return beta;}
        void exchange(Replica& other);
        Replica& operator=(const Replica& other){
            state = other.state;
            dE = other.dE;
            energy = other.energy;
            beta = other.beta;
            flips = other.flips;
            return *this;
        }
        void sweep (Problem& p);
};
class ReplicaChain {
    private:
        std::vector<Replica> replicas;
        uint32_t Nt;
        double beta0;
        double beta1;
        std::vector<std::vector<double>> dE_bank;
        std::vector<StateVector> state_bank;
    public:
        ReplicaChain(Problem& p,
                     uint32_t _Nt,
                     double _beta0,
                     double _beta1);
        ReplicaChain(): Nt(0), beta0(0), beta1(0) {};
        ReplicaChain& operator=(const ReplicaChain& other){
            replicas = other.replicas;
            dE_bank = other.dE_bank;
            state_bank = other.state_bank;
            size_t i = 0;
            for (Replica& r: replicas) {
                r.set_state(&state_bank[i]);
                r.set_dE(&dE_bank[i++]);
            }
            Nt = other.Nt;
            beta0 = other.beta0;
            beta1 = other.beta1;
            return *this;
        }
        void do_sweeps (Problem& p); 
        const std::vector<Replica>& get_replicas () {return replicas;}; 
        bool re_consider (Problem& p, Replica& r1, Replica& r2);
        void replica_exchanges (Problem& p);
        StateVector* solution(double& ene);
        
};
class ChainEnsemble {
    private:
        ReplicaChain chain_1;
        ReplicaChain chain_2;
        uint32_t Nc;
        uint32_t sweeps;
    public:
        ChainEnsemble(Problem& p,
                      uint32_t _Nt,
                      uint32_t _Nc,
                      double _beta0,
                      double _beta1,
                      uint32_t _sweeps);
       void solve(Problem& p, double& ene, StateVector& result);
};
std::ostream& operator<<(std::ostream& out, ReplicaChain& rc);




#endif