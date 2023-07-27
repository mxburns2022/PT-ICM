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
typedef std::vector<bool> StateVector; // for convenience, use boolean vectors
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
        size_t burn_in; // number of burn-in steps before sampling
    public:
        Problem(std::string path, uint64_t seed = 0, size_t _burn_in = UINT64_MAX);
        double energy (const StateVector& s) const;
        double maxcut (const StateVector& s) const;
        uint64_t vcount () const {return N;};
        double rand_double () {return unif(rng_engine);}
        size_t get_burn_in() const {return burn_in;}
        bool is_sparse () const {return N;}
        double calculate_dE (std::vector<double>& dE, const StateVector& s) const;
        std::mt19937_64& get_rng() {return rng_engine;};
        // Incremental Isakov step
        void incremental_dE(std::vector<double>& dE, const StateVector& s, uint32_t index) const;
};

class Replica {
    private:
        StateVector state; // Current replica state in the problem state bank
        std::vector<double> dE; // deltaE pointer, will indicate location in the problem dE bank
        uint32_t flips; // Number of flips (optional counter)
        double energy; // Current configuration energy
        double beta; // Inverse temperature of this replica
        size_t offset; // 

        // Compare the boltzmann exponent with a log uniform random number
        bool mh_consider (Problem& p, uint32_t index);
        // Accept a proposed MH move
        void mh_accept (Problem& p, uint32_t index);
    public:
        Replica(): energy(0), beta(0), flips(0) {};
        Replica(Problem& p, 
                double _beta, 
                const StateVector& init_state, 
                const std::vector<double>& init_dE, 
                double _energy);
        inline double get_energy() const {return energy;}
        inline const StateVector& get_state() const {return state;}
        inline size_t get_flips() const {return flips;}
        inline void set_state( const StateVector& new_state) { state = new_state;}
        inline void set_dE( std::vector<double> new_dE) { dE = new_dE;}
        inline void set_ene(double new_ene) { energy = new_ene;}
        inline double get_beta() const {return beta;}
        inline void set_beta(const double new_beta) {beta = new_beta;}
        void exchange(Replica& other);
        void log();
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
        size_t exchange_index;
        std::vector<double> beta_values;
        uint32_t sweep_counter;
        size_t num_swaps;
        std::vector<uint32_t> swaps; // store the replica swap order (encoded as flattened indices in an NtxNt array)
        std::vector<size_t> indices; // store the replica indices in order of ascending beta
        std::vector<size_t> correlations; // number of sweeps for each temperature
        // Logs to store data for specific heat, entropy, susceptibility values
        std::vector<std::vector<double>> energy_log;
        std::vector<std::vector<double>> magnetiziation_log;
        std::vector<size_t> flip_logs;
        std::vector<size_t> cons_logs;
        bool log_data;
    public:
        ReplicaChain(Problem& p,
                     uint32_t _Nt,
                     double _beta0,
                     double _beta1,
                     size_t sweeps,
                     bool _log_data = false);
        ReplicaChain(Problem& p,
                     uint32_t _Nt,
                     const std::vector<double>& temperatures,
                     const std::vector<size_t>& corr,
                     size_t sweeps,
                     bool _log_data = false);
        ReplicaChain(): Nt(0) {};
        ReplicaChain& operator=(const ReplicaChain& other){
            replicas = other.replicas;
            Nt = other.Nt;
            num_swaps = other.num_swaps;
            flip_logs = other.flip_logs;
            cons_logs = other.cons_logs;
            beta_values = other.beta_values;
            sweep_counter = other.sweep_counter;
            indices = other.indices;
            swaps = other.swaps;
            energy_log = other.energy_log;
            magnetiziation_log = other.magnetiziation_log;
            correlations = other.correlations;
            log_data = other.log_data;
            exchange_index = other.exchange_index;
            return *this;
        }
        inline const std::vector<std::vector<double>> get_energy_log() const {return energy_log;}
        inline const std::vector<std::vector<double>> get_mag_log() const {return magnetiziation_log;}
        void do_sweeps (Problem& p); 
        size_t get_swaps () const {return num_swaps;}; 
        const std::vector<Replica>& get_replicas () const {return replicas;}; 
        const std::vector<size_t>& get_indices () const {return indices;}; 
        bool re_consider (Problem& p, Replica& r1, Replica& r2);
        void replica_exchanges (Problem& p);
        void print_A_swap() const {
            std::vector<double> output(Nt);
            printf("[");
            for (size_t i = 0; i < Nt-1; i++) {
                printf("%f, ", static_cast<double>(flip_logs[i]) / cons_logs[i]);
            }
            printf("%f]\n", static_cast<double>(flip_logs.back()) / cons_logs.back());
            std::cout << flip_logs << std::endl;
            std::cout << cons_logs << std::endl;
        }
        const StateVector& solution(double& ene);
        
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
                      uint32_t _sweeps,
                      bool record = false);
        ChainEnsemble(Problem& p,
                      uint32_t _Nt,
                      uint32_t _Nc,
                      const std::vector<double>& temps,
                      const std::vector<size_t>& corr,
                      uint32_t _sweeps,
                      bool record = false);
       void solve(Problem& p, double& ene, StateVector& result);
       inline const ReplicaChain& get_chain1  () const {return chain_1;};
       inline const ReplicaChain& get_chain2  () const {return chain_2;};
};
std::ostream& operator<<(std::ostream& out, ReplicaChain& rc);




#endif