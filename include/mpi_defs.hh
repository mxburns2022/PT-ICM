#include <mpi.h>
#include <random>
#include <memory>
#include <cstring>
#include <cassert>
#include <algorithm>
#include "problem.hh"
#include "replica.hh"
#ifndef MPIDEFS
#define MPIDEFS
namespace pticm {
    struct Parameters {
        // Command line parameters
        uint64_t seed;
        uint64_t Nt;
        uint64_t pts;
        uint64_t burn_in;
        uint64_t sweeps;
        double T0;
        double T1;
        double cutoff;
        bool inverse;
        bool log_data;
        bool icm;
        uint32_t verbosity;
    };
    std::vector<double> get_cae(Problem& p, uint64_t Nt, uint64_t trials);
    std::ostream& operator<<(std::ostream& out, std::vector<uint8_t> vec);
    std::ostream& operator<<(std::ostream& out, std::list<Edge> vec);
    struct PTWorker {
        int wid; // MPI worker rank
        int numw; // number of MPI workers
        uint64_t worker_r;
        uint64_t last_r;
        uint64_t loc_r;
        MPI_Status stat;
        // Fields used only by maanger (id = 0)
        double* beta_values; // vector storing all beta values (sorted high to low)
        uint64_t swaps;
        uint32_t Nc;
        double characteristic_scale; // cutoff temperature for ICM
        double* energy_values_global;  // energy values from all replicas (updated each iteration)
        uint32_t* chain_1_exchange_order; // Sequence from 1...M-1 for pair exchanges for each chain (shuffled)
        uint32_t* chain_2_exchange_order;
        uint8_t* chain_index_global; // the chain index for all replicas
        uint32_t exchange_index; // next pair exchange to attempt
        uint32_t* beta_indices; // Map from replica index to beta index. beta_values[beta_indices[i]] gives the current inverse temp of replica i
        int* replica_workers; // map from replica index to current MPI worker rank
        // Record the number of swap attempts vs swap successes
        std::vector<uint32_t> swap_attempts; 
        std::vector<uint32_t> swap_successes;




        uint8_t* chain_index_local; // store which chain each replica belongs to (avoids confusion)
        double* energy_values_local; // energy values from replicas in this worker
        double* beta_values_local;  // beta values for replicas in this worker
        Problem p; // problem instance
        std::vector<Replica> replicas;
        Parameters params;
    };

} // namespace pticm

#endif