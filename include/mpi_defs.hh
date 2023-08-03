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

struct Parameters {
    // Command line parameters
    uint64_t seed;
    uint64_t Nt;
    uint64_t pts;
    uint64_t burn_in;
    double T0;
    double T1;
    bool inverse;
    bool log_data;
    uint32_t verbosity;
};

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
    double* energy_values_global;  // energy values from all replicas (updated each iteration)
    uint32_t* chain_1_exchange_order; // Sequence from 1...M-1 for pair exchanges for each chain (shuffled)
    uint32_t* chain_2_exchange_order;
    uint8_t* chain_index_global; // the chain index for all replicas
    uint32_t exchange_index; // next pair exchange to attempt
    uint32_t* beta_indices; // Map from replica index to beta index. beta_values[beta_indices[i]] gives the current inverse temp of replica i
    int* replica_workers; // map from replica index to current MPI worker rank



    uint8_t* chain_index_local; // store which chain each replica belongs to (avoids confusion)
    double* energy_values_local; // energy values from replicas in this worker
    double* beta_values_local;  // beta values for replicas in this worker
    Problem p; // problem instance
    std::vector<Replica> replicas;
    Parameters params;
};
int parse_args (int argc, 
                 char** argv,
                 std::string& problem_path,
                 Parameters& params);
void free_worker (PTWorker& w);
void create_datatype(MPI_Datatype& param_type);
void beta_linspace(double* temps, double T0, double T1, uint64_t num_temp, bool inverse);
void initialize_temperatures(PTWorker& worker_data);
void read_problem(PTWorker& manager_data, std::string problem_path);
void reconstruct_problem(PTWorker& worker_data);
void initialize_replicas(PTWorker& worker_data);
void initialize_replica_exchanges(PTWorker& manager_data);
void gather_energy_values(PTWorker& manager_data);
void sent_energy_values(PTWorker& work_data);
bool try_replica_exchange(PTWorker& manager_data, uint8_t chain_index, uint32_t exchange_location);
void do_replica_exchange_remote(PTWorker& manager_data, uint32_t worker_rank, double beta_old, double beta_new, uint8_t chain_index);
void do_replica_exchange_local(PTWorker& worker_data, double beta_old, double beta_new, uint8_t chain_index);
void check_for_replica_exchange(PTWorker& worker_data);
void solve(PTWorker& w, std::string& problem_path);

#endif