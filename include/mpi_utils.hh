#ifndef MPI_UTILS
#define MPI_UTILS
#include "mpi_defs.hh"

namespace pticm {
    void free_worker (PTWorker& w);
    uint64_t get_flips(PTWorker& w);
    void gather_energy_values(PTWorker& w);
    void check_correct_beta_values(PTWorker& w);
    void create_datatype(MPI_Datatype& param_type);
    Replica& get_replica(PTWorker& w, 
                        double beta, 
                        uint8_t chain);
    void gather_best_energy_values(PTWorker& w, 
                                double& best_ene_global, 
                                StateVector& best_state);
} //namespace pticm
#endif