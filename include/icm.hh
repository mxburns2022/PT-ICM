#ifndef ICM
#define ICM
#include "mpi_defs.hh"
namespace pticm {
    void check_icm (PTWorker& w);
    void icm (PTWorker& w, 
            size_t target_replica);
    void bfs_icm(Problem& p, 
                int64_t seed_val, 
                bool* state_1, 
                bool* state_2);
    int64_t xor_icm (PTWorker& w, 
                    bool* state_1, 
                    bool* state_2);
}; // namespace pticm
#endif