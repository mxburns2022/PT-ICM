#ifndef INITIAL
#define INITIAL
#include "mpi_defs.hh"

namespace pticm {
    void read_problem(PTWorker& w, std::string problem_path);
    void initialize_replicas(PTWorker& w);
    void initialize_replica_exchanges(PTWorker& w);
    void initialize_temperatures(PTWorker& w);
    
} // namespace pticm
#endif