#ifndef SOLVER
#define SOLVER
#include "mpi_defs.hh"

namespace pticm {
    void solve(PTWorker& w, std::string& problem_path, double& best_ene, size_t& flips, StateVector& best_state);
} // namespace pticm
#endif