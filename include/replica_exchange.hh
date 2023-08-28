#ifndef REPLICA_EXCHANGE
#define REPLICA_EXCHANGE
#include "mpi_defs.hh"
namespace pticm
{
    void try_replica_exchange(PTWorker& w,
        uint8_t chain_index);
    void do_replica_exchange_local(PTWorker& w,
        double beta_old,
        double beta_new,
        uint8_t chain_index);
    void check_for_replica_exchange(PTWorker& w);
} // namespace pticm
#endif