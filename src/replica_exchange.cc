#include "pt_icm.hh"
#include "replica_exchange.hh"
namespace pticm {
    // attempt replica exchange using metropolis-hastings criterion
    void try_replica_exchange(PTWorker& w,
                            uint8_t chain_index) {
        assert(w.wid == 0);
        uint64_t Nt = w.params.Nt;
        
        uint64_t chain_exchange_index;
        if (chain_index == 0) {
            chain_exchange_index = w.chain_1_exchange_order[w.exchange_index];
        } else {
            chain_exchange_index = w.chain_1_exchange_order[w.exchange_index] + Nt;
        }
        w.swap_attempts[chain_exchange_index] ++;
        uint32_t replica_1_ind = w.beta_indices[chain_exchange_index];
        uint32_t replica_2_ind = w.beta_indices[chain_exchange_index+1];
        uint32_t worker_1 = w.replica_workers[replica_1_ind];
        uint32_t worker_2 = w.replica_workers[replica_2_ind];

        double beta_low = w.beta_values[chain_exchange_index];
        double beta_high = w.beta_values[chain_exchange_index+1];
        double dBeta = beta_low - beta_high;
        double dE = w.energy_values_global[replica_1_ind] -\
                    w.energy_values_global[replica_2_ind]; 
        double log_boltzmann = dBeta * dE;
        double rndcomp = std::log(w.p.rand_double());
        if (w.params.verbosity >= 3) {
            std::cout << "\n\n(%d %d)\ndB: " << dBeta << std::endl;
            std::cout << "dE: " << dE << std::endl;
            std::cout << "B1: " << beta_low << " B2:" << beta_high  << std::endl;
            std::cout << "I1: " << replica_1_ind << " I2:" << replica_2_ind  << std::endl;
            std::cout << "E1: " << w.energy_values_global[replica_1_ind] << " E2:" << w.energy_values_global[replica_2_ind]  << std::endl;
            std::cout << "log(dB*dE): " <<  log_boltzmann << std::endl;
            std::cout<< "RND: " <<  rndcomp << std::endl;
        }
        // getchar();
        // attempt replica exchange using metropolis-hastings criterion
        if ((dE * dBeta) > 0 || rndcomp < log_boltzmann) {
            w.swap_successes[chain_exchange_index] ++;
            if (w.params.verbosity >= 2) {
                printf("\nSWAPPING temperatures %f <-> %f, dE = %f, in replicas %d and %d\n", beta_low, beta_high, dE, worker_1, worker_2); 
                fflush(stdout);
            }
            w.swaps ++;
            bool swap = true;
            //NOTE: do_replica_exchange_local will automatically perform both swaps if the same worker has both replicas
            if (worker_1 == 0) {
                // perform the replica exchange in replica 1 locally
                do_replica_exchange_local(w, beta_low, beta_high, chain_index);
            } else {
                // send the old and new beta values to worker 1
                MPI_Send(&swap, 1, MPI_CXX_BOOL, worker_1, 0, MPI_COMM_WORLD);
                MPI_Send(&beta_high, 1, MPI_DOUBLE, worker_1, 0, MPI_COMM_WORLD);
                MPI_Send(&beta_low, 1, MPI_DOUBLE, worker_1, 0, MPI_COMM_WORLD);
                MPI_Send(&chain_index, 1, MPI_UNSIGNED_CHAR, worker_1, 0, MPI_COMM_WORLD);
            }
            // send info to worker with the first replica
            // send info to worker with the second replica
            if (worker_1 != worker_2) {
                if (worker_2 == 0) {
                    // perform the replica exchange in replica 2 locally
                    do_replica_exchange_local(w, beta_high, beta_low, chain_index);
                } else {
                    // send the old and new beta values
                    MPI_Send(&swap, 1, MPI_CXX_BOOL, worker_2, 0, MPI_COMM_WORLD);
                    MPI_Send(&beta_low, 1, MPI_DOUBLE, worker_2, 0, MPI_COMM_WORLD);
                    MPI_Send(&beta_high, 1, MPI_DOUBLE, worker_2, 0, MPI_COMM_WORLD);
                    MPI_Send(&chain_index, 1, MPI_UNSIGNED_CHAR, worker_2, 0, MPI_COMM_WORLD);
                }
            }
            // swap the beta indices
            std::swap(w.beta_indices[chain_exchange_index], w.beta_indices[chain_exchange_index + 1]);
        } else {
            // tell the workers that no swap will occur
            bool swap = false;
            if (worker_1 != 0)
                MPI_Send(&swap, 1, MPI_CXX_BOOL, worker_1, 0, MPI_COMM_WORLD);

            if (worker_1 != worker_2) {
                if (worker_2 != 0)
                    MPI_Send(&swap, 1, MPI_CXX_BOOL, worker_2, 0, MPI_COMM_WORLD);
            }
        }
        w.exchange_index++;
        if (w.exchange_index == w.params.Nt - 1) {
            w.exchange_index = 0;
            std::shuffle(w.chain_1_exchange_order, w.chain_1_exchange_order + w.params.Nt - 1, w.p.get_rng());
            std::shuffle(w.chain_2_exchange_order, w.chain_2_exchange_order + w.params.Nt - 1, w.p.get_rng());
        }
        for (uint64_t i = 1; i < w.numw; i++) {
            bool dummy = false;
            if (i == worker_1 || i == worker_2) {
                continue;
            }
            MPI_Send(&dummy, 1, MPI_CXX_BOOL, i, 0, MPI_COMM_WORLD);
        }
        
    };

    // swap the replica beta values locally if the exchange was accepted
    void do_replica_exchange_local(PTWorker& w, 
                                double beta_old, 
                                double beta_new,
                                uint8_t chain_index) {
    
        bool check = false;
        for (uint64_t j = 0; j < w.loc_r; j++) {
            // will perform both swaps if the worker contains both beta_old and beta_new replicas in `chain_index`
            if (std::abs(w.beta_values_local[j] - beta_new) < 1e-10 && w.chain_index_local[j] == chain_index) {
                Replica& r = w.replicas[j];
                w.beta_values_local[j] = beta_old;
                r.set_beta(beta_old);
                check = true;
            }
            else if (std::abs(w.beta_values_local[j] - beta_old) < 1e-10 && w.chain_index_local[j] == chain_index) {
                Replica& r = w.replicas[j];
                w.beta_values_local[j] = beta_new;
                r.set_beta(beta_new);
                check = true;
            }
        }   
        assert(check);
    };

    // check for messages from manager regarding replica exchanges
    void check_for_replica_exchange(PTWorker& w) {
        bool exchange = false;
        MPI_Recv(&exchange, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, &w.stat);
        if (exchange) {
            double beta_old, beta_new;
            uint8_t chain_index;
            MPI_Recv(&beta_old, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &w.stat);
            MPI_Recv(&beta_new, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &w.stat);
            MPI_Recv(&chain_index, 1, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &w.stat);
            do_replica_exchange_local(w, beta_old, beta_new, chain_index);   
        }
    };
} // namespace pticm 