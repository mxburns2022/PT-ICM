#include "solver.hh"
#include "problem.hh"
#include "icm.hh"
#include "replica_exchange.hh"
#include "mpi_utils.hh"
#include "utils.hh"
#include "initialization.hh"
#include "definitions.hh"
#include <boost/timer/progress_display.hpp>
namespace bt = boost::timer; 
namespace pticm {
    void solve(PTWorker& w, std::string& problem_path, double& best_ene, size_t& flips, StateVector& best_state) {
        MPI_Datatype parameter_type;
        create_datatype(parameter_type);
        // send the relavent parameters to each worker
        MPI_Bcast((void*)&w.params, 1, parameter_type, 0, MPI_COMM_WORLD);
        uint64_t Nt = w.params.Nt;
        // Number of replicas to be handled by each worker
        w.worker_r = Nt * 2 / w.numw;
        uint64_t regular = w.worker_r * w.numw;
        w.last_r = w.worker_r;
        if (Nt != regular) {
            w.last_r += (Nt * 2) - regular; 
        }
        w.loc_r = w.wid != w.numw - 1 ? w.worker_r : w.last_r;
        MPI_Barrier(MPI_COMM_WORLD);
        // instantiate the problem
        read_problem(w, problem_path);

        // initialize the temperatures (either preprocessed CAE or linearly spaced)
        write_update(w, "Initializing Temperatures", 1);
        initialize_temperatures(w);

        // initialize replicas and send info to workers
        write_update(w, "Initializing Replicas", 1);
        initialize_replicas(w);

        // Initialize the replica exchange data structures
        write_update(w, "Initializing Exchanges", 1);
        initialize_replica_exchanges(w);
        // ensure that each worker is present before starting loop
        MPI_Barrier(MPI_COMM_WORLD);
        write_update(w, "Starting MCMC loop", 1);
        double progress;
        size_t display_index = 0;
        std::vector<char> display;
        double itracker;
        double trackerinc;

        if (w.wid == 0 && w.params.verbosity >= 1){
            display =  std::vector<char>(51); // enable bound checking, no buffer overflow allowed
            display.back() = '\0';
            for (size_t i = 0; i < display.size()-1; i++) {
                display[i] = ' ';
            }
            itracker = 0;
            trackerinc = 1.0 / w.params.pts;
            progress = 0.00;            

        }
        // MAIN MCMC LOOP
        for (size_t i = 0; i < w.params.pts; i++) {
            // sweep each replica (here order doesn't matter)
            for (size_t s = 0; s < w.params.sweeps; s++){
                uint64_t replica_ind = 0;
                for (Replica& r: w.replicas) {
                    r.sweep(w.p);
                    // store energy values for use in replica exchange
                    w.energy_values_local[replica_ind++] = r.get_energy();
                }
            }
            
            // send all energy values to manager
            gather_energy_values(w);
            
            // perform Nt - 1 replica exchange proposals for each chain
            for (size_t m = 0; m < w.params.Nt - 1; m++) {
                if (w.wid == 0){
                    try_replica_exchange(w, 0); // chain 0 proposal
                    try_replica_exchange(w, 1); // chain 1 proposal
                } else {
                    check_for_replica_exchange(w); // check for exchange in either chain
                    check_for_replica_exchange(w);
                }
            }

            
            // if ICM is enabled
            if (w.params.icm) {
                // perform ICM on all relevant replicas

                for (int icm_target = w.params.Nt-1; icm_target >= static_cast<int>(w.params.Nt - w.Nc); icm_target--) {
                    if (w.wid == 0)
                    {
                        icm(w, icm_target);
                    }
                    else
                    {
                        check_icm(w);
                    }
                }
            }
            
            // if verbosity set to 1 or above, add nice colored progress bar
            if (w.wid == 0 && w.params.verbosity >= 1){
                itracker+=trackerinc;
                bool changed = false;
                while ((progress + 0.021) <= itracker) {
                    progress += 0.02;
                    changed = true;
                    display.at(display_index++) = '#';
                }
                if (changed) {
                    printf("\r%s%ld%c Done [%s%s%s]", display_index >= 5 ? " " : "  ",
                                                      display_index * 2, 
                                                      '%', 
                                                      GREEN, 
                                                      display.data(), 
                                                      RESET);
                    fflush (stdout);
                }
            }
        }
        if (w.wid == 0 && w.params.verbosity >= 1){
            display.at(display_index++) = '#';
            printf("\r%ld%c Done [%s%s%s]\n", display_index * 2, '%', GREEN, display.data(), RESET);
        }
        write_update(w, "Exiting MCMC loop", 1);
        MPI_Barrier(MPI_COMM_WORLD);
        best_ene = __DBL_MAX__;
        gather_best_energy_values(w, best_ene, best_state);
        // std::vector<double> final_energies (w.energy_values_local, w.energy_values_local + w.loc_r);
        flips = get_flips(w);
        if (w.wid == 0) {
            std::vector<double> swap_rates = w.swap_successes / w.swap_attempts;
            std::cout << swap_rates << std::endl;
            std::cout << w.swap_attempts << std::endl;
            std::cout << w.swap_successes << std::endl;

            std::cout << w.p.energy(best_state) << std::endl;
            std::cout << best_ene << " " << w.p.maxcut(best_state) << " " << flips << " " << w.swaps << " " << w.params.seed << std::endl;
            std::cout << best_state << std::endl;
            std::cout << std::vector<double> (w.beta_values, w.beta_values + w.params.Nt) << "\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    };
} // namespace pticm