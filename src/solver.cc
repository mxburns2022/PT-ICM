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
    void solve(PTWorker& w, std::string& problem_path) {
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
        // local problem declaration
        w.loc_r = w.wid != w.numw - 1 ? w.worker_r : w.last_r;
        
        read_problem(w, problem_path);
        write_update(w, "Initializing Temperatures", 1);
        initialize_temperatures(w);
        write_update(w, "Initializing Replicas", 1);
        initialize_replicas(w);
        write_update(w, "Initializing Exchanges", 1);
        initialize_replica_exchanges(w);
        write_update(w, "Starting MCMC loop", 1);
        MPI_Barrier(MPI_COMM_WORLD);
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
        for (size_t i = 0; i < w.params.pts; i++) {

            for (size_t s = 0; s < w.params.sweeps; s++){
                uint64_t replica_ind = 0;
                for (Replica& r: w.replicas) {
                    r.sweep(w.p);
                    w.energy_values_local[replica_ind++] = r.get_energy();
                }
            }
            gather_energy_values(w);
            for (size_t m = 0; m < w.params.Nt - 1; m++) {
                if (w.wid == 0){
                    try_replica_exchange(w, 0);
                    try_replica_exchange(w, 1);
                } else {

                    check_for_replica_exchange(w);
                    check_for_replica_exchange(w);
                }
            }
            if (w.params.icm) {
                for (size_t icm_target = 0; icm_target < w.Nc; icm_target++) {
                    if (w.wid == 0) {
                        icm(w, icm_target);
                    } else {
                        check_icm(w);
                    }
                }
            }
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
        double best_ene = __DBL_MAX__;
        StateVector best_state;
        gather_best_energy_values(w, best_ene, best_state);
        // std::vector<double> final_energies (w.energy_values_local, w.energy_values_local + w.loc_r);
        uint64_t flips = get_flips(w);
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
        free_worker(w);

    };
} // namespace pticm