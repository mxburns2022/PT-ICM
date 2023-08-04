#include "solver.hh"
#include "problem.hh"
#include "icm.hh"
#include "replica_exchange.hh"
#include "mpi_utils.hh"
#include "utils.hh"
#include "initialization.hh"
#include "definitions.hh"

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
        initialize_temperatures(w);
        initialize_replicas(w);
        initialize_replica_exchanges(w);
        //check_correct_beta_values(w);
        std::cout << "Worker " << w.wid << " entering main loop\n";
        std::cout << w.wid << " " << w.Nc << " " << w.params.icm << std::endl; 
        for (size_t i = 0; i < w.params.pts; i++) {

            uint64_t replica_ind = 0;
            for (Replica& r: w.replicas) {
                r.sweep(w.p);
                w.energy_values_local[replica_ind++] = r.get_energy();
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
        }
        MPI_Barrier(MPI_COMM_WORLD);
        double best_ene = __DBL_MAX__;
        StateVector best_state;
        gather_best_energy_values(w, best_ene, best_state);
        // std::vector<double> final_energies (w.energy_values_local, w.energy_values_local + w.loc_r);
        uint64_t flips = get_flips(w);
        if (w.wid == 0) {
            std::cout << w.p.energy(best_state) << std::endl;
            std::cout << best_ene << " " << w.p.maxcut(best_state) << " " << flips << " " << w.swaps << " " << w.params.seed << std::endl;
            std::cout << best_state << std::endl;
            std::cout << std::vector<double> (w.beta_values, w.beta_values + w.params.Nt) << "\n";
        }
        free_worker(w);

    };
} // namespace pticm