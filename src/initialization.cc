#include "pt_icm.hh"
#include "mpi_utils.hh"
#include "utils.hh"
namespace pticm {
    void read_problem(PTWorker& w, std::string problem_path) {
        size_t M, N;
        bool is_sparse;
        std::unique_ptr<double[]> weights; // arrays to store J matrix info
        std::unique_ptr<int[]> src;
        std::unique_ptr<int[]> dst;
        if (w.wid == 0) {
            Problem& p = w.p;
            
            p = Problem(problem_path,
                        w.params.seed, 
                        w.params.burn_in);
            w.characteristic_scale = p.get_char_scale();
            N = p.vcount();
            M = p.ecount();
            is_sparse = p.is_sparse();
        }
        MPI_Bcast(&is_sparse, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        if (is_sparse) {
            weights = std::make_unique<double[]>(M);
            src = std::make_unique<int[]>(M);
            dst = std::make_unique<int[]>(M);
            if (w.wid == 0) {
                w.p.serialize(weights.get(), src.get(), dst.get());
            }
            MPI_Bcast(weights.get(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(src.get(), M, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(dst.get(), M, MPI_INT, 0, MPI_COMM_WORLD);
        } else {
            weights = std::make_unique<double[]>(N * N);
            if (w.wid == 0) {
                w.p.serialize(weights.get(), NULL, NULL);
            }
            MPI_Bcast(weights.get(), N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if (w.wid != 0){
            if (is_sparse) {
                w.p = Problem(weights.get(), src.get(), dst.get(), w.params.seed + w.wid, w.params.burn_in, N, M);
            } else {
                w.p = Problem(weights.get(), w.params.seed + w.wid, w.params.burn_in, N);
            }
        }
    };

    void initialize_replicas(PTWorker& w) {
        // std::unique_ptr<bool[]> init_state = std::make_unique<bool[]>(w.p.vcount());
        // if (w.wid == 0) {

        std::bernoulli_distribution init_gen(0.5);
        // }
        // MPI_Bcast(init_state.get(), w.p.vcount(), MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        w.replicas = std::vector<Replica>(w.loc_r);
        StateVector init_vector (w.p.vcount());
        fflush(stdout);
        for (size_t i = 0; i < w.loc_r; i++) {
            for (uint64_t j = 0; j < w.p.vcount(); j++) {
                init_vector[j] = init_gen(w.p.get_rng());
            }
            std::vector<double> dE(w.p.vcount(), 0);
            double init_ene = w.p.calculate_dE(dE, init_vector); 
            w.replicas[i] = Replica(w.p, w.beta_values_local[i], init_vector, dE, init_ene);
        }
    };
    void initialize_replica_exchanges(PTWorker& w) {
        if (w.wid == 0) {
            size_t Nt = w.params.Nt;

            w.energy_values_global = new double [Nt * 2];
            w.chain_1_exchange_order = new uint32_t [Nt - 1];
            w.chain_2_exchange_order = new uint32_t [Nt - 1];
            w.beta_indices = new uint32_t [Nt * 2];
            w.swap_attempts = std::vector<uint32_t> (Nt * 2, 0);
            w.swap_successes = std::vector<uint32_t> (Nt * 2, 0);
            for (uint64_t i = 0; i < Nt - 1; i++) {
                w.chain_1_exchange_order[i] = i;
                w.chain_2_exchange_order[i] = i;
            }
            for (uint64_t i = 0; i < Nt * 2; i++) {
                w.beta_indices[i] = i;
            }
            w.exchange_index = 0;
            std::shuffle(w.chain_1_exchange_order, w.chain_1_exchange_order + Nt - 1, w.p.get_rng());
            std::shuffle(w.chain_2_exchange_order, w.chain_2_exchange_order + Nt - 1, w.p.get_rng());
        }
        w.energy_values_local = new double [w.loc_r];
    };

    // initialize worker data
    void initialize_temperatures(PTWorker& w) {
        Parameters& par = w.params;
        if (w.wid == 0) {
            w.beta_values = new double[par.Nt * 2];
            if (w.params.preprocess) {
                get_beta(w, 2, 1e-1, 0.3, 5000);
            } else {
                beta_linspace(w.beta_values, w.params.T0, w.params.T1, w.params.Nt, w.params.inverse);
            }
        }
        MPI_Bcast((void*)&w.params.Nt, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        uint64_t Nt = w.params.Nt;
        // Number of replicas to be handled by each worker
        if (Nt * 2 < w.numw) {
            if (w.wid == 0)
                fprintf(stderr,
                "%sERROR:%s %d workers created for %ld total replicas. Please call with fewer processes\n", RED, RESET,w.numw, Nt * 2);
            MPI_Finalize();
            exit(0);
        } 
        w.worker_r = Nt * 2 / w.numw;
        uint64_t regular = w.worker_r * w.numw;
        w.last_r = w.worker_r;
        if ((Nt * 2) != regular) {
            w.last_r += (Nt * 2) - regular; 
        }
        w.beta_values_local = new double[w.loc_r];
        w.chain_index_local = new uint8_t[w.loc_r];
        // std::cout <<  w.wid << " " << w.loc_r << std::endl;
        // getchar();
        // local problem declaration
        w.loc_r = w.wid != w.numw - 1 ? w.worker_r : w.last_r;
        if (w.wid == 0) {
            w.swaps = 0;
            w.replica_workers = new int[par.Nt * 2];
            w.chain_index_global = new uint8_t [par.Nt * 2];
            uint64_t ind = 0;
            for (uint64_t i = 0; i < w.numw - 1; i++) {
                for (uint64_t j = 0; j < w.worker_r; j++, ind++) {
                    w.replica_workers[ind] = i;
                    w.chain_index_global[ind] = ind >= par.Nt; 
                }
            }
            for (uint64_t j = 0; j < w.last_r; j++, ind++) {
                w.replica_workers[ind] = w.numw - 1;
                w.chain_index_global[ind] = ind >= par.Nt; 
            }
            // Determine the starting beta values for each replica

            // beta_linspace(w.beta_values, w.params.T0, w.params.T1, par.Nt, w.params.inverse);
            w.Nc = 0;
            for (size_t i = 0; i < par.Nt; i++) {
                if (w.beta_values[i] > w.characteristic_scale) {
                    break;
                }
                w.Nc ++;
            }
            std::memcpy(w.beta_values + par.Nt, w.beta_values, par.Nt * sizeof(double));
            // copy to the beta values for worker 0
            std::memcpy(w.chain_index_local, w.chain_index_global, w.loc_r * sizeof(uint8_t));
            std::memcpy(w.beta_values_local, w.beta_values, w.loc_r * sizeof(double));
            // send local beta values to each other worker
            for (uint64_t i = 1; i < w.numw - 1; i++) {
                double* start_pointer = i * w.worker_r + w.beta_values;
                MPI_Send(start_pointer, w.worker_r, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                uint8_t* start_chain_pointer = i * w.worker_r + w.chain_index_global;
                MPI_Send(start_chain_pointer, w.worker_r, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
            }
            if (w.numw > 1) {
                // send to final
                double* start_pointer = (w.numw - 1) * w.worker_r + w.beta_values;
                MPI_Send(start_pointer, w.last_r, MPI_DOUBLE, w.numw - 1, 0, MPI_COMM_WORLD);
                uint8_t* start_chain_pointer = (w.numw - 1) * w.worker_r + w.chain_index_global;
                MPI_Send(start_chain_pointer, w.last_r, MPI_UNSIGNED_CHAR, w.numw - 1, 0, MPI_COMM_WORLD);
            }
        } else {
            MPI_Recv(w.beta_values_local, w.loc_r, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &w.stat);
            MPI_Recv(w.chain_index_local, w.loc_r, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &w.stat);
        }
        MPI_Bcast(&w.Nc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(&w.params.icm, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    };
} // namespace pticm