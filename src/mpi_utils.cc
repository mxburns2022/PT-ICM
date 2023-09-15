#include "pt_icm.hh"

namespace pticm {
    MPI_Status stat;
    // linearly space the temperature values
    void beta_linspace(double* temps, double T0, double T1, uint64_t num_temp, bool inverse) {
        double step, curr_temp;
        // assume T0/T1 given as inverse temperatures
        if (inverse) {
            step = T1 < T0 ? (T0 - T1) / (num_temp - 1) : (T1 - T0) / (num_temp - 1);
            curr_temp = T1 < T0 ? T1 : T0;
        } else {
            // assume T0/T1 given as temperatures
            step = T1 < T0 ? (1 / T1 - 1 / T0) / (num_temp - 1) : (1 / T0 - 1 / T1) / (num_temp - 1);
            curr_temp = T1 < T0 ? 1/T0 : 1/T1;
        }
        for (uint64_t i = 0; i < num_temp; i++) {
            temps[i] = curr_temp;
            curr_temp += step;
        }
        // std::sort(temps, temps+num_temp); // ensure ordered low to high
    }

    void create_datatype(MPI_Datatype& param_type) {
        // Create a datatype to store parameter information
        const int nitems=2;
        int          param_blocklengths[] = {1,
                                             1,
                                             1,
                                             1,
                                             1,
                                             1,
                                             1,
                                             1,
                                             1,
                                             1,
                                             1,
                                             1,
                                             1};
        MPI_Datatype param_types[] = {MPI_UNSIGNED_LONG, // seed
                                    MPI_UNSIGNED_LONG, // Nt
                                    MPI_UNSIGNED_LONG, // pts
                                    MPI_UNSIGNED_LONG, // burn_in
                                    MPI_UNSIGNED_LONG, // Sweeps
                                    MPI_UNSIGNED, // verbosity
                                    MPI_DOUBLE,        // T0
                                    MPI_DOUBLE,        // T1
                                    MPI_DOUBLE,        // cutoff
                                    MPI_CXX_BOOL,      // inverse
                                    MPI_CXX_BOOL,      // preprocess
                                    MPI_CXX_BOOL,      // log_data
                                    MPI_CXX_BOOL      // icm
                            };
        MPI_Aint     param_offsets[13];
        uint32_t runner = 0;
        param_offsets[runner++] = offsetof(Parameters, seed);
        param_offsets[runner++] = offsetof(Parameters, Nt);
        param_offsets[runner++] = offsetof(Parameters, pts);
        param_offsets[runner++] = offsetof(Parameters, burn_in);
        param_offsets[runner++] = offsetof(Parameters, sweeps);
        param_offsets[runner++] = offsetof(Parameters, verbosity);
        param_offsets[runner++] = offsetof(Parameters, T0);
        param_offsets[runner++] = offsetof(Parameters, T1);
        param_offsets[runner++] = offsetof(Parameters, cutoff);
        param_offsets[runner++] = offsetof(Parameters, inverse);
        param_offsets[runner++] = offsetof(Parameters, preprocess);
        param_offsets[runner++] = offsetof(Parameters, log_data);
        param_offsets[runner++] = offsetof(Parameters, icm);
        MPI_Type_create_struct(runner,
                                param_blocklengths,
                                param_offsets,
                                 param_types,
                                 &param_type);
        MPI_Type_commit(&param_type);
    }; 
    // get the current energy values from each replica across all workers
    void gather_energy_values(PTWorker& w) {
        if (w.wid == 0) {
            for (uint64_t worker = 1; worker < w.numw - 1; worker++) {
                MPI_Recv(w.energy_values_global + w.worker_r * worker,
                        w.worker_r, 
                        MPI_DOUBLE, 
                        worker, 
                        0, 
                        MPI_COMM_WORLD, 
                        &w.stat);
            }
            if (w.numw > 1)
                MPI_Recv(w.energy_values_global + w.worker_r * (w.numw - 1),
                        w.last_r, 
                        MPI_DOUBLE, 
                        w.numw - 1, 
                        0, 
                        MPI_COMM_WORLD, 
                        &w.stat);
            std::memcpy(w.energy_values_global, w.energy_values_local, w.loc_r * sizeof(double));
        } else {
            MPI_Send(w.energy_values_local, w.loc_r, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

    };
    // get the best recorded energy value/state across all replicas
    void gather_best_energy_values(PTWorker& w, double& best_ene_global, StateVector& best_state) {
        double best_ene_local = __DBL_MAX__;
        StateVector state;
        size_t ind = 0;
        for (Replica& r: w.replicas) {
            double r_ene = r.get_best_ene();
            if (r_ene < best_ene_local) {
                best_ene_local = r_ene;
                state = r.get_best_state();
            }
            ind ++;
        }
        std::unique_ptr<bool[]> state_temp = std::make_unique<bool[]>(w.p.vcount());
        if (w.wid == 0) {
            best_ene_global = best_ene_local;
            best_state = state;
            double ene_temp;
            for (uint64_t worker = 1; worker < w.numw; worker++) {
                MPI_Recv(&ene_temp,
                        1, 
                        MPI_DOUBLE, 
                        worker, 
                        0, 
                        MPI_COMM_WORLD, 
                        &w.stat);
                MPI_Recv(state_temp.get(),
                        w.p.vcount(), 
                        MPI_CXX_BOOL, 
                        worker, 
                        0, 
                        MPI_COMM_WORLD, 
                        &w.stat);
                if (ene_temp < best_ene_global) {
                    best_state = StateVector(state_temp.get(), state_temp.get() + w.p.vcount());
                    best_ene_global = ene_temp;
                }
            }
        } else {
            for (size_t i = 0; i < w.p.vcount(); i++) {
                state_temp[i] = state[i];
            }
            MPI_Send(&best_ene_local, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(state_temp.get(), w.p.vcount(), MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD);
        }

    };


    // get the flip counts from each replica
    uint64_t get_flips(PTWorker& w) {
        uint64_t flip_sum = std::accumulate(w.replicas.begin(), w.replicas.end(), 0, [] (uint64_t sum, const Replica& r) {return sum + r.get_flips();}); 
        uint64_t recv_buff = 0;   
        std::unique_ptr<uint64_t[]> sumvals;
        if (w.wid == 0) {
            std::unique_ptr<uint64_t[]> sumvals = std::make_unique<uint64_t[]>(w.params.Nt * 2);
        }         
        MPI_Reduce(&flip_sum, &recv_buff, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        return recv_buff;   
    };
    // ensure that the beta indexing scheme is correct
    void check_correct_beta_values(PTWorker& w) {
        bool correct = 1;
        if (w.wid == 0) {
            correct = 1;
            std::unique_ptr<double[]> true_beta_values = std::make_unique<double[]>(w.params.Nt * 2);
            int i;
            for (i = 1; i < w.numw-1; i++) {
                MPI_Recv(true_beta_values.get()+(i * w.worker_r), w.worker_r, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &w.stat );
            }
            if (w.numw > 1) {
                MPI_Recv(true_beta_values.get()+(i * w.worker_r), w.last_r, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &w.stat );
            }
            std::memcpy(true_beta_values.get(), w.beta_values_local, w.loc_r * sizeof(double));
            std::vector<double> true_beta_vec (w.beta_values, w.beta_values + w.params.Nt*2);
            std::vector<double> rec_beta_vec (w.params.Nt*2);
            for (i = 0; i < w.params.Nt*2; i++) {
                uint32_t recorded_ind = w.beta_indices[i];
                rec_beta_vec[i] = true_beta_values[recorded_ind];
            }
            for (i = 0; i < w.params.Nt*2; i++) {

                if (std::abs(rec_beta_vec[i] - true_beta_vec[i]) > 1e-10) {
                    std::cout << "Beta tracking error " << i <<" " << rec_beta_vec[i] << " " << true_beta_vec[i] << "\n";
                    correct = false;
                    break;
                }
            }
        } else {
            MPI_Send(w.beta_values_local, w.loc_r, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        MPI_Bcast(&correct, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        if (!correct) {
            MPI_Finalize();
            exit(1);
        }
    }
    // get replica based on the beta value
    Replica& get_replica(PTWorker& w, double beta, uint8_t chain) {
        for (size_t i = 0; i < w.loc_r; i++) {
            Replica& r = w.replicas[i];
            if (std::abs(beta - r.get_beta()) < 1e-10 && chain == w.chain_index_local[i]) {
                return r;
            }
        }
    
        printf("Worker %d %d %d %f %d\n", w.wid, w.loc_r,  chain, beta, __LINE__);
        fflush(stdout);
        exit(-1);
    }

    // free worker memory
    void free_worker (PTWorker& w){
        if (w.wid == 0) {
            
            delete [] w.beta_values;
            
            delete [] w.energy_values_global;
            
            delete [] w.chain_1_exchange_order;
            
            delete [] w.chain_2_exchange_order;
            
            delete [] w.beta_indices;
            
            delete [] w.replica_workers;
            
        }
        delete [] w.energy_values_local;
        delete [] w.beta_values_local;
    };
} // namespace pticm