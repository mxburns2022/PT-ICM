#include "mpi_defs.hh"

MPI_Status stat;
void beta_linspace(double* temps, double T0, double T1, uint64_t num_temp, bool inverse) {
    double step = inverse ? (T1 - T0) / (num_temp - 1) : (1 / T1 - 1 / T0) / (num_temp - 1);
    double curr_temp = inverse ? T0 : 1 / T0;
    for (uint64_t i = 0; i < num_temp; i++) {
        temps[i] = curr_temp;
        curr_temp += step;
    }
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
                                         1};
    MPI_Datatype param_types[] = {MPI_UNSIGNED_LONG, // seed
                                  MPI_UNSIGNED_LONG, // Nt
                                  MPI_UNSIGNED_LONG, // pts
                                  MPI_UNSIGNED_LONG, // burn_in
                                  MPI_UNSIGNED, // verbosity
                                  MPI_DOUBLE,        // T0
                                  MPI_DOUBLE,        // T1
                                  MPI_CXX_BOOL,      // inverse
                                  MPI_CXX_BOOL      // log_data
                           };
    MPI_Aint     param_offsets[9];
    uint32_t runner = 0;
    param_offsets[runner++] = offsetof(Parameters, seed);
    param_offsets[runner++] = offsetof(Parameters, Nt);
    param_offsets[runner++] = offsetof(Parameters, pts);
    param_offsets[runner++] = offsetof(Parameters, burn_in);
    param_offsets[runner++] = offsetof(Parameters, verbosity);
    param_offsets[runner++] = offsetof(Parameters, T0);
    param_offsets[runner++] = offsetof(Parameters, T1);
    param_offsets[runner++] = offsetof(Parameters, inverse);
    param_offsets[runner++] = offsetof(Parameters, log_data);
    MPI_Type_create_struct(runner,
                           param_blocklengths,
                           param_offsets,
                           param_types,
                           &param_type);
    MPI_Type_commit(&param_type);
}; 
// initialize worker data
void initialize_temperatures(PTWorker& w) {
    Parameters& par = w.params;
    w.beta_values_local = new double[w.loc_r];
    w.chain_index_local = new uint8_t[w.loc_r];
    // std::cout <<  w.wid << " " << w.loc_r << std::endl;
    // getchar();
    // Distribute the replicas to each worker
    if (w.wid == 0) {
        std::cout << "Distributing replicas\n";
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
        w.beta_values = new double[par.Nt * 2];
        beta_linspace(w.beta_values, w.params.T0, w.params.T1, par.Nt, w.params.inverse);
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

        // std::cout << w.wid << " chain indices: " <<std::vector<uint8_t> (w.chain_index_global, w.chain_index_global + w.params.Nt * 2) << std::endl;
    } else {
        MPI_Recv(w.beta_values_local, w.loc_r, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &w.stat);
        MPI_Recv(w.chain_index_local, w.loc_r, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &w.stat);
        // std::cout << w.wid << " chain indices: " <<std::vector<uint8_t> (w.chain_index_local, w.chain_index_local + w.loc_r) << std::endl;
    }
};
void read_problem(PTWorker& w, std::string problem_path) {
    size_t M, N;
    bool is_sparse;
    std::unique_ptr<double[]> weights; // arrays to store J matrix info
    std::unique_ptr<int[]> src;
    std::unique_ptr<int[]> dst;
    if (w.wid == 0) {
        std::cout << "Constructing original problem... ";
        Problem& p = w.p;
        p = Problem(problem_path,
                    w.params.seed, 
                    w.params.burn_in);
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
void reconstruct_problem(PTWorker& w) {

};
void initialize_replicas(PTWorker& w) {
    std::unique_ptr<bool[]> init_state = std::make_unique<bool[]>(w.p.vcount());
    if (w.wid == 0) {

        std::bernoulli_distribution init_gen(0.5);
        for (uint64_t i = 0; i < w.p.vcount(); i++) {
            init_state[i] = init_gen(w.p.get_rng());
        }
    }
    MPI_Bcast(init_state.get(), w.p.vcount(), MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    std::vector<bool> init_vector (init_state.get(), init_state.get() + w.p.vcount());
    std::vector<double> dE(w.p.vcount());
    double init_ene = w.p.calculate_dE(dE, init_vector); 
    w.replicas = std::vector<Replica>(w.loc_r);
    for (size_t i = 0; i < w.loc_r; i++) {
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
    w.energy_values_local = new double[w.loc_r];
};
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
        // std::cout << std::vector<double>(w.energy_values_global, w.energy_values_global+w.params.Nt * 2) << std::endl; 
        // std::cout << std::vector<double>(w.beta_values, w.beta_values+w.params.Nt * 2) << std::endl; 
        // getchar();
    } else {
        MPI_Send(w.energy_values_local, w.loc_r, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

};
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
    }
    std::unique_ptr<bool[]> state_temp = std::make_unique<bool[]>(w.p.vcount());
    if (w.wid == 0) {
        // std::cout << best_ene_local << std::endl;
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
                // std::cout << best_ene_global << std::endl;
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
void sent_energy_values(PTWorker& w) {

};
void try_replica_exchange(PTWorker& w,
                          uint8_t chain_index) {
    assert(w.wid == 0);
    uint64_t Nt = w.params.Nt;
    uint64_t chain_exchange_index = w.exchange_index + (chain_index > 0 ? Nt : 0);
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
    if (w.params.verbosity >= 2) {
        std::cout << "\n\n(%d %d)\ndB: " << dBeta << std::endl;
        std::cout << "dE: " << dE << std::endl;
        std::cout << "B1: " << beta_low << " B2:" << beta_high  << std::endl;
        std::cout << "I1: " << replica_1_ind << " I2:" << replica_2_ind  << std::endl;
        std::cout << "E1: " << w.energy_values_global[replica_1_ind] << " E2:" << w.energy_values_global[replica_2_ind]  << std::endl;
        std::cout << "log(dB*dE): " <<  log_boltzmann << std::endl;
        std::cout<< "RND: " <<  rndcomp << std::endl;
    }
    // getchar();
    if (rndcomp > log_boltzmann) {
        if (w.params.verbosity >= 2) {
            printf("\nSWAPPING temperatures %f <-> %f in replicas %d and %d\n", beta_low, beta_high, worker_1, worker_2); 
            fflush(stdout);
        }
        w.swaps ++;
        bool swap = true;
        if (worker_1 == 0) {
            do_replica_exchange_local(w, beta_low, beta_high, chain_index);
        } else {
            // send the old and new beta values
            MPI_Send(&swap, 1, MPI_CXX_BOOL, worker_1, 0, MPI_COMM_WORLD);
            MPI_Send(&beta_high, 1, MPI_DOUBLE, worker_1, 0, MPI_COMM_WORLD);
            MPI_Send(&beta_low, 1, MPI_DOUBLE, worker_1, 0, MPI_COMM_WORLD);
            MPI_Send(&chain_index, 1, MPI_UNSIGNED_CHAR, worker_1, 0, MPI_COMM_WORLD);
        }
        // send info to worker with the first replica
        // send info to worker with the second replica
        if (worker_1 != worker_2) {
            if (worker_2 == 0) {
                do_replica_exchange_local(w, beta_high, beta_low, chain_index);
            } else {
                // send the old and new beta values
                MPI_Send(&swap, 1, MPI_CXX_BOOL, worker_2, 0, MPI_COMM_WORLD);
                MPI_Send(&beta_low, 1, MPI_DOUBLE, worker_2, 0, MPI_COMM_WORLD);
                MPI_Send(&beta_high, 1, MPI_DOUBLE, worker_2, 0, MPI_COMM_WORLD);
                MPI_Send(&chain_index, 1, MPI_UNSIGNED_CHAR, worker_2, 0, MPI_COMM_WORLD);
            }
        }
        std::swap(w.beta_indices[chain_exchange_index], w.beta_indices[chain_exchange_index + 1]);
    } else {
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
void do_replica_exchange_local(PTWorker& w, 
                               double beta_old, 
                               double beta_new,
                               uint8_t chain_index) {
    //  std::cout << w.wid << " " << beta_old << " " << beta_new << std::endl;
    //  std::cout << std::vector<double>(w.beta_values_local, w.beta_values_local + w.loc_r) << "\n"; 
     bool check = false;
     for (uint64_t j = 0; j < w.loc_r; j++) {
        // printf("(%d, %f) ", w.chain_index_local[j], w.beta_values_local[j]);
        if (std::abs(w.beta_values_local[j] - beta_new) < 1e-10 && w.chain_index_local[j] == chain_index) {
            Replica& r = w.replicas[j];
            // printf("\nWorker %d replacing %f (%f) with %f\n", w.wid, r.get_beta(), w.beta_values_local[j], beta_old);
            w.beta_values_local[j] = beta_old;
            r.set_beta(beta_old);
            check = true;
        }
        else if (std::abs(w.beta_values_local[j] - beta_old) < 1e-10 && w.chain_index_local[j] == chain_index) {
            Replica& r = w.replicas[j];
            // printf("\nWorker %d replacing %f (%f) with %f\n", w.wid, r.get_beta(), w.beta_values_local[j], beta_new);
            w.beta_values_local[j] = beta_new;
            r.set_beta(beta_new);
            check = true;
        }
    }   
    assert(check);
};
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
void check_correct_beta_values(PTWorker& w) {
    bool correct = 1;
    if (w.wid == 0) {
        // std::cout << "\nSTARTING CHECK\n";
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
        // std::cout << rec_beta_vec << std::endl << std::endl;
        // std::cout << true_beta_vec << std::endl;
        for (i = 0; i < w.params.Nt*2; i++) {

            if (std::abs(rec_beta_vec[i] - true_beta_vec[i]) > 1e-10) {
                std::cout << "Beta tracking error " << i <<" " << rec_beta_vec[i] << " " << true_beta_vec[i] << "\n";
                correct = false;
                break;
            }
        }
        // std::cout << "DONE CHECK\n";
    } else {
        MPI_Send(w.beta_values_local, w.loc_r, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    MPI_Bcast(&correct, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    if (!correct) {
        std::cout << correct << std::endl;
        MPI_Finalize();
        exit(1);
    }

}
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
    initialize_temperatures(w);
    read_problem(w, problem_path);
    initialize_replicas(w);
    initialize_replica_exchanges(w);
    //check_correct_beta_values(w);
    std::cout << "Worker " << w.wid << " entering main loop\n";
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
            check_correct_beta_values(w);
        }
    }

    double best_ene = __DBL_MAX__;
    StateVector best_state;
    gather_best_energy_values(w, best_ene, best_state);
    // std::vector<double> final_energies (w.energy_values_local, w.energy_values_local + w.loc_r);
    uint64_t flips = get_flips(w);
    if (w.wid == 0) {
       std::cout << best_ene << " " << w.p.maxcut(best_state) << " " << flips << " " << w.swaps << std::endl;
       std::cout << best_state << std::endl;
    }
    free_worker(w);

};
void free_worker (PTWorker& w){
    if (w.wid == 0) {
        delete w.beta_values;
        delete w.energy_values_global;
        delete w.chain_1_exchange_order;
        delete w.chain_2_exchange_order;
        delete w.beta_indices;
        delete w.replica_workers;
    }
    delete w.energy_values_local;
    delete w.beta_values_local;

};