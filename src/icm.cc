#include "icm.hh"
#include "mpi_utils.hh"
#include <stack>

namespace pticm {
    int64_t xor_icm(PTWorker& w, bool* state_1, bool* state_2) {
        size_t Nt = w.params.Nt;
        std::vector<int64_t> oneinds;
        oneinds.reserve(w.p.vcount());
        for (int64_t j = 0; j < static_cast<int64_t>(w.p.vcount()); j++) {
            size_t val = state_1[j] ^ state_2[j];
            if (val == 1) {
                oneinds.push_back(j);
            }
        }
        size_t cluster_sites = oneinds.size();
        // Flip all states in R1 if cluster site counts are sufficiently large
        if (cluster_sites > w.p.vcount() / 2) {
            oneinds.clear();
            for (size_t j = 0; j < w.p.vcount(); j++) {
                state_1[j] = !state_1[j];
                if (state_1[j] != state_2[j]) {
                    oneinds.push_back(static_cast<int64_t>(j));
                }
            }
        }
        std::vector<int64_t> seed_val(1);
        if (oneinds.size() == 0) {
            return -1;
        }
        std::sample(oneinds.begin(), oneinds.end(), seed_val.begin(), 1, w.p.get_rng());
        return seed_val.front();
    }
    void bfs_icm(Problem& p, int64_t seed_val, bool* state_1, bool* state_2) {
        std::stack<int64_t> ind_stack;
        ind_stack.push(seed_val);
        std::vector<bool> added(p.vcount(), false);

        size_t numflipped = 0;
        added[seed_val] = true;
        while (!ind_stack.empty()) {
            int64_t ind = ind_stack.top();
            ind_stack.pop();

            state_1[ind] = !state_1[ind];
            state_2[ind] = !state_2[ind];

            numflipped++;
            if (p.is_sparse()) {
                const std::list<Edge>& elist = p.get_sparse_neighbors(ind);
                for (Edge e : elist) {
                    if (state_1[e.dest] != state_2[e.dest] && !added.at(e.dest)) {
                        added[e.dest] = true;
                        ind_stack.push(static_cast<int64_t>(e.dest));
                    }
                }
            }
            else {
                const std::vector<double>& elist = p.get_dense_neighbors(ind);
                for (size_t i = 0; i < p.vcount(); i++) {
                    if (std::abs(elist[i]) > 1e-10 && state_1[i] != state_2[i] && !added[i]) {
                        ind_stack.push(static_cast<int64_t>(i));
                        added[i] = true;
                    }
                }
            }
        }

        // if(std::abs(start_ene - (r1_ene + r2_ene)) > 1e-10){
        //     printf("ERROR, cluster move not isoenergetic\n\t%f\t%f\t%f\n", start_ene, r1_ene + r2_ene, p.energy(state_1)+p.energy(state_2));
        //     exit(-1);
        // }
        // std::cout << r1.get_beta() << " " << r2.get_beta() << " " << static_cast<double>(numflipped) / p.vcount() << std::endl;
    }

    void icm(PTWorker& w, size_t target_replica) {
        // get the indices for the replicas with equal temperatures in the two chains
        uint32_t replica_1 = w.beta_indices[target_replica];
        uint32_t replica_2 = w.beta_indices[target_replica + w.params.Nt];
        double beta_val = w.beta_values[target_replica];
        uint32_t worker_1 = w.replica_workers[replica_1];
        uint32_t worker_2 = w.replica_workers[replica_2];
        int8_t nack = -1;

        std::unique_ptr<bool[]> state_1 = std::make_unique<bool[]>(w.p.vcount());
        std::unique_ptr<bool[]> state_2 = std::make_unique<bool[]>(w.p.vcount());
        // Send acknowledgements to the workers holding the replicas
        // ack = 0/1 means only expecting one replica (from chain 0 or 1), ack = 2 means expecting 2
        for (int i = 1; i < w.numw; i++) {
            // Send negative acknowledgements to all other workers
            if (i != worker_1 && i != worker_2)
            {
                MPI_Send(&nack, 1, MPI_SIGNED_CHAR, i, 0, MPI_COMM_WORLD);
            }
        }
        // return;
        if (worker_1 != worker_2 && worker_1 != 0 && worker_2 != 0) {
            int8_t ack = 2;
            MPI_Send(&ack, 1, MPI_SIGNED_CHAR, worker_1, 0, MPI_COMM_WORLD);
            MPI_Send(&beta_val, 1, MPI_DOUBLE, worker_1, 0, MPI_COMM_WORLD);
            MPI_Send(&worker_2, 1, MPI_UNSIGNED, worker_1, 0, MPI_COMM_WORLD);

            ack = 3;
            MPI_Send(&ack, 1, MPI_SIGNED_CHAR, worker_2, 0, MPI_COMM_WORLD);
            MPI_Send(&beta_val, 1, MPI_DOUBLE, worker_2, 0, MPI_COMM_WORLD);
            MPI_Send(&worker_1, 1, MPI_UNSIGNED, worker_2, 0, MPI_COMM_WORLD);
            return;
            // MPI_Recv(state_1.get(), w.p.vcount(), MPI_CXX_BOOL, worker_1, 0, MPI_COMM_WORLD, &w.stat);
            // MPI_Recv(state_2.get(), w.p.vcount(), MPI_CXX_BOOL, worker_2, 0, MPI_COMM_WORLD, &w.stat);
        }
        else if (worker_1 == worker_2 && worker_2 != 0) {
            int8_t ack = 4;
            // the worker will perform the swap locally
            MPI_Send(&ack, 1, MPI_SIGNED_CHAR, worker_1, 0, MPI_COMM_WORLD);
            MPI_Send(&beta_val, 1, MPI_DOUBLE, worker_1, 0, MPI_COMM_WORLD);
            return;
        }
        else if (worker_1 != 0) {
            int8_t ack = 0;
            MPI_Send(&ack, 1, MPI_SIGNED_CHAR, worker_1, 0, MPI_COMM_WORLD);
            MPI_Send(&beta_val, 1, MPI_DOUBLE, worker_1, 0, MPI_COMM_WORLD);
            MPI_Recv(state_1.get(), w.p.vcount(), MPI_CXX_BOOL, worker_1, 0, MPI_COMM_WORLD, &w.stat);
            StateVector& s2 = get_replica(w, beta_val, 1).get_state();
            std::copy(s2.begin(), s2.end(), state_2.get());
        }
        else if (worker_2 != 0) {
            StateVector& s1 = get_replica(w, beta_val, 0).get_state();
            std::copy(s1.begin(), s1.end(), state_1.get());
            int8_t ack = 1;

            MPI_Send(&ack, 1, MPI_SIGNED_CHAR, worker_2, 0, MPI_COMM_WORLD);
            MPI_Send(&beta_val, 1, MPI_DOUBLE, worker_2, 0, MPI_COMM_WORLD);
            MPI_Recv(state_2.get(), w.p.vcount(), MPI_CXX_BOOL, worker_2, 0, MPI_COMM_WORLD, &w.stat);
        }
        else {
            StateVector& s1 = get_replica(w, beta_val, 0).get_state();
            std::copy(s1.begin(), s1.end(), state_1.get());

            StateVector& s2 = get_replica(w, beta_val, 1).get_state();
            std::copy(s2.begin(), s2.end(), state_2.get());
        }
        int64_t seed_val = xor_icm(w, state_1.get(), state_2.get());

        if (seed_val == -1) {
            return;
        }
        double ene1, ene2, ene1_2, ene2_2;
        if (w.params.verbosity >= 2) {
            ene1 = w.p.energy(std::vector<bool>(state_1.get(), state_1.get() + w.p.vcount()));
            ene2 = w.p.energy(std::vector<bool>(state_2.get(), state_2.get() + w.p.vcount()));

        }
        bfs_icm(w.p, seed_val, state_1.get(), state_2.get());
        if (w.params.verbosity >= 2) {
            ene1_2 = w.p.energy(std::vector<bool>(state_1.get(), state_1.get() + w.p.vcount()));
            ene2_2 = w.p.energy(std::vector<bool>(state_2.get(), state_2.get() + w.p.vcount()));

        }
        // if worker_1 == worker_2, it will have performed the ICM locally

        if (worker_1 != 0) {

            // printf("Sender %d %d\n", __LINE__, worker_2);
            MPI_Send(state_1.get(), w.p.vcount(), MPI_CXX_BOOL, worker_1, 0, MPI_COMM_WORLD);

            Replica& r2 = get_replica(w, beta_val, 1);
            StateVector s2 = std::vector<bool>(state_2.get(), state_2.get() + w.p.vcount());
            r2.set_state(s2);
            std::vector<double> dEtemp2(w.p.vcount());
            double energy = w.p.calculate_dE(dEtemp2, s2);
            //  if(std::abs(r2.get_energy() - w.p.energy(r2.get_state())) > 1e-10) {
            //     printf("Worker %d %d\n", w.wid, __LINE__);
            //     fflush(stdout);
            //     exit(-1);
            // }
        }
        else if (worker_2 != 0) {

            MPI_Send(state_2.get(), w.p.vcount(), MPI_CXX_BOOL, worker_2, 0, MPI_COMM_WORLD);
            Replica& r1 = get_replica(w, beta_val, 0);
            StateVector s1 = std::vector<bool>(state_1.get(), state_1.get() + w.p.vcount());
            r1.set_state(s1);
            std::vector<double> dEtemp1(w.p.vcount());
            double energy = w.p.calculate_dE(dEtemp1, s1);
            r1.set_energy(energy);
            r1.set_dE(dEtemp1);
        }
        else {
            Replica& r1 = get_replica(w, beta_val, 0);
            StateVector s1 = std::vector<bool>(state_1.get(), state_1.get() + w.p.vcount());
            r1.set_state(s1);
            std::vector<double> dEtemp1(w.p.vcount());
            double energy = w.p.calculate_dE(dEtemp1, s1);
            r1.set_energy(energy);
            r1.set_dE(dEtemp1);


            Replica& r2 = get_replica(w, beta_val, 1);
            StateVector s2 = std::vector<bool>(state_2.get(), state_2.get() + w.p.vcount());
            r2.set_state(s2);
            std::vector<double> dEtemp2(w.p.vcount());
            energy = w.p.calculate_dE(dEtemp2, s2);

            r2.set_energy(energy);
            r2.set_dE(dEtemp2);
            // if(std::abs(r1.get_energy() - w.p.energy(r1.get_state())) > 1e-10) {
            //     printf("Worker %d %d\n", w.wid, __LINE__);
            //     fflush(stdout);
            //     exit(-1);
            // }
            // if(std::abs(r2.get_energy() - w.p.energy(r2.get_state())) > 1e-10) {
            //     printf("Worker %d %d\n", w.wid, __LINE__);
            //     fflush(stdout);
            //     exit(-1);
            // }
        }


    }
    void check_icm(PTWorker& w) {
        int8_t incoming;
        MPI_Recv(&incoming, 1, MPI_SIGNED_CHAR, 0, 0, MPI_COMM_WORLD, &w.stat);
        if (incoming == -1) {
            return;
        }
        double beta_value;
        MPI_Recv(&beta_value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &w.stat);
        if (incoming < 2) {
            // copy the replica state
            Replica& r = get_replica(w, beta_value, incoming);
            std::unique_ptr<bool[]> sarr = std::make_unique<bool[]>(w.p.vcount());
            StateVector& s = r.get_state();
            std::copy(s.begin(), s.end(), sarr.get());
            // send the current state
            MPI_Send(sarr.get(), w.p.vcount(), MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD);
            // get the new state
            MPI_Recv(sarr.get(), w.p.vcount(), MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, &w.stat);
            StateVector s1 = std::vector<bool>(sarr.get(), sarr.get() + w.p.vcount());
            r.set_state(s1);
            std::vector<double> dEtemp1(w.p.vcount());
            double energy = w.p.calculate_dE(dEtemp1, s1);
            r.set_energy(energy);
            r.set_dE(dEtemp1);
            // if(std::abs(r.get_energy() - w.p.energy(r.get_state())) > 1e-10) {
            //     printf("Worker %d %d\n", w.wid, __LINE__);
            //     fflush(stdout);
            //     exit(-1);
            // }
        }
        else if (incoming == 2) {
            // send the data to index 0
            Replica& r = get_replica(w, beta_value, 0);
            std::unique_ptr<bool[]> sarr1 = std::make_unique<bool[]>(w.p.vcount());
            std::unique_ptr<bool[]> sarr2 = std::make_unique<bool[]>(w.p.vcount());
            StateVector& s1 = r.get_state();
            std::copy(s1.begin(), s1.end(), sarr1.get());
            uint32_t other_worker;
            MPI_Recv(&other_worker, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &w.stat);
            MPI_Recv(sarr2.get(), w.p.vcount(), MPI_CXX_BOOL, other_worker, 0, MPI_COMM_WORLD, &w.stat);

            int64_t seed_val = xor_icm(w, sarr1.get(), sarr2.get());

            MPI_Send(sarr2.get(), w.p.vcount(), MPI_CXX_BOOL, other_worker, 0, MPI_COMM_WORLD);
            StateVector s1new = std::vector<bool>(sarr1.get(), sarr1.get() + w.p.vcount());
            r.set_state(s1new);
            std::vector<double> dEtemp1(w.p.vcount());
            double energy = w.p.calculate_dE(dEtemp1, s1new);
            r.set_energy(energy);
            r.set_dE(dEtemp1);
            // if(std::abs(r.get_energy() - w.p.energy(r.get_state())) > 1e-10) {
            //     printf("Worker %d %d\n", w.wid, __LINE__);
            //     fflush(stdout);
            //     exit(-1);
            // }
        }
        else if (incoming == 3) {
            // get the rank of the other worker
            uint32_t other_worker;
            MPI_Recv(&other_worker, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &w.stat);
            // copy the replica state
            Replica& r = get_replica(w, beta_value, 1);
            std::unique_ptr<bool[]> sarr = std::make_unique<bool[]>(w.p.vcount());
            StateVector& s = r.get_state();
            std::copy(s.begin(), s.end(), sarr.get());
            // send the current state
            MPI_Send(sarr.get(), w.p.vcount(), MPI_CXX_BOOL, other_worker, 0, MPI_COMM_WORLD);
            // get the new state
            MPI_Recv(sarr.get(), w.p.vcount(), MPI_CXX_BOOL, other_worker, 0, MPI_COMM_WORLD, &w.stat);
            StateVector s1 = std::vector<bool>(sarr.get(), sarr.get() + w.p.vcount());
            r.set_state(s1);
            std::vector<double> dEtemp1(w.p.vcount());
            double energy = w.p.calculate_dE(dEtemp1, s1);
            r.set_energy(energy);
            r.set_dE(dEtemp1);
            // if(std::abs(r.get_energy() - w.p.energy(r.get_state())) > 1e-10) {
            //     printf("Worker %d %d\n", w.wid, __LINE__);
            //     fflush(stdout);
            //     exit(-1);
            // }

        }
        else if (incoming == 4) {
            // perform the ICM locally
            Replica& r1 = get_replica(w, beta_value, 0);
            std::unique_ptr<bool[]> sarr1 = std::make_unique<bool[]>(w.p.vcount());
            StateVector& s1 = r1.get_state();
            std::copy(s1.begin(), s1.end(), sarr1.get());
            Replica& r2 = get_replica(w, beta_value, 1);
            std::unique_ptr<bool[]> sarr2 = std::make_unique<bool[]>(w.p.vcount());
            StateVector& s2 = r2.get_state();
            std::copy(s2.begin(), s2.end(), sarr2.get());

            int64_t seed_val = xor_icm(w, sarr1.get(), sarr2.get());

            bfs_icm(w.p, seed_val, sarr1.get(), sarr2.get());

            StateVector s1new = std::vector<bool>(sarr1.get(), sarr1.get() + w.p.vcount());
            r1.set_state(s1new);
            std::vector<double> dEtemp1(w.p.vcount());
            double energy = w.p.calculate_dE(dEtemp1, s1new);
            r1.set_energy(energy);
            r1.set_dE(dEtemp1);

            StateVector s2new = std::vector<bool>(sarr2.get(), sarr2.get() + w.p.vcount());
            r2.set_state(s2new);
            std::vector<double> dEtemp2(w.p.vcount());
            energy = w.p.calculate_dE(dEtemp2, s2new);

            r2.set_energy(energy);
            r2.set_dE(dEtemp2);
            // if(std::abs(r1.get_energy() - w.p.energy(r1.get_state())) > 1e-10) {
            //     printf("Worker %d %d\n", w.wid, __LINE__);
            //     fflush(stdout);
            //     exit(-1);
            // }
            // if(std::abs(r2.get_energy() - w.p.energy(r2.get_state())) > 1e-10) {
            //     printf("Worker %d %d\n", w.wid, __LINE__);
            //     fflush(stdout);
            //     exit(-1);
            // }
        }

    }
} // namespace pticm