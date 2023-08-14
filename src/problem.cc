#include "pt_icm.hh"

#include "utils.hh"
namespace pticm {
    void Problem::serialize(double* weights, int* src, int* dest) {
        if (sparse) {
            uint64_t ind = 0;
            int srcnode = 0;
            for (auto& i : formulation_sparse) {
                for (Edge e: i) {
                    if (srcnode < e.dest){
                        weights[ind] = e.weight;
                        src[ind] = srcnode;
                        dest[ind] = e.dest;
                        ind++;
                        
                    }
                }
                srcnode++;
            }
            assert(ind == M);
        } else {
            printf("DENSE\n");
            double* curr_val = weights;
            for (auto& i : formulation) {
                std::memcpy(curr_val, i.data(), N * sizeof(double));
                curr_val += N;
            }
            getchar();
        }
    };

    double Problem::calculate_dE (std::vector<double>& dE, const StateVector& s) const {
        double ene = 0.0;
        double ene_term;
        assert(sign(false) == -1);
        assert(sign(true) == 1); 
        for (uint64_t i = 0; i < N; i++) {
            int ival = sign(s[i]);
            if (sparse) {
                const std::list<Edge>& iedges = formulation_sparse.at(i);
                for (auto e: iedges) {
                    ene_term = ival * sign(s[e.dest]) * e.weight;
                    ene += ene_term / 2;
                    dE[i] -= 2 * ene_term;
                }
            } else {
                const std::vector<double>& iedges = formulation.at(i);
                for (uint64_t j = 0; j < N; j++) {
                    if (i == j) {
                        continue;
                    }
                    ene_term = ival * sign(s[j]) * iedges.at(j);
                    ene += ene_term / 2;
                    dE[i] -= 2 * ene_term;
                }
            }
        }
    // #ifndef NDEBUG
    //     StateVector state = s;

    //     for (uint64_t i = 0; i < N; i++) {
    //         state[i] = !state[i];
    //         double ene_new = energy(state);
    //         double actual_de = ene_new - ene;
    //         if(std::abs(actual_de - dE[i]) > 1e-10) {
    //             printf("%lu: %f %f\n", i, actual_de, dE[i]);
    //             assert(false);
    //         }
    //         state[i] = !state[i];
    //     } 
    //     std::vector<double> temp = dE;
    //     for (uint64_t i = 0; i < N; i++) {
    //         double temp_ene = ene + dE[i];
    //         incremental_dE(dE, state, i);
    //         state[i] = !state[i];
    //         for (uint64_t j = 0; j < N; j++) {
    //             if (j == i) {
    //                 continue;
    //             }
    //             state[j] = !state[j];
    //             double ene_new = energy(state);
    //             double actual_de = ene_new - temp_ene;
    //             if(std::abs(actual_de - dE[j]) > 1e-10) {
    //                 printf("%lu: %f %f\n", j, actual_de, dE[j]);
    //                 assert(false);
    //             }
    //             state[j] = !state[j];
    //         };
    //         dE = temp;
    //         state[i] = !state[i];
    //     } 
    // #endif
        return ene;
    };

    double Problem::energy (bool* s) const {
        double ene = 0.0;
        for (uint64_t i = 0; i < N; i++) {
            int8_t ival = sign(s[i]);
            if (sparse) {
                const std::list<Edge>& iedges = formulation_sparse.at(i);
                for (auto e: iedges) {
                    ene += 0.5 * ival * sign(s[e.dest]) * e.weight;
                }
            } else {
                const std::vector<double>& iedges = formulation.at(i);
                for (uint64_t j = i+1; j < N; j++) {
                    ene += ival * sign(s[j]) * iedges.at(j);
                }
            }
        }
        return ene;
    };;
    double Problem::energy (const StateVector& s) const {
        double ene = 0.0;
        for (uint64_t i = 0; i < N; i++) {
            int8_t ival = sign(s[i]);
            if (sparse) {
                const std::list<Edge>& iedges = formulation_sparse.at(i);
                for (auto e: iedges) {
                    ene += 0.5 * ival * sign(s[e.dest]) * e.weight;
                }
            } else {
                const std::vector<double>& iedges = formulation.at(i);
                for (uint64_t j = i+1; j < N; j++) {
                    ene += ival * sign(s[j]) * iedges.at(j);
                }
            }
        }
        return ene;
    };
    double Problem::maxcut(const StateVector& s) const {
        double cut = 0.0;
        for (uint64_t i = 0; i < N; i++) {
            int8_t ival = s[i];
            if (sparse) {
                const std::list<Edge>& iedges = formulation_sparse.at(i);
                for (auto e: iedges) {
                    cut += ival != s[e.dest] ? 0.5 * e.weight: 0;
                }
            } else {
                const std::vector<double>& iedges = formulation.at(i);
                for (uint64_t j = i+1; j < N; j++) {
                    cut += ival != s[j] ? iedges[j]: 0;
                }
            }
        }
        return cut;
    }
    // Incremental Isakov step
    void Problem::incremental_dE(std::vector<double>& dE, const StateVector& s, uint32_t index) const {
        double local_field, ene_term;
        dE[index] *= -1;
        int8_t ival = sign(s[index]);
        local_field = 0;
        if (sparse) {
            const std::list<Edge>& iedges = formulation_sparse.at(index);
            for (auto e: iedges) {
                dE.at(e.dest) += 4 * ival * sign(s[e.dest]) * e.weight;
            }
        } else {
            const std::vector<double>& iedges = formulation.at(index);
            for (uint64_t j = 0; j < N; j++) {
                if (j == index) {
                    continue;
                }
                dE[j] += 4 * ival * sign(s[j]) * iedges.at(j);
            }
        }
    }

    Problem::Problem(std::string path, uint64_t seed, uint64_t _burn_in): rng_engine(seed), burn_in(_burn_in) {
        std::fstream infile;
        infile.open(path, std::fstream::in);

        assert(infile.is_open());
        infile >> N >> M; // Read nodes and edges
        uint64_t u, v;
        double w;
        // Sparse if less than 10% of possible edges exist
        sparse = (static_cast<double>(M) / (N*(N-1) / 2)) < 0.1;
        if (sparse) {
            formulation_sparse = SparseMatrix(N);
        } else {
            formulation = Matrix(N);
            for (uint64_t i = 0; i < N; i++) {
                formulation[i] = std::vector<double>(N, 0.0);
            }
        }
        
        int offset = (path.find(".gset") == std::string::npos) ? -1 : 0;
        
        std::string line;
        std::istringstream linestream;
        double sum = 0;
        double sqsum = 0;
        while(std::getline(infile, line)) {
            if (line.size() <= 1) {
                continue;
            }

            linestream = std::istringstream(line);
            linestream >> u >> v;
            u += offset;
            v += offset;
            if (linestream.eof()) {
                w = 1.0;
            } else {
                linestream >> w;
            }
            sum += w;
            sqsum += w * w; // add the value squared
            if (sparse) {
                formulation_sparse[u].push_back({v, w});
                formulation_sparse[v].push_back({u, w});
            } else {
                formulation[u][v] = w;
                formulation[v][u] = w;
            }
        }
        double variance = sqsum / M - (sum / M * sum / M); // store the variance
        characteristic_scale = std::max(sum / M, variance);
        infile.close();
    };
} // namespace pticm