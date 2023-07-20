#include "pt_icm.hh"

double Problem::calculate_dE (std::vector<double>* dE, StateVector* s) const {
    double ene = 0.0;
    double ene_term;
    for (size_t i = 0; i < N; i++) {
        int8_t ival = s->at(i);
        if (sparse) {
            const std::list<Edge>& iedges = formulation_sparse.at(i);
            for (auto e: iedges) {
                ene_term = ival * s->at(e.dest) * e.weight;
                ene += ene_term / 2;
                dE->at(i) -= 2 * ene_term;
            }
        } else {
            const std::vector<double>& iedges = formulation.at(i);
            for (size_t j = 0; j < N; j++) {
                if (i == j) {
                    continue;
                }
                ene_term = ival * s->at(j) * iedges.at(j);
                ene += ene_term / 2;
                dE->at(i) -= 2 * ene_term;
            }
        }
    }
#ifndef NDEBUG
    StateVector state = *s;

    for (size_t i = 0; i < N; i++) {
        state[i] *= -1;
        double ene_new = energy(state);
        double actual_de = ene_new - ene;
        if(std::abs(actual_de - dE->at(i)) > 1e-10) {
            printf("%lu: %f %f\n", i, actual_de, dE[i]);
            assert(false);
        }
        state[i] *= -1;
    } 
    std::vector<double> temp = *dE;
    for (size_t i = 0; i < N; i++) {
        double temp_ene = ene + dE->at(i);
        incremental_dE(dE, &state, i);
        state[i] *= -1;
        for (size_t j = 0; j < N; j++) {
            if (j == i) {
                continue;
            }
            state[j] *= -1;
            double ene_new = energy(state);
            double actual_de = ene_new - temp_ene;
            if(std::abs(actual_de - dE->at(j)) > 1e-10) {
                printf("%lu: %f %f\n", j, actual_de, dE->at(j));
                assert(false);
            }
            state[j] *= -1;
        };
        *dE = temp;
        state[i] *= -1;

    } 
#endif
    return ene;
};
double Problem::energy (const StateVector s) const {
    double ene = 0.0;
    for (size_t i = 0; i < N; i++) {
        int8_t ival = s[i];
        if (sparse) {
            const std::list<Edge>& iedges = formulation_sparse.at(i);
            for (auto e: iedges) {
                ene += 0.5 * ival * s.at(e.dest) * e.weight;
            }
        } else {
            const std::vector<double>& iedges = formulation.at(i);
            for (size_t j = i+1; j < N; j++) {
                ene += ival * s.at(j) * iedges.at(j);
            }
        }
    }
    return ene;
};
double Problem::energy (const StateVector* s) const {
    double ene = 0.0;
    for (size_t i = 0; i < N; i++) {
        int8_t ival = s->at(i);
        if (sparse) {
            const std::list<Edge>& iedges = formulation_sparse.at(i);
            for (auto e: iedges) {
                ene += 0.5 * ival * s->at(e.dest) * e.weight;
            }
        } else {
            const std::vector<double>& iedges = formulation.at(i);
            for (size_t j = i+1; j < N; j++) {
                ene += ival * s->at(j) * iedges.at(j);
            }
        }
    }
    return ene;
};
// Incremental Isakov step
void Problem::incremental_dE(std::vector<double>* dE, StateVector* s, uint32_t index) const {
    double local_field, ene_term;
    dE->at(index) *= -1;
    int8_t ival = s->at(index);
    local_field = 0;
    if (sparse) {
        const std::list<Edge>& iedges = formulation_sparse.at(index);
        for (auto e: iedges) {
            dE->at(e.dest) += 4 * ival * s->at(e.dest) * e.weight;
        }
    } else {
        const std::vector<double>& iedges = formulation.at(index);
        for (size_t j = 0; j < N; j++) {
            if (j == index) {
                continue;
            }
            dE->at(j) += 4 * ival * s->at(j) * iedges.at(j);
        }
    }
}

Problem::Problem(std::string path, uint64_t seed ): rng_engine(seed) {
    std::fstream infile;
    infile.open(path, std::fstream::in);

    assert(infile.is_open());
    infile >> N >> M; // Read nodes and edges
    size_t u, v;
    double w;
    // Sparse if less than 10% of possible edges exist
    sparse = (static_cast<double>(M) / (N*(N-1) / 2)) < 0.1;
    if (sparse) {
        std::cout << "Sparse" << std::endl;
        formulation_sparse = SparseMatrix(N);
    } else {
        std::cout << "Dense" << std::endl;
        formulation = Matrix(N);
        for (size_t i = 0; i < N; i++) {
            formulation[i] = std::vector<double>(N, 0.0);
        }
    }

    int offset = (path.find(".gset") == std::string::npos) ? -1 : 0;
    
    std::string line;
    std::istringstream linestream;
    char_ene = 0;
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
        char_ene = std::max(std::abs(w), char_ene);
        if (sparse) {
            formulation_sparse[u].push_back({v, w});
            formulation_sparse[v].push_back({u, w});
        } else {
            formulation[u][v] = w;
            formulation[v][u] = w;
        }
    }
};