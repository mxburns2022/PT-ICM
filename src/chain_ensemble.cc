#include "pt_icm.hh"
#include <stack>
#include "mpi.h"
ChainEnsemble::ChainEnsemble(Problem& p,
              uint32_t _Nt,
              uint32_t _Nc,
              double _beta0,
              double _beta1,
              bool record,
              bool inverse): Nc(_Nc)  {
    log_data = record;
    correlations = std::vector<size_t>(Nt, 1);
    num_swaps = 0;
    beta_values = std::vector<double>(Nt);
    exchange_index = 0;
    double beta_step = Nt > 1 ? (_beta1 - _beta0) / (Nt - 1) : 0;
    double beta = _beta0;
    swap_logs = std::vector<size_t>(Nt, 0);
    cons_logs = std::vector<size_t>(Nt, 0);
    indices = std::vector<size_t>(Nt);
    swaps = std::vector<uint32_t>(Nt-1);
    StateVector init = StateVector(p.vcount(), false);
    energy_log = std::vector<std::vector<double>>(Nt);
    magnetization_log = std::vector<std::vector<double>>(Nt);

    std::vector<double> dE(p.vcount(), 0);
    for (size_t i = 0; i < init.size(); i++) {
        if (p.rand_double() > 0.5) {
            init[i] = true;
        }
    }
    double init_ene = p.calculate_dE(dE, init);
    double beta_value;
    for (size_t i = 0; i < Nt-1; i++) {
        swaps[i] = i;
        indices[i] = i;
        if (inverse) {
            beta_value = beta;
        } else {
            beta_value = 1 / beta;
        }
        beta_values[i] = beta_value;
        beta += beta_step;
    }
    if (inverse) {
        beta_value = beta;
    } else {
        beta_value = 1 / beta;
    }
    beta_values[Nt-1] = beta_value;
    indices[Nt-1] = Nt-1;
}
ChainEnsemble::ChainEnsemble(Problem& p,
                             uint32_t _Nt,
                             uint32_t _Nc,
                             const std::vector<double>& temps,
                             const std::vector<size_t>& corr,
                             bool record,
                             bool inverse): Nc(_Nc)  {
        chain_1 = ReplicaChain(p, _Nt, temps, corr, record, inverse);
        chain_2 = ReplicaChain(p, _Nt, temps, corr, record, inverse);       
        log_data = record;
        correlations = corr;
        swap_logs = std::vector<size_t>(Nt, 0);
        cons_logs = std::vector<size_t>(Nt, 0);
        exchange_index = 0;
        num_swaps = 0;
        beta_values = std::vector<double>(Nt);
        indices = std::vector<size_t>(Nt);
        swaps = std::vector<uint32_t>(Nt-1);
        StateVector init = StateVector(p.vcount(), false);
        std::vector<double> dE(p.vcount(), 0);
        for (size_t i = 0; i < init.size(); i++) {
            if (p.rand_double() > 0.5) {
                init[i] = true;
            }
        }
        double init_ene = p.calculate_dE(dE, init);
        double beta_value;
        energy_log = std::vector<std::vector<double>>(Nt);
        magnetization_log = std::vector<std::vector<double>>(Nt);

        for (size_t i = 0; i < Nt-1; i++) {
            swaps[i] = i;
            indices[i] = i;
            if (inverse) {
                beta_value = temps[i];
            } else {
                beta_value = 1 / temps[i];
            }
            beta_values[i] = beta_value;
        }
        if (inverse) {
            beta_value = temps[Nt-1];
        } else {
            beta_value = 1 / temps[Nt-1];
        }
        beta_values[Nt-1] = beta_value;
        indices[Nt-1] = Nt-1;                                     
        };
void ChainEnsemble::solve(Problem& p, double& ene, StateVector& result, size_t pts) {
    
    // if (c1_ene < c2_ene) {
    //     ene = c1_ene;
    //     result = chain_1_sol;
    //     return;
    // }
    // ene = c2_ene;
    // result = chain_2_sol;
    return;
}
void ChainEnsemble::bfs_icm(Problem& p, size_t seed_val, const std::vector<bool>& xorvals, Replica& r1, Replica& r2) {
    std::stack<size_t> ind_stack;
    ind_stack.push(seed_val);
    StateVector& state_1 = r1.get_state();
    StateVector& state_2 = r2.get_state();
    std::vector<bool> added(p.vcount(), false);
    size_t numflipped = 0;
    double r1_ene = r1.get_energy();
    double r2_ene = r2.get_energy();
    double start_ene = r1_ene + r2_ene;
    added[seed_val] = true;
    std::vector<double>& dE_1 = r1.get_dE();
    std::vector<double>& dE_2 = r2.get_dE();
    while(!ind_stack.empty()) {
        size_t ind = ind_stack.top();
        ind_stack.pop();
        r1_ene += dE_1[ind];
        r2_ene += dE_2[ind];
        p.incremental_dE(dE_1, state_1, ind);
        p.incremental_dE(dE_2, state_2, ind);
        state_1[ind] = !state_1[ind];
        state_2[ind] = !state_2[ind];
        numflipped ++;
        if (p.is_sparse()) {
            const std::list<Edge>& elist = p.get_sparse_neighbors(ind);
            for (Edge e: elist) {
                if (xorvals[e.dest] && !added[e.dest]) {
                    ind_stack.push(e.dest);
                    added[e.dest] = true;
                }
            }
        } else {
            const std::vector<double>& elist = p.get_dense_neighbors(ind);
            for (size_t i = 0; i < p.vcount(); i++) {
                if (std::abs(elist[i]) > 1e-10 && xorvals[i] && !added[i]) {
                    ind_stack.push(i);
                    added[i] = true;
                }
            }
        }
    }
    r1.set_ene(r1_ene);
    r2.set_ene(r2_ene);

    if(std::abs(start_ene - (r1_ene + r2_ene)) > 1e-10){
        printf("ERROR, cluster move not isoenergetic\n\t%f\t%f\t%f\n", start_ene, r1_ene + r2_ene, p.energy(state_1)+p.energy(state_2));
        exit(-1);
    }
    // std::cout << r1.get_beta() << " " << r2.get_beta() << " " << static_cast<double>(numflipped) / p.vcount() << std::endl;
}
void ChainEnsemble::perform_icm(Problem& p) {    
    std::vector<Replica>& replicas_1 = chain_1.get_replicas();
    std::vector<Replica>& replicas_2 = chain_2.get_replicas();
    const std::vector<size_t>& inds_1 = chain_1.get_indices();
    const std::vector<size_t>& inds_2 = chain_2.get_indices();
    const std::vector<double>& betas = chain_1.get_beta_values();
    size_t Nt = replicas_1.size();
    for ( size_t i = 0; i < Nt; i++) {
        std::vector<bool> xorval(p.vcount());
        std::vector<size_t> oneinds;
        oneinds.reserve(p.vcount());
        if (betas[i] < 1) {
            continue;
        }
        size_t ind1 = inds_1[i];
        size_t ind2 = inds_2[i];
        Replica& r1 = replicas_1[ind1];
        Replica& r2 = replicas_2[ind2];
        StateVector& state_1 = r1.get_state();
        StateVector& state_2 = r2.get_state();
        for (size_t j = 0; j < p.vcount(); j++) {
            size_t val = state_1[j] ^ state_2[j];
            xorval[j] = val;
            if (val == 1) {
                oneinds.push_back(j);
            }
        }
        size_t cluster_sites = oneinds.size();
        // Flip all states in R1 if cluster site counts are sufficiently large
        if (cluster_sites > p.vcount() / 2) {
            oneinds.clear();
            std::vector<double>& dE = r1.get_dE();
            for (size_t j = 0; j < p.vcount(); j++) {
                state_1[j] = !state_1[j];
                xorval[j] = !xorval[j];
                if (xorval[j]) {
                    oneinds.push_back(j);
                }
            }
        }
        std::vector<size_t> seed_val(1);
        if (oneinds.size() == 0) {
            continue;
        }
        std::sample(oneinds.begin(), oneinds.end(), seed_val.begin(), 1, p.get_rng());
        bfs_icm(p, seed_val.front(), xorval, r1, r2);
    }
};