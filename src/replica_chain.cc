#include "pt_icm.hh"
#include <float.h>
#include <boost/format.hpp>

ReplicaChain::ReplicaChain(Problem& p,
                           uint32_t _Nt,
                           double beta0,
                           double beta1,
                           bool _log_data,
                           bool space_inverse): Nt(_Nt),
                                            sweep_counter(0) {
    log_data = _log_data;
    correlations = std::vector<size_t>(Nt, 1);
    num_swaps = 0;
    beta_values = std::vector<double>(Nt);
    exchange_index = 0;
    double beta_step = Nt > 1 ? (beta1 - beta0) / (Nt - 1) : 0;
    double beta = beta0;
    swap_logs = std::vector<size_t>(Nt, 0);
    cons_logs = std::vector<size_t>(Nt, 0);
    indices = std::vector<size_t>(Nt);
    replicas = std::vector<Replica>(Nt);
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
        if (space_inverse) {
            beta_value = beta;
        } else {
            beta_value = 1 / beta;
        }
        replicas[i] = Replica(p, beta_value, init, dE, init_ene);
        beta_values[i] = beta_value;
        beta += beta_step;
    }
    if (space_inverse) {
        beta_value = beta;
    } else {
        beta_value = 1 / beta;
    }
    beta_values[Nt-1] = beta_value;
    indices[Nt-1] = Nt-1;
    replicas[Nt-1] = Replica(p, beta_value, init, dE, init_ene);
};

ReplicaChain::ReplicaChain(Problem& p,
                           uint32_t _Nt,
                           const std::vector<double>& temperatures,
                           const std::vector<size_t>& corr,
                           bool _log_data,
                           bool inverse): Nt(_Nt), 
                                           sweep_counter(0) {
    log_data = _log_data;
    correlations = corr;
    swap_logs = std::vector<size_t>(Nt, 0);
    cons_logs = std::vector<size_t>(Nt, 0);
    exchange_index = 0;
    num_swaps = 0;
    beta_values = std::vector<double>(Nt);
    indices = std::vector<size_t>(Nt);
    replicas = std::vector<Replica>(Nt);
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
            beta_value = temperatures[i];
        } else {
            beta_value = 1 / temperatures[i];
        }
        beta_values[i] = beta_value;
        replicas[i] = Replica(p, beta_values[i], init, dE, init_ene);
    }
    if (inverse) {
        beta_value = temperatures[Nt-1];
    } else {
        beta_value = 1 / temperatures[Nt-1];
    }
    beta_values[Nt-1] = beta_value;
    indices[Nt-1] = Nt-1;
    replicas[Nt-1] = Replica(p, beta_values[Nt-1], init, dE, init_ene);                        
}
void ReplicaChain::allocate_logs(size_t sweeps) {
    if (!log_data) {
        return;
    }
    for (std::vector<double>& i: energy_log) {
        i.reserve(sweeps);
    }
    for (std::vector<double>& i: magnetization_log) {
        i.reserve(sweeps);
    }
}

void ReplicaChain::do_sweeps (Problem& p) {

    if (exchange_index == Nt - 1) {
        exchange_index = 0;
        std::shuffle(swaps.begin(), swaps.end(), p.get_rng());
    }
    sweep_counter++;
    for (size_t ind : indices) {
            Replica& r = replicas[ind];
            double start_ene = r.get_energy();
            for (size_t local_moves = 0; local_moves < correlations[ind]; local_moves++)
            {
                r.sweep(p);
            }
            // std::cout << ind << " " << r.get_beta() << " " << r.get_energy() << std::endl;
    }
    size_t i = swaps[exchange_index++];
    size_t index_low = indices[i];
    size_t index_high = indices[i+1];
    Replica& low = replicas[indices[i]];
    Replica& high = replicas[indices[i+1]];
    cons_logs[i]+=1;
    if (re_consider(p, low, high)) {
        num_swaps ++;
        double beta1 = low.get_beta();
        double beta2 = high.get_beta();
        low.set_beta(beta2);
        high.set_beta(beta1);
        std::swap(indices[i], indices[i+1]);
        swap_logs[i]+=1;

#ifndef NDEBUG
        double old_ene1 = p.energy(low.get_state());
        double old_ene2 = p.energy(high.get_state());
        if (std::abs(old_ene1 - replicas[indices[i+1]].get_energy()) > 1e-10) {
            std::cout << old_ene1 << " " <<  replicas[i].get_energy() << std::endl;
            std::cout << old_ene2 << " " <<  replicas[i+1].get_energy()  << std::endl;
            exit(-1);
        }
#endif
    
    }
    if (log_data && sweep_counter > p.get_burn_in()) {
        for (size_t ind = 0; ind < Nt; ind++) {
            Replica& r = replicas[indices[ind]];
            double m = 0.0;
            for (bool b: r.get_state()) {
                m += b ? 1.0 : 0.0;
            }
            magnetization_log[ind].push_back(m / p.vcount());
            energy_log[ind].push_back(r.get_energy());
        }
        
    }
// }

    
}


bool ReplicaChain::re_consider (Problem& p, Replica& r1, Replica& r2) {
    double dE = r1.get_energy() - r2.get_energy();
    double dB = r1.get_beta() - r2.get_beta();
    // If a higher temp (lower Beta) is at a lower energy, always exchange
    // B1 > B2 and E1 > E2
    // if (dE <= 0 && dB <= 0 || dE >= 0 && dB >= 0) {
    //     return true;
    // }
    double boltzmann = std::exp(dB * dE);
    double rng_comp = p.rand_double();
    return rng_comp < boltzmann; 
}
void ReplicaChain::replica_exchanges (Problem& p) {
    StateVector temp;
    // Strategy (NN)_a from Malakis & Pap. 2013
    if (exchange_index == Nt - 1) {
        exchange_index = 0;
        std::shuffle(swaps.begin(), swaps.end(), p.get_rng());
    }
    // size_t i = swaps[exchange_index];
    for (size_t i : swaps) {

        size_t index_low = indices[i];
        size_t index_high = indices[i+1];
        Replica& low = replicas[indices[i]];
        Replica& high = replicas[indices[i+1]];
        if (re_consider(p, low, high)) {
            double beta1 = low.get_beta();
            double beta2 = high.get_beta();
            low.set_beta(beta2);
            high.set_beta(beta1);
            std::swap(indices[i], indices[i+1]);
    #ifndef NDEBUG
            double old_ene1 = p.energy(low.get_state());
            double old_ene2 = p.energy(high.get_state());
            if (std::abs(old_ene1 - replicas[indices[i+1]].get_energy()) > 1e-10) {
                std::cout << old_ene1 << " " <<  replicas[i].get_energy() << std::endl;
                std::cout << old_ene2 << " " <<  replicas[i+1].get_energy() << std::endl;
                exit(-1);
            }
    #endif
        }
    }
    exchange_index ++;
    
}
const StateVector ReplicaChain::solution(double& ene) {
    StateVector result;
    for (Replica& r : replicas) {
        if (r.get_energy() < ene) {
            ene = r.get_best_ene();
            result = r.get_best_state();
        }
    }
    return result;
}
void ReplicaChain::solve(Problem& p,
                         size_t pts,
                         double& ene,
                         StateVector& result,
                         Initialization init_type) {
    for (size_t i = 0; i < pts; i ++) {
        do_sweeps(p);
    }
    ene = DBL_MAX;
    for (Replica& r : replicas) {
        if (r.get_energy() < ene) {
            ene = r.get_best_ene();
            result = r.get_best_state();
        }
    }

}