#include "pt_icm.hh"
ReplicaChain::ReplicaChain(Problem& p,
                           uint32_t _Nt,
                           double _beta0,
                           double _beta1): Nt(_Nt),
                                           beta0(_beta0),
                                           beta1(_beta1) {
    double beta_step = Nt > 1 ? (beta1 - beta0) / (Nt - 1) : 0;
    double beta = beta0;
    indices = std::vector<size_t>(Nt);
    replicas = std::vector<Replica>(Nt);
    for (size_t i = 0; i < Nt; i++) {
        indices[i] = i;
        replicas[i] = Replica(p, beta);
        beta += beta_step;
    }
};
void ReplicaChain::do_sweeps (Problem& p) {
    size_t i = 0;
    for (Replica& r: replicas) {
        double start_ene = r.get_energy();
        r.sweep(p);
        i++;
    }
}
bool ReplicaChain::re_consider (Problem& p, Replica& r1, Replica& r2) {
    double boltzmann = std::exp((r1.get_beta() - r2.get_beta()) * (r1.get_energy() - r2.get_energy()));
    double rng_comp = p.rand_double();
    return rng_comp < boltzmann; 
}
void ReplicaChain::replica_exchanges (Problem& p) {
    StateVector temp;
    for (size_t i = 0; i < Nt-1; i++) {
        size_t index_low = indices[i];
        size_t index_high = indices[i+1];
        Replica& low = replicas[indices[i]];
        Replica& high = replicas[indices[i+1]];
        if (re_consider(p, low, high)) {
            double old_ene1 = p.energy(low.get_state());
            double old_ene2 = p.energy(high.get_state());
            double beta1 = low.get_beta();
            double beta2 = high.get_beta();
            low.set_beta(beta2);
            high.set_beta(beta1);
            std::swap(indices[i], indices[i+1]);
#ifndef NDEBUG
            if (std::abs(old_ene1 - replicas[indices[i+1]].get_energy()) > 1e-10) {
                std::cout << old_ene1 << " " <<  replicas[i].get_energy() << std::endl;
                std::cout << old_ene2 << " " <<  replicas[i+1].get_energy() << std::endl;
                exit(-1);
            }
#endif
        }
    }
}
const StateVector& ReplicaChain::solution(double& ene) {
    ene = replicas.back().get_energy();
    return replicas.back().get_state();
}
