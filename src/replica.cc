#include "pt_icm.hh"
#include <float.h>

void Replica::sweep (Problem& p) {
    for (uint64_t i = 0; i < p.vcount(); i++) {
        if (mh_consider(p, i)) {
            mh_accept(p, i);
        }
    }
};
void  Replica::mh_accept (Problem& p, uint32_t index) {
    energy += dE[index];
    p.incremental_dE(dE, state, index);
    state[index] = !state[index];
    if (energy < best_ene) {
        best_ene = energy;
        best_state = state;
    }
    // if (std::abs(energy - p.energy(state)) > 1e-10) {
    //     std::cout << "Incorrect Update " << energy << " " <<  p.energy(state) << std::endl;
    //     exit(0);
    // }
    flips += 1;
}
bool Replica::mh_consider (Problem& p, uint32_t index) {
    if (dE[index] <= 0) {
        return true;
    }
    double boltzmann = std::exp(-beta * dE[index]);
    double rng_comp = p.rand_double();
    return rng_comp < boltzmann; 

}
Replica::Replica(Problem& p, 
                 double _beta, 
                 const StateVector& init_state, 
                 const std::vector<double>& init_dE,
                 double _energy): beta(_beta), 
                                  flips(0), 
                                  state(init_state),
                                  dE(init_dE),
                                  energy(_energy){
    best_ene = DBL_MAX;
    best_state = init_state;
};
void Replica::exchange(Replica& other) {
    std::swap(other.energy, energy);
    std::swap(other.dE, dE);
    std::swap(other.state, state);
}