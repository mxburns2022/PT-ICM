#include "pt_icm.hh"
void Replica::sweep (Problem& p) {
    for (size_t i = 0; i < p.vcount(); i++) {
        if (mh_consider(p, i)) {
            mh_accept(p, i);
        }
    }
};
void  Replica::mh_accept (Problem& p, uint32_t index) {
    energy += dE->at(index);
    p.incremental_dE(dE, state, index);
    state->at(index) *= -1;
#ifndef NDEBUG
    if (std::abs(energy - p.energy(state)) > 1e-10) {
        std::cout << energy << " " <<  p.energy(state) << std::endl;
        exit(-1);
    }
#endif
    flips += 1;
}
bool Replica::mh_consider (Problem& p, uint32_t index) {
    if (dE->at(index) <= 0) {
        return true;
    }
    double boltzmann = std::exp(-beta * dE->at(index));
    double rng_comp = p.rand_double();
    return rng_comp < boltzmann; 

}
Replica::Replica(Problem& p, double _beta, StateVector* _state, std::vector<double>* _dE) : beta(_beta), flips(0), state(_state), dE(_dE) {
    for (size_t i = 0; i < state->size(); i++) {
        if (p.rand_double() > 0.5) {
            state->at(i) = -1;
        }
    }
    energy = p.calculate_dE(dE, state);

};
void Replica::exchange(Replica& other) {
    std::swap(other.energy, energy);
    std::swap(other.dE, dE);
    std::swap(other.state, state);
}