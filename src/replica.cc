#include "pt_icm.hh"
#include <float.h>
namespace pticm {
    // perform a sweep
    void Replica::sweep (Problem& p) {
        for (uint64_t i = 0; i < p.vcount(); i++) {
            if (mh_consider(p, i)) {
                mh_accept(p, i);
            }
            if (log_energy && sweep_counter > burn_in) {
                ene_log.push_back(energy);
            }

        }
        sweep_counter ++;
    }
    // accept the state change: increments energy, updates dE, flips spin
    void  Replica::mh_accept (Problem& p, uint32_t index) {
        energy += dE[index];
        p.incremental_dE(dE, state, index);// NOTE: This must occur before flipping the spin
        state[index] = !state[index];
        // if best seen so far, log
        if (energy < best_ene) {
            best_ene = energy;
            best_state = state;
        }
        flips += 1;
    }
    // consider a M-H proposal
    bool Replica::mh_consider (Problem& p, uint32_t index) {
        if (dE[index] <= 0) {
            return true;
        }
        double boltzmann = std::exp(-beta * dE[index]);
        double rng_comp = p.rand_double();
        return rng_comp < boltzmann; 
    }
    // Replica CTOR
    Replica::Replica(Problem& p, 
                    double _beta, 
                    const StateVector& init_state, 
                    const std::vector<double>& init_dE,
                    double _energy,
                    bool _log_ene,
                    size_t _burn_in): energy(_energy),
                                      beta(_beta), 
                                      flips(0), 
                                      state(init_state),
                                      dE(init_dE),
                                      log_energy(_log_ene),
                                      burn_in(_burn_in),
                                      sweep_counter(0) {
        best_ene = DBL_MAX;
        best_state = init_state;
    }
} // namespace pticm