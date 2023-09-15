#ifndef REPLICA
#define REPLICA
#include "definitions.hh"
#include "problem.hh"
namespace pticm {
    class Replica {
        private:
            double energy; // Current configuration energy
            double beta; // Inverse temperature of this replica
            uint32_t flips; // Number of flips (optional counter)
            StateVector state; // Current replica state in the problem state bank
            std::vector<double> dE; // deltaE pointer, will indicate location in the problem dE bank
            uint64_t offset; // 
            double best_ene;
            size_t best_ene_sweep; // store the sweep when the best ene was found
            StateVector best_state;
            // Compare the boltzmann exponent with a log uniform random number
            bool mh_consider (Problem& p, uint32_t index);
            // Accept a proposed MH move
            void mh_accept (Problem& p, uint32_t index);
            std::list<double> ene_log;
            bool log_energy;
            size_t burn_in;
            size_t sweep_counter;
        public:
            Replica(): energy(0), beta(0), flips(0) {};
            Replica(Problem& p, 
                    double _beta, 
                    const StateVector& init_state, 
                    const std::vector<double>& init_dE, 
                    double _energy,
                    bool _log_ene = false,
                    size_t _burn_in = 0);
            inline double get_energy() const {return energy;}
            inline void set_energy(double new_ene) {energy = new_ene;}
            inline const StateVector& get_state() const {return state;}
            inline StateVector& get_state() {return state;}
            inline uint64_t get_flips() const {return flips;}
            inline void set_state( const StateVector& new_state) { state = new_state;}
            inline void set_dE( std::vector<double> new_dE) { dE = new_dE;}
            inline void set_ene(double new_ene) { energy = new_ene;}
            inline double get_beta() const {return beta;}
            inline void set_beta(const double new_beta) {beta = new_beta;}
            inline std::vector<double>& get_dE() {return dE;}
            inline const StateVector& get_best_state() const { return best_state; }
            inline size_t get_best_sweep() const { return best_ene_sweep; }
            inline double get_best_ene() const { return best_ene; }
            inline const std::list<double>& get_ene_log() const { return ene_log; }

            void log();
            Replica& operator=(const Replica& other){
                state = other.state;
                dE = other.dE;
                energy = other.energy;
                best_ene = other.best_ene;
                best_ene_sweep = other.best_ene_sweep;
                ene_log = other.ene_log;
                burn_in = other.burn_in;
                log_energy = other.log_energy;
                beta = other.beta;
                flips = other.flips;
                sweep_counter = other.sweep_counter;
                return *this;
            }
            void sweep (Problem& p);
    };
} // namespace pticm
#endif