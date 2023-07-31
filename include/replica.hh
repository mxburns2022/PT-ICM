#ifndef REPLICA
#define REPLICA
#include "definitions.hh"
#include "problem.hh"
class Replica {
    private:
        StateVector state; // Current replica state in the problem state bank
        std::vector<double> dE; // deltaE pointer, will indicate location in the problem dE bank
        uint32_t flips; // Number of flips (optional counter)
        double energy; // Current configuration energy
        double beta; // Inverse temperature of this replica
        size_t offset; // 
        double best_ene;
        StateVector best_state;
        // Compare the boltzmann exponent with a log uniform random number
        bool mh_consider (Problem& p, uint32_t index);
        // Accept a proposed MH move
        void mh_accept (Problem& p, uint32_t index);
    public:
        Replica(): energy(0), beta(0), flips(0) {};
        Replica(Problem& p, 
                double _beta, 
                const StateVector& init_state, 
                const std::vector<double>& init_dE, 
                double _energy);
        inline double get_energy() const {return energy;}
        inline const StateVector& get_state() const {return state;}
        inline size_t get_flips() const {return flips;}
        inline void set_state( const StateVector& new_state) { state = new_state;}
        inline void set_dE( std::vector<double> new_dE) { dE = new_dE;}
        inline void set_ene(double new_ene) { energy = new_ene;}
        inline double get_beta() const {return beta;}
        inline void set_beta(const double new_beta) {beta = new_beta;}
        inline const StateVector& get_best_state() const { return best_state; }
        inline double get_best_ene() const { return best_ene; }

        void exchange(Replica& other);
        void log();
        Replica& operator=(const Replica& other){
            state = other.state;
            dE = other.dE;
            energy = other.energy;
            beta = other.beta;
            flips = other.flips;
            return *this;
        }
        void sweep (Problem& p);
};
#endif