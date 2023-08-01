#ifndef CHAIN
#define CHAIN
#include "replica.hh"
// Type of PT initialization from Roma 2009
enum class Initialization {
    TypeA,
    TypeB,
    TypeC
};
class ReplicaChain {
    protected:
        std::vector<Replica> replicas;
        uint32_t Nt;
        size_t exchange_index;
        std::vector<double> beta_values;
        uint32_t sweep_counter;
        size_t num_swaps;
        std::vector<uint32_t> swaps; // store the replica swap order (encoded as flattened indices in an NtxNt array)
        std::vector<size_t> indices; // store the replica indices in order of ascending beta
        std::vector<size_t> correlations; // number of sweeps for each temperature
        // Logs to store data for specific heat, entropy, susceptibility values
        std::vector<std::vector<double>> energy_log;
        std::vector<std::vector<double>> magnetization_log;
        std::vector<size_t> swap_logs;
        std::vector<size_t> cons_logs;
        bool log_data;
    public:
        // Ctor to uniformly space (inverse) temperatures
        ReplicaChain(Problem& p,
                     uint32_t _Nt,
                     double _beta0,
                     double _beta1,
                     bool _log_data = false,
                     bool space_inverse = true); // uniformly space inverse temps
        // Ctor to initialize temperatures from given inverse list 
        ReplicaChain(Problem& p,
                     uint32_t _Nt,
                     const std::vector<double>& temperatures,
                     const std::vector<size_t>& corr,// temperatures are inverted ("beta" values)
                     bool _log_data = false,
                     bool inverse = true);
        // empty ctor
        ReplicaChain(): Nt(0) {};
        // assignment operator
        ReplicaChain& operator=(const ReplicaChain& other){
            replicas = other.replicas;
            Nt = other.Nt;
            num_swaps = other.num_swaps;
            swap_logs = other.swap_logs;
            cons_logs = other.cons_logs;
            beta_values = other.beta_values;
            sweep_counter = other.sweep_counter;
            indices = other.indices;
            swaps = other.swaps;
            energy_log = other.energy_log;
            magnetization_log = other.magnetization_log;
            correlations = other.correlations;
            log_data = other.log_data;
            exchange_index = other.exchange_index;
            return *this;
        }
        // Accessor functions
        inline const std::vector<std::vector<double>> get_energy_log() const {return energy_log;}
        inline const std::vector<std::vector<double>> get_mag_log() const {return magnetization_log;}
        inline const std::vector<double> get_beta_values() const {return beta_values;}
        const std::vector<Replica>& get_replicas () const {return replicas;}; 
        std::vector<Replica>& get_replicas () {return replicas;}; 
        const std::vector<size_t>& get_indices () const {return indices;}; 
        size_t get_swaps () const {return num_swaps;}; 
        size_t get_flips () const { 
            return std::accumulate(replicas.begin(), replicas.end(), 0, [] (size_t val, const Replica& r1) {return val + r1.get_flips();});
        }
        void allocate_logs (size_t sweeps); // allocate log space
        void update_log (); // Add energy and magnetization data to the logs
    
        // core MCMC functions 
        void solve (Problem& p, size_t pts, double& ene, StateVector& result, Initialization init_type = Initialization::TypeA); // solve problem
        void do_sweeps (Problem& p); 
        bool re_consider (Problem& p, Replica& r1, Replica& r2);
        void replica_exchanges (Problem& p);

        // print data
        void print_A_swap() const {
            std::vector<double> output(Nt);
            printf("[");
            for (size_t i = 0; i < Nt-1; i++) {
                printf("%f, ", static_cast<double>(swap_logs[i]) / cons_logs[i]);
            }
            printf("%f]\n", static_cast<double>(swap_logs.back()) / cons_logs.back());
            std::cout << swap_logs << std::endl;
            std::cout << cons_logs << std::endl;
        }
        
        // return solution
        const StateVector solution(double& ene);
};
#endif