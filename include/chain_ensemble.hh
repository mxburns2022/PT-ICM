#ifndef ENSEMBLE
#define ENSEMBLE
#include "definitions.hh"
#include "replica_chain.hh"
class Problem;
class Replica;
class ChainEnsemble {
    private:
        ReplicaChain chain_1;
        ReplicaChain chain_2;
        uint32_t Nc;
        uint32_t sweeps;
        uint32_t Nt;
        uint64_t exchange_index;
        std::vector<double> beta_values;
        uint32_t sweep_counter;
        uint64_t num_swaps;
        std::vector<uint32_t> swaps; // store the replica swap order (encoded as flattened indices in an NtxNt array)
        std::vector<uint64_t> indices; // store the replica indices in order of ascending beta
        std::vector<uint64_t> correlations; // number of sweeps for each temperature
        // Logs to store data for specific heat, entropy, susceptibility values
        std::vector<std::vector<double>> energy_log;
        std::vector<std::vector<double>> magnetization_log;
        std::vector<uint64_t> swap_logs;
        std::vector<uint64_t> cons_logs;
        bool log_data;
        void perform_icm(Problem& p);
        void bfs_icm(Problem& p, uint64_t seed, const std::vector<bool>& xorvals, Replica& r1, Replica& r2);
    public:
        ChainEnsemble(Problem& p,
                      uint32_t _Nt,
                      uint32_t _Nc,
                      double _beta0,
                      double _beta1,
                      bool record = false,
                      bool inverse = true);
        ChainEnsemble(Problem& p,
                      uint32_t _Nt,
                      uint32_t _Nc,
                      const std::vector<double>& temps,
                      const std::vector<uint64_t>& corr,
                      bool record = false,
                      bool inverse = true);
        void solve(Problem& p, double& ene, StateVector& result, uint64_t pts);
        inline const ReplicaChain& get_chain1  () const {return chain_1;};
        inline const ReplicaChain& get_chain2  () const {return chain_2;};
        uint64_t get_swaps () const {return chain_1.get_swaps() + chain_2.get_swaps();}; 
        uint64_t get_flips () const { return chain_1.get_flips() + chain_2.get_flips();};
        
};
#endif