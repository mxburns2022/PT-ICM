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
        void perform_icm(Problem& p);
        void bfs_icm(Problem& p, size_t seed, const std::vector<bool>& xorvals, Replica& r1, Replica& r2);
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
                      const std::vector<size_t>& corr,
                      bool record = false,
                      bool inverse = true);
        void solve(Problem& p, double& ene, StateVector& result, size_t pts);
        inline const ReplicaChain& get_chain1  () const {return chain_1;};
        inline const ReplicaChain& get_chain2  () const {return chain_2;};
        size_t get_swaps () const {return chain_1.get_swaps() + chain_2.get_swaps();}; 
        size_t get_flips () const { return chain_1.get_flips() + chain_2.get_flips();};
        
};
#endif