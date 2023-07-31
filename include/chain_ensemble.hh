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
    public:
        ChainEnsemble(Problem& p,
                      uint32_t _Nt,
                      uint32_t _Nc,
                      double _beta0,
                      double _beta1,
                      uint32_t _sweeps,
                      bool record = false);
        ChainEnsemble(Problem& p,
                      uint32_t _Nt,
                      uint32_t _Nc,
                      const std::vector<double>& temps,
                      const std::vector<size_t>& corr,
                      uint32_t _sweeps,
                      bool record = false);
       void solve(Problem& p, double& ene, StateVector& result);
       inline const ReplicaChain& get_chain1  () const {return chain_1;};
       inline const ReplicaChain& get_chain2  () const {return chain_2;};
};
#endif