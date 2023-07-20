#include "pt_icm.hh"

ChainEnsemble::ChainEnsemble(Problem& p,
              uint32_t _Nt,
              uint32_t _Nc,
              double _beta0,
              double _beta1,
              uint32_t _sweeps): Nc(_Nc), 
                                sweeps(_sweeps)  {
    chain_1 = ReplicaChain(p, _Nt, _beta0, _beta1);
    chain_2 = ReplicaChain(p, _Nt, _beta0, _beta1);
}
void ChainEnsemble::solve(Problem& p, double& ene, StateVector& result) {
    std::cout << "Started Solving" << std::endl;
    for (size_t i = 0; i < sweeps; i++) {
        chain_1.do_sweeps(p);
        // std::cout << "Before Exchange: " << chain_1 << std::endl;
        chain_1.replica_exchanges(p);
        // std::cout << "After Exchange: " << chain_1 << std::endl;
    }
    std::cout << "Done Solving" << std::endl;
    double c1_ene, c2_ene;
    const StateVector*  chain_1_sol = chain_1.solution(c1_ene); 
    const StateVector* chain_2_sol = chain_2.solution(c2_ene); 
    if (c1_ene < c2_ene) {
        ene = c1_ene;
        result = *chain_1_sol;
        return;
    }
    ene = c2_ene;
    result = *chain_2_sol;
    return;
}