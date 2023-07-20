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
    size_t update_iters = sweeps / 100;
    size_t update_cnt = 0;
    size_t progress_bar = 0;
    char progress[101];
    for (size_t i = 0; i < 100; i++)
        progress[i] = ' ';
    progress[100] = '\0';
    printf("%ld\% Complete [%s]", progress_bar, progress);
    std::fflush(stdout);
    for (size_t i = 0; i < sweeps; i++, update_cnt++) {
        chain_1.do_sweeps(p);
        chain_1.replica_exchanges(p);
        if (update_cnt == update_iters) {
            progress[progress_bar] = '#';
            progress_bar ++;
            printf("\r%ld%s Complete [%s]", progress_bar,"%", progress);
        std::fflush(stdout);
            update_cnt = 0;
        }
    }
    progress[99] = '#';
    printf("\r100%s Complete [%s]\n","%", progress);
    std::cout << "Done Solving" << std::endl;
    double c1_ene, c2_ene;
    const StateVector  chain_1_sol = chain_1.solution(c1_ene); 
    const StateVector chain_2_sol = chain_2.solution(c2_ene); 
    if (c1_ene < c2_ene) {
        ene = c1_ene;
        result = chain_1_sol;
        return;
    }
    ene = c2_ene;
    result = chain_2_sol;
    return;
}