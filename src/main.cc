#include "pt_icm.hh"


int main(int argc, char** argv) {
    assert(argc > 1);
    size_t seed = 0;
    if (argc > 2) {
        seed = std::stoul(argv[2]);
    }
    std::string problem_path = argv[1];
    Problem p(problem_path, seed);
    ChainEnsemble tester(p, 30, 5, 0.05, 1., 800);
    double final_ene;
    StateVector solution;
    tester.solve(p, final_ene, solution);
    std::cout << final_ene << std::endl;
    std::cout << solution << std::endl;
    return 0;
}