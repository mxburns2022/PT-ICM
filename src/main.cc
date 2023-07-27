#include "pt_icm.hh"
#include "temps.hh"
#include "utils.hh"
#include <boost/format.hpp>
#include <time.h>

int main(int argc, char** argv) {
    assert(argc > 1);
    size_t seed = time(NULL);
    if (argc > 2) {
        seed = std::stoul(argv[2]);
    }
    std::string problem_path = argv[1];
    size_t burn_in = 2000;
    size_t sweeps = 4000;
    Problem p(problem_path, seed, 1000);
    size_t Nt = 19;
    std::vector<double> beta = get_cae(p, Nt, 0);
    std::vector<size_t> corr_sweeps = std::vector<size_t>(Nt, 1);
    std::cout << beta << std::endl;
    ChainEnsemble tester(p, Nt, 5, beta, corr_sweeps, 3000, true);
    double final_ene;
    StateVector solution;
    tester.solve(p, final_ene, solution);
    std::cout << final_ene << std::endl;
    std::cout << p.maxcut(solution) << std::endl;
    std::vector<double> specific_heat;
    std::vector<double> u;
    specific_heat.reserve(Nt);
    const ReplicaChain& chain1 = tester.get_chain1();
    double m_avg, m2_avg, e_avg, e2_avg, evar, mvar;
    std::fstream csvfile;
    csvfile.open("results.csv", std::fstream::out);
    csvfile << "T,C,X,u,m\n";
    const std::vector<std::vector<double>>& ene_samples = chain1.get_energy_log();
    const std::vector<std::vector<double>>& mag_samples = chain1.get_mag_log();
    printf("%ld\n", chain1.get_swaps());
    chain1.print_A_swap();
    if(!csvfile.is_open()) {
        fprintf(stderr, "Error opening CSV file\n");
        exit(-1);
    }
    size_t j = 0;

    for (size_t ind = 0; ind < Nt; ind++) {
        double sample_size = mag_samples[ind].size();
        m_avg = 0.0;
        m2_avg = 0.0;
        e_avg = 0.0;
        e2_avg = 0.0;
        double beta = chain1.get_replicas()[chain1.get_indices()[ind]].get_beta();
        for (size_t i = 0; i < sample_size; i++) {
            m_avg += mag_samples.at(ind).at(i);
            m2_avg += std::pow(mag_samples.at(ind).at(i), 2.0);
            e_avg += ene_samples.at(ind).at(i);
            e2_avg += ene_samples.at(ind).at(i) * ene_samples.at(ind).at(i);
        }
        m_avg /= sample_size; // std::reduce(mag_samples.begin(), mag_samples.end(), 0.0, std::plus<>()) / sample_size;
        m2_avg /= sample_size;// std::reduce(mag_samples.begin(), mag_samples.end(), 0.0, [] (double left, double right) {return left + std::pow(right, 2);})  / sample_size;
        e_avg /= sample_size; //std::reduce(ene_samples.begin(), ene_samples.end(), 0.0, std::plus<>())  / (sample_size);
        e2_avg /= sample_size; //std::reduce(ene_samples.begin(), ene_samples.end(), 0.0, [] (double left, double right) {return left + std::pow(right, 2);})  / sample_size;
        evar = e2_avg - std::pow(e_avg, 2);
        mvar = m2_avg - std::pow(m_avg, 2);
        specific_heat.push_back(std::pow(beta, 2) / p.vcount() * evar);
        j++;
        csvfile << boost::format("%f,%f,%f,%f,%f\n") % (1 / beta) % (specific_heat.back()) % (beta * p.vcount() * mvar) % ( e_avg / p.vcount()) % (m_avg);
    }
    csvfile.close();
    return 0;
}