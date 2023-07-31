#include "pt_icm.hh"
#include "temps.hh"
#include "utils.hh"
#include <time.h>
#include <boost/format.hpp>

void parse_args (int argc, 
                 char** argv,
                 size_t& seed,
                 std::string& problem_path,
                 size_t& Nt,
                 size_t& pts,
                 size_t& burn_in,
                 double& beta0,
                 double& beta1,
                 bool& inverse,
                 bool& log_data) {
    inverse = false;
    log_data = false;
    beta0 = 1.6;
    beta1 = 0.1;
    pts = 1000;
    burn_in = 0;
    Nt = 20;
    seed = 0;
    for (size_t i = 1; i < argc; i++) {
        std::string param_name = argv[i];
        if (param_name == "--seed" || param_name == "-s") {
            assert(i != argc - 1);
            seed = std::stoull(argv[++i]);
        }
        else if (param_name == "--problem" || param_name == "-p") {
            assert(i != argc - 1);
            problem_path = argv[++i];
            
        }
        else if (param_name == "--burn_in" || param_name == "-b") {
            assert(i != argc - 1);
            burn_in = std::stoull(argv[++i]);
            
        }
        else if (param_name == "--replicas" || param_name == "-r") {
            assert(i != argc - 1);
            Nt = std::stoull(argv[++i]);
            
        }
        else if (param_name == "--pts" || param_name == "-t") {
            assert(i != argc - 1);
            pts = std::stoull(argv[++i]);
            
        }
        else if (param_name == "--beta0" || param_name == "-b0") {
            assert(i != argc - 1);
            beta0 = std::stod(argv[++i]);
            
        }
        else if (param_name == "--beta1" || param_name == "-b1") {
            assert(i != argc - 1);
            beta1 = std::stod(argv[++i]);
        }
        else if (param_name == "--inverse" || param_name == "-i") {
            inverse = true;
        }
        else if (param_name == "--log" || param_name == "-l") {
            log_data = true;
        } else {
            printf("Unrecognized option \"%s\"\n", param_name.c_str());
            exit(-1);
        }

    }
}


int main(int argc, char** argv) {
    size_t seed;
    std::string problem_path;
    size_t Nt;
    size_t pts;
    size_t burn_in;
    double beta0;
    double beta1;
    bool inverse;
    bool log_data;

    parse_args(argc,
               argv,
               seed,
               problem_path,
               Nt,
               pts,
               burn_in,
               beta0,
               beta1,
               inverse,
               log_data);
    // std::vector<double> beta = get_cae(p, Nt, 0);
    // std::vector<size_t> corr_sweeps = std::vector<size_t>(Nt, 1);
    // std::cout << beta << std::endl
    Problem p(problem_path, seed, burn_in);
    ChainEnsemble chaine(p, Nt, 5, beta0, beta1, log_data, inverse);
    double final_ene;
    StateVector solution;
    chaine.solve(p, final_ene, solution, burn_in + pts);
    size_t last_slash = problem_path.rfind("/"); 
    std::string graphname = problem_path.substr(last_slash+1);
    // chaine.get_chain1().print_A_swap();
    std::cout << "graph,energy,cut,flips,swaps,pts,Nt,nodes,seed,beta0,beta1,inverse,log_data" << std::endl <<
        graphname << "," <<
        final_ene << "," <<
        p.maxcut(solution) << "," <<
        chaine.get_flips() << "," <<
        chaine.get_swaps() << "," <<
        pts << "," <<
        Nt << "," <<
        p.vcount() << "," <<
        seed << "," <<
        beta0 << "," <<
        beta1 << "," <<
        inverse << "," <<
        log_data << std::endl;
    const ReplicaChain& chain1 = chaine.get_chain1();
    if (log_data){
        std::vector<double> specific_heat;
        std::vector<double> u;
        specific_heat.reserve(Nt);
        double m_avg, m2_avg, e_avg, e2_avg, evar, mvar;
        std::fstream csvfile;
        csvfile.open("results.csv", std::fstream::out);
        csvfile << "T,C,X,u,m\n";
        const std::vector<std::vector<double>>& ene_samples = chain1.get_energy_log();
        const std::vector<std::vector<double>>& mag_samples = chain1.get_mag_log();
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
    }
    
    return 0;
}