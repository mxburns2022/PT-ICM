#include "pt_icm.hh"
#include "temps.hh"
#include "utils.hh"
#include <time.h>
#include <sstream>
#include <cmath>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
int parse_args (int argc, 
                 char** argv,
                 size_t& seed,
                 std::string& problem_path,
                 size_t& Nt,
                 size_t& pts,
                 size_t& burn_in,
                 double& T0,
                 double& T1,
                 bool& inverse,
                 bool& log_data) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("seed", po::value<size_t>()->default_value(time(NULL)), "random seed")
        ("problem", po::value<std::string>(), "path to problem file (GSET)")
        ("burn_in", po::value<size_t>()->default_value(500), "number of MCMC burn-in steps before logging")
        ("replicas", po::value<size_t>()->default_value(20), "number of independent temperatures to simulate")
        ("pts", po::value<size_t>()->default_value(2000), "number of parallel tempering steps (one sweep per replica, one replica exchange attempt, one ICM)")
        ("temp0", po::value<double>()->default_value(0.1), "Starting temperature")
        ("temp1", po::value<double>()->default_value(1.6), "Ending temperature")
        ("inverse", po::value<bool>()->default_value(false), "Interpret temperatures as inverted")
        ("log", po::value<bool>()->default_value(false), "Log data and calculate observables (u and c)")
    ;
    po::variables_map opts;
    po::store(po::parse_command_line(argc, argv, desc), opts);
    if (opts.count("help")) {
        std::cout << desc << std::endl;
        return -1;
    }
    seed = opts["seed"].as<size_t>();
    problem_path = opts["problem"].as<std::string>();
    std::cout << problem_path << std::endl;
    Nt = opts["replicas"].as<size_t>();
    pts = opts["pts"].as<size_t>();
    burn_in = opts["burn_in"].as<size_t>();
    T0 = opts["temp0"].as<double>();
    T1 = opts["temp1"].as<double>();
    inverse = opts["inverse"].as<bool>();
    log_data = opts["log"].as<bool>();
    return 0;

}

void beta_linspace(double* temps, double T0, double T1, size_t num_temp, bool inverse) {
    std::cout << inverse << std::endl;

    double step = inverse ? (T1 - T0) / (num_temp - 1) : (1 / T1 - 1 / T0) / (num_temp - 1);
    std::cout << step << std::endl;
    double curr_temp = inverse ? T0 : 1 / T0;
    std::cout << curr_temp << std::endl;
    for (size_t i = 0; i < num_temp; i++) {
        if (inverse) {
            temps[i] = 1 / curr_temp;
        } else {
            temps[i] = curr_temp;
        }
        curr_temp += step;
    }
}

int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);
    int world_size, worker_id;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    // inital values
    size_t seed;
    std::string problem_path;
    size_t Nt;
    size_t pts;
    size_t burn_in;
    size_t N;
    size_t M;
    double T0;
    double T1;
    bool inverse;
    bool log_data;
    size_t help_disp = 1;
    // have worker 0 parse the arguments
    if (worker_id == 0) {
        help_disp = parse_args(argc,
                argv,
                seed,
                problem_path,
                Nt,
                pts,
                burn_in,
                T0,
                T1,
                inverse,
                log_data);
    }
    // broadcast the parsing result, exit if help message displayed
    MPI_Bcast(&help_disp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (help_disp) {
        MPI_Finalize();
        return 0;
    }
    // send the relavent parameters to each 
    MPI_Bcast((void*)&pts, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)&burn_in, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)&Nt, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)&seed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)&inverse, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)&log_data, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    
    size_t replicas_per_process = Nt * 2 / world_size;
    size_t regular = replicas_per_process * world_size;
    size_t last_process = replicas_per_process;
    if (Nt != regular) {
        last_process += (Nt * 2) - regular; 
    }
    Problem p;
    size_t local_jobs = worker_id != world_size - 1 ? replicas_per_process : last_process;
    std::unique_ptr<int[]> replica_workers;
    std::unique_ptr<double[]> beta_values;
    std::unique_ptr<double[]> beta_values_local = std::make_unique<double[]>(local_jobs);
    std::unique_ptr<double[]> weights;
    std::unique_ptr<int[]> src;
    std::unique_ptr<int[]> dest;
    std::unique_ptr<bool[]> init_state ;
    bool is_sparse;
    std::mt19937 rng;
    // Distribute the replicas to each worker
    if (worker_id == 0) {
        rng = std::mt19937(seed);
        replica_workers = std::make_unique<int[]>(Nt * 2);
        size_t ind = 0;
        for (size_t i = 0; i < world_size - 1; i++) {
            for (size_t j = 0; j < replicas_per_process; j++, ind++) {
                std::cout << ind << std::endl;
                replica_workers[ind] = i;
            }
        }
        for (size_t j = 0; j < last_process; j++, ind++) {
            replica_workers[ind] = world_size - 1;
        }
        // Determine the starting beta values for each replica
        beta_values = std::make_unique<double[]>(Nt * 2);
        beta_linspace(beta_values.get(), T0, T1, Nt, inverse);
        std::memcpy(beta_values.get() + Nt, beta_values.get(), Nt * sizeof(double));
        // copy to the beta values for worker 0
        std::memcpy(beta_values_local.get(), beta_values.get(), local_jobs * sizeof(double));
        // send local beta values to each other worker
        for (size_t i = 1; i < world_size - 1; i++) {
            double* start_pointer = i * replicas_per_process + beta_values.get();
            MPI_Send(start_pointer, replicas_per_process, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        // send to final
        double* start_pointer = (world_size - 1) * replicas_per_process + beta_values.get();
        MPI_Send(start_pointer, last_process, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD);
        
        std::bernoulli_distribution init_gen(0.5);
        init_state = std::make_unique<bool[]>(p.vcount());
        std::cout << p.vcount() << std::endl;
        for (size_t i = 0; i < p.vcount(); i++) {
            init_state[i] = init_gen(rng);
        }

        p = Problem(problem_path, seed, burn_in);

        is_sparse = p.is_sparse();
        if (is_sparse) {
            weights = std::make_unique<double[]>(p.ecount());
            src = std::make_unique<int[]>(p.ecount());
            dest = std::make_unique<int[]>(p.ecount());
        } else {
            weights = std::make_unique<double[]>(p.vcount() * p.vcount());
        }
        p.serialize(weights.get(), src.get(), dest.get());
        N = p.vcount();
        M = p.ecount();
        std::cout << N << " " << M << std::endl;
    } else {
        MPI_Recv(beta_values_local.get(), local_jobs, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, NULL);

    }
    MPI_Bcast(&is_sparse, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    if (is_sparse) {
        if (worker_id != 0) {
            weights = std::make_unique<double[]>(M);
            src = std::make_unique<int[]>(M);
            dest = std::make_unique<int[]>(M);
            std::cout << M << std::endl;
        }
        
        MPI_Bcast(weights.get(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(src.get(), M, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(dest.get(), M, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        if (worker_id != 0)
            weights = std::make_unique<double[]>(N*N);
        MPI_Bcast(weights.get(), N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    if (worker_id != 0) {

        if (is_sparse) {
            std::cout << N << std::endl;
            p = Problem(weights.get(), src.get(), dest.get(), seed, burn_in, N, M);
        } else {
            p = Problem(weights.get(), seed, burn_in, N);
        }
        init_state = std::make_unique<bool[]>(p.vcount());
    }
    std::cout << "Variable Count" << p.vcount() << std::endl;
    MPI_Bcast(init_state.get(), p.vcount(), MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    StateVector init = StateVector(init_state.get(), init_state.get()+p.vcount());
    std::vector<Replica> local_replicas = std::vector<Replica> (local_jobs);
    std::vector<double> dE(p.vcount());
    double init_ene = p.calculate_dE(dE, init); 
    size_t ind;
    for (Replica& r: local_replicas) {
        r = Replica(p, beta_values_local[ind++], init, dE, init_ene);
    }
    std::unique_ptr<double[]> energy_values;
    std::unique_ptr<size_t[]> chain_1_exchange_order;
    std::unique_ptr<size_t[]> chain_2_exchange_order;
    std::unique_ptr<size_t[]> indices;
    std::uniform_real_distribution<double> unif(0,1);
    size_t exchange_index = 0;
    if (worker_id == 0) {
        energy_values = std::make_unique<double[]> (Nt * 2);
        chain_1_exchange_order = std::make_unique<size_t[]> (Nt - 1);
        chain_2_exchange_order = std::make_unique<size_t[]> (Nt - 1);
        indices = std::make_unique<size_t[]> (Nt * 2);
        for (size_t i = 0; i < Nt - 1; i++) {
            chain_1_exchange_order[i] = i;
            chain_2_exchange_order[i] = i;
        }
        for (size_t i = 0; i < Nt * 2; i++) {
            indices[i] = i;
        }
        std::shuffle(chain_1_exchange_order.get(), chain_1_exchange_order.get() + Nt - 1, rng);
        std::shuffle(chain_2_exchange_order.get(), chain_2_exchange_order.get() + Nt - 1, rng);
    }
    std::unique_ptr<double[]> local_energy_values = std::make_unique<double[]> (local_jobs);
    size_t chain_1_exchange, chain_2_exchange;
    // Main MCMC loop
    for (size_t i = 0; i < pts; i++) {
        size_t replica_ind = 0;
        for (Replica& r: local_replicas) {
            r.sweep(p);
            local_energy_values[replica_ind++] = r.get_energy();
        }
        if (worker_id == 0) {
            chain_1_exchange = chain_1_exchange_order[exchange_index];
            chain_2_exchange = Nt + chain_2_exchange_order[exchange_index];
            for (size_t worker = 1; worker < world_size - 1; worker++) {
                MPI_Recv(energy_values.get() + replicas_per_process * worker, replicas_per_process, MPI_DOUBLE, worker, 0, MPI_COMM_WORLD, NULL);
            }
            MPI_Recv(energy_values.get() + replicas_per_process * (world_size - 1), last_process, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, NULL);
            // Perform PT comparison for chain 1
            double dBeta = beta_values[chain_1_exchange] - beta_values[chain_1_exchange + 1];
            double dE = energy_values[indices[chain_1_exchange]] - energy_values[indices[chain_1_exchange + 1]]; 
            double log_boltzmann_1 = dBeta * dE;
            double new_b1 = beta_values[chain_1_exchange+1];
            double new_e1 = energy_values[indices[chain_1_exchange+1]];
            double new_b2 = beta_values[chain_1_exchange];
            double new_e2 = energy_values[indices[chain_1_exchange]];
            if (std::log(unif(rng)) > log_boltzmann_1) {
                bool swap = true;
                double new_b1 = beta_values[chain_1_exchange+1];
                double new_b2 = beta_values[chain_1_exchange];
                MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_1_exchange], 0, MPI_COMM_WORLD);
                // send the old and new beta values along with the new energy 
                // send info to worker with the first replica
                MPI_Send(&new_b1, 1, MPI_DOUBLE, replica_workers[chain_1_exchange], 0, MPI_COMM_WORLD);
                MPI_Send(&new_b2, 1, MPI_DOUBLE, replica_workers[chain_1_exchange], 0, MPI_COMM_WORLD);
                // send info to worker with the second replica
                if (replica_workers[chain_1_exchange+1] != replica_workers[chain_1_exchange]) {
                    MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_1_exchange+1], 0, MPI_COMM_WORLD);
                    MPI_Send(&new_b2, 1, MPI_DOUBLE, replica_workers[chain_1_exchange+1], 0, MPI_COMM_WORLD);
                    MPI_Send(&new_b1, 1, MPI_DOUBLE, replica_workers[chain_1_exchange+1], 0, MPI_COMM_WORLD);
                }
                std::swap(replica_workers[chain_1_exchange], replica_workers[chain_1_exchange+1]);
                std::swap(indices[chain_1_exchange], indices[chain_1_exchange+1]);
            } else {
                bool swap = false;
                MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_1_exchange], 0, MPI_COMM_WORLD);
                MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_1_exchange+1], 0, MPI_COMM_WORLD);
            }
            for (size_t i = 1; i < world_size; i++) {
                bool dummy = false;
                if (i == replica_workers[chain_1_exchange] || i == replica_workers[chain_1_exchange + 1]) {
                    continue;
                }
                MPI_Send(&dummy, 1, MPI_CXX_BOOL, i, 0, MPI_COMM_WORLD);
            }
            dBeta = beta_values[chain_2_exchange] - beta_values[chain_2_exchange + 1];
            dE = energy_values[indices[chain_2_exchange]] - energy_values[indices[chain_2_exchange + 1]]; 
            log_boltzmann_1 = dBeta * dE;
            new_b1 = beta_values[chain_2_exchange+1];
            new_e1 = energy_values[indices[chain_2_exchange+1]];
            new_b2 = beta_values[chain_2_exchange];
            new_e2 = energy_values[indices[chain_2_exchange]];
            if (std::log(unif(rng)) > log_boltzmann_1) {
                bool swap = true;
                double new_b1 = beta_values[chain_2_exchange+1];
                double new_b2 = beta_values[chain_2_exchange];
                MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_2_exchange], 0, MPI_COMM_WORLD);
                // send the old and new beta values along with the new energy 
                // send info to worker with the first replica
                MPI_Send(&new_b1, 1, MPI_DOUBLE, replica_workers[chain_2_exchange], 0, MPI_COMM_WORLD);
                MPI_Send(&new_b2, 1, MPI_DOUBLE, replica_workers[chain_2_exchange], 0, MPI_COMM_WORLD);
                // send info to worker with the second replica
                if (replica_workers[chain_2_exchange+1] != replica_workers[chain_2_exchange]) {
                    MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_2_exchange+1], 0, MPI_COMM_WORLD);
                    MPI_Send(&new_b2, 1, MPI_DOUBLE, replica_workers[chain_2_exchange+1], 0, MPI_COMM_WORLD);
                    MPI_Send(&new_b1, 1, MPI_DOUBLE, replica_workers[chain_2_exchange+1], 0, MPI_COMM_WORLD);
                }
                std::swap(replica_workers[chain_2_exchange], replica_workers[chain_2_exchange+1]);
                std::swap(indices[chain_2_exchange], indices[chain_2_exchange+1]);
            } else {
                bool swap = false;
                MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_2_exchange], 0, MPI_COMM_WORLD);
                MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_2_exchange+1], 0, MPI_COMM_WORLD);
            }
            for (size_t i = 1; i < world_size; i++) {
                bool dummy = false;
                if (i == replica_workers[chain_2_exchange] || i == replica_workers[chain_2_exchange + 1]) {
                    continue;
                }
                MPI_Send(&dummy, 1, MPI_CXX_BOOL, i, 0, MPI_COMM_WORLD);
            }
            if (++exchange_index == Nt-1) {
                exchange_index = 0;
                std::shuffle(chain_1_exchange_order.get(), chain_1_exchange_order.get() + Nt-1, rng);
                std::shuffle(chain_2_exchange_order.get(), chain_2_exchange_order.get() + Nt-1, rng);
            }
        } else {
            MPI_Send(local_energy_values.get(), local_jobs, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            bool exchange = false;
            MPI_Recv(&exchange, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, NULL);
            if (exchange) {
                double old_beta, new_beta;
                MPI_Recv(&old_beta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, NULL);
                MPI_Recv(&new_beta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, NULL);
                for (size_t j = 0; j < local_jobs; j++) {
                    if (beta_values_local[j] == old_beta) {
                        Replica& r = local_replicas[j];
                        r.set_beta(new_beta);
                    }
                    else if (beta_values_local[j] == new_beta) {
                        Replica& r = local_replicas[j];
                        r.set_beta(old_beta);
                    }
                }
                
            }
            MPI_Recv(&exchange, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, NULL);
            if (exchange) {
                double old_beta, new_beta;
                MPI_Recv(&old_beta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, NULL);
                MPI_Recv(&new_beta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, NULL);
                for (size_t j = 0; j < local_jobs; j++) {
                    if (beta_values_local[j] == old_beta) {
                        Replica& r = local_replicas[j];
                        r.set_beta(new_beta);
                    }
                    else if (beta_values_local[j] == new_beta) {
                        Replica& r = local_replicas[j];
                        r.set_beta(old_beta);
                    }
                }
                
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
    }
    std::vector<double> local_betas (beta_values_local.get(), beta_values_local.get() + local_jobs);
    std::vector<double> local_ene (local_energy_values.get(), local_energy_values.get() + local_jobs);
    std::stringstream sstr;
    sstr << worker_id << "\n\t" << local_betas << "\n\t" << local_ene << std::endl;
    std::cout << sstr.str();

    MPI_Finalize();
    // double final_ene;
    // StateVector solution;
    // chaine.solve(p, final_ene, solution, burn_in + pts);
    // size_t last_slash = problem_path.rfind("/"); 
    // std::string graphname = problem_path.substr(last_slash+1);
    // // chaine.get_chain1().print_A_swap();
    // std::cout << "graph,energy,cut,flips,swaps,pts,Nt,nodes,seed,beta0,beta1,inverse,log_data" << std::endl <<
    //     graphname << "," <<
    //     final_ene << "," <<
    //     p.maxcut(solution) << "," <<
    //     chaine.get_flips() << "," <<
    //     chaine.get_swaps() << "," <<
    //     pts << "," <<
    //     Nt << "," <<
    //     p.vcount() << "," <<
    //     seed << "," <<
    //     beta0 << "," <<
    //     beta1 << "," <<
    //     inverse << "," <<
    //     log_data << std::endl;
    // const ReplicaChain& chain1 = chaine.get_chain1();
    // if (log_data){
    //     std::vector<double> specific_heat;
    //     std::vector<double> u;
    //     specific_heat.reserve(Nt);
    //     double m_avg, m2_avg, e_avg, e2_avg, evar, mvar;
    //     std::fstream csvfile;
    //     csvfile.open("results.csv", std::fstream::out);
    //     csvfile << "T,C,X,u,m\n";
    //     const std::vector<std::vector<double>>& ene_samples = chain1.get_energy_log();
    //     const std::vector<std::vector<double>>& mag_samples = chain1.get_mag_log();
    //     if(!csvfile.is_open()) {
    //         fprintf(stderr, "Error opening CSV file\n");
    //         exit(-1);
    //     }
    //     size_t j = 0;

    //     for (size_t ind = 0; ind < Nt; ind++) {
    //         double sample_size = mag_samples[ind].size();
    //         m_avg = 0.0;
    //         m2_avg = 0.0;
    //         e_avg = 0.0;
    //         e2_avg = 0.0;
    //         double beta = chain1.get_replicas()[chain1.get_indices()[ind]].get_beta();
    //         for (size_t i = 0; i < sample_size; i++) {
    //             m_avg += mag_samples.at(ind).at(i);
    //             m2_avg += std::pow(mag_samples.at(ind).at(i), 2.0);
    //             e_avg += ene_samples.at(ind).at(i);
    //             e2_avg += ene_samples.at(ind).at(i) * ene_samples.at(ind).at(i);
    //         }
    //         m_avg /= sample_size; // std::reduce(mag_samples.begin(), mag_samples.end(), 0.0, std::plus<>()) / sample_size;
    //         m2_avg /= sample_size;// std::reduce(mag_samples.begin(), mag_samples.end(), 0.0, [] (double left, double right) {return left + std::pow(right, 2);})  / sample_size;
    //         e_avg /= sample_size; //std::reduce(ene_samples.begin(), ene_samples.end(), 0.0, std::plus<>())  / (sample_size);
    //         e2_avg /= sample_size; //std::reduce(ene_samples.begin(), ene_samples.end(), 0.0, [] (double left, double right) {return left + std::pow(right, 2);})  / sample_size;
    //         evar = e2_avg - std::pow(e_avg, 2);
    //         mvar = m2_avg - std::pow(m_avg, 2);
    //         specific_heat.push_back(std::pow(beta, 2) / p.vcount() * evar);
    //         j++;
    //         csvfile << boost::format("%f,%f,%f,%f,%f\n") % (1 / beta) % (specific_heat.back()) % (beta * p.vcount() * mvar) % ( e_avg / p.vcount()) % (m_avg);
    //     }
    //     csvfile.close();
    // }
    
    return 0;
}