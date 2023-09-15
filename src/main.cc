#include <time.h>
#include <sstream>
#include <cmath>
// #include <boost/format.hpp>
// #include <boost/program_options.hpp>
#include "solver.hh"
#include "mpi_utils.hh"
#include "argparse/argparse.hpp"
#include <filesystem>
#include <fstream>

// namespace po = boost::program_options;
// parse command line (self explanatory)
int parse_args(int argc,
               char **argv,
               std::string &problem_path,
               std::string &problem_folder,
               std::string &csv_path,
               size_t& restarts,
               pticm::Parameters & params)
{
    params.inverse = false;
    params.log_data = false;
    params.icm = false;
    params.Nt = 20;
    params.seed = time(NULL);
    params.burn_in = 500;
    params.pts = 1000;
    params.T0 = 0.6;
    params.T1 = 2.8;
    params.verbosity = 0;
    argparse::ArgumentParser program("pt_icm", "1.0");
    program.add_argument("-s", "--seed").help("random seed").default_value((uint64_t) time(NULL)).scan<'u', uint64_t>();
    program.add_argument("-p", "--problem").help("path to problem file (GSET)").required();
    program.add_argument("-b", "--burn_in").help("number of MCMC burn-in steps before logging").default_value((uint64_t)500).scan<'u', uint64_t>();
    program.add_argument("--sweeps").help("replica sweeps per parallel tempering step").default_value((uint64_t)1).scan<'u', uint64_t>();
    program.add_argument("-r", "--replicas").help("number of independent temperatures to simulate").default_value((uint64_t)20).scan<'u', uint64_t>();
    program.add_argument("--pts").help("number of parallel tempering steps (one sweep per replica, one replica exchange attempt, one ICM)").default_value((uint64_t)2000).scan<'u', uint64_t>();
    program.add_argument("-T0", "--temp0").help("Starting temperature").default_value((double)0.6).scan<'f', double>();
    program.add_argument("-T1", "--temp1").help("Ending temperature").default_value((double)2.8).scan<'f', double>();
    program.add_argument("-i", "--inverse").help("Interpret temperatures as inverted").default_value(false).implicit_value(true);
    program.add_argument("--preprocess").help("Enable CAE temperature selection").default_value(false).implicit_value(true);
    program.add_argument("--log").help("Log data and calculate observables (u and c)").default_value(false).implicit_value(true);
    program.add_argument("--icm").help("Enable isoenergetic cluster moves (ICM)").default_value(false).implicit_value(true);
    program.add_argument("--verbose").help("Verbosity level (0 for minimal messages, 1 for periodic updates, 2 for debugging)").default_value((uint32_t)0).scan<'u', uint32_t>();

    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& err) {
        std::cout << err.what() << std::endl;
        std::cout << program;
        return false;
    }

    // Show help
    if (program["--help"] == true) {
        std::cout << program;
        return false;
    }
    
    problem_path = program.get<std::string>("--problem");
    params.seed = program.get<uint64_t>("--seed");
    params.Nt = program.get<uint64_t>("--replicas");
    params.sweeps = program.get<uint64_t>("--sweeps");
    params.pts = program.get<uint64_t>("--pts");
    params.burn_in = program.get<uint64_t>("--burn_in");
    params.T0 = program.get<double>("--temp0");
    params.T1 = program.get<double>("--temp1");
    params.icm = program.get<bool>("--icm");
    params.preprocess = program.get<bool>("--preprocess");
    params.log_data = program.get<bool>("--log");
    params.inverse = program.get<bool>("--inverse");
    params.verbosity = program.get<uint32_t>("--verbose");

    return 0;
}

void print_args(pticm::Parameters& params) {
    std::cout << "Seed: " << params.seed << std::endl;
    std::cout << "Replicas: " << params.Nt << std::endl;
    std::cout << "Sweeps: " << params.sweeps << std::endl;
    std::cout << "PT Steps: " << params.pts << std::endl;
    std::cout << "Burn-in: " << params.burn_in << std::endl;
    std::cout << "T0: " << params.T0 << std::endl;
    std::cout << "T1: " << params.T1 << std::endl;
    std::cout << "Inverse: " << params.inverse << std::endl;
    std::cout << "Preprocess: " << params.preprocess << std::endl;
    std::cout << "Log Data: " << params.log_data << std::endl;
    std::cout << "ICM: " << params.icm << std::endl;
    std::cout << "Verbosity: " << params.verbosity << std::endl;
}

int main(int argc, char** argv) {
    // startup time
    MPI_Init(NULL, NULL);
    MPI_Status stat;

    pticm::PTWorker w;
    MPI_Comm_size(MPI_COMM_WORLD, &w.numw);
    MPI_Comm_rank(MPI_COMM_WORLD, &w.wid);
    // inital values
    std::string problem_path;
    std::string problem_folder;
    std::string csv;
    size_t restarts;
    uint64_t help_disp = 1;
    uint64_t folder_solve = 0;
    // have worker 0 parse the arguments

    if (w.wid == 0) {
        help_disp = parse_args(argc,
                               argv,
                               problem_path,
                               problem_folder,
                               csv,
                               restarts,
                               w.params);
        folder_solve = (problem_folder == "");
        print_args(w.params);
    }
    std::cout << problem_path << std::endl;

    // broadcast the parsing result, exit if help message displayed
    MPI_Bcast(&restarts, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&help_disp, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&folder_solve, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (help_disp) {
        MPI_Finalize();
        return 0;
    }

    MPI_Datatype parameter_type;
    pticm::create_datatype(parameter_type);
    // send the relavent parameters to each worker
    double best_ene;
    pticm::StateVector best_state;
    size_t num_flips;
    MPI_Bcast((void *)&w.params, 1, parameter_type, 0, MPI_COMM_WORLD);
    if (folder_solve) {
        // solve the problem
        pticm::solve(w, problem_path, best_ene, num_flips, best_state);
        pticm::free_worker(w);
    } else {
        std::fstream csvfile;
        if (w.wid == 0) {
            csvfile.open(csv, std::fstream::out);
            csvfile << "problem,iter,pts,replicas,tlist,pswap,icm,preprocess,burn_in,inverse,sweeps,ene,state,cut,flips\n";
        }
        std::filesystem::directory_iterator dir_iter;
        std::vector<std::string> problem_list;
        size_t num_problems = 0;
        size_t original_seed = w.params.seed;
        if (w.wid == 0)
        {
            dir_iter = std::filesystem::directory_iterator(problem_folder);
            for (auto& i: dir_iter) {
                std::string problem_instance = i.path();
                problem_list.push_back(problem_instance);
                num_problems++;
            }
        }
        MPI_Bcast(&num_problems, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        std::string problem_instance;
        for (size_t p_index = 0; p_index < num_problems; p_index++)
        {
            w.params.seed = original_seed;
            if (w.wid == 0)
            {
                problem_instance = problem_list[p_index];
            }
            pticm::Parameters &par = w.params;
            for (size_t i = 0; i < restarts; i++) {

                pticm::solve(w, problem_instance, best_ene, num_flips, best_state);
                if (w.wid == 0) {
                    std::vector<double> beta_values(w.beta_values, w.beta_values + 2*par.Nt);
                    std::vector<double> pswap = (w.swap_successes / w.swap_attempts);
                    csvfile << problem_instance << ","
                            << i << ","
                            << par.pts << ","
                            << par.Nt << ","
                            << beta_values << ","
                            << pswap << ","
                            << par.icm << ","
                            << par.preprocess << ","
                            << par.burn_in << ","
                            << par.inverse << ","
                            << par.sweeps << ","
                            << best_ene << ","
                            << "\"" << best_state << "\"" << ","
                            << w.p.maxcut(best_state) << ","
                            << num_flips << "\n";
                }
                w.params.seed++;
                MPI_Bcast(&w.params.seed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
                pticm::free_worker(w);
            }
        }
        if (w.wid == 0)
            csvfile.close();
        
    }
    MPI_Finalize();

    //__IGNORE BELOW__ Keeping for future reference, from an older version of the code that calculated thermodynamic observables
    // // local data structures
    // std::unique_ptr<int[]> replica_workers;
    // std::unique_ptr<double[]> beta_values;
    // std::unique_ptr<double[]> beta_values_local = std::make_unique<double[]>(w.loc_r);
    // std::unique_ptr<double[]> weights;
    // std::unique_ptr<int[]> src;
    // std::unique_ptr<int[]> dest;
    // std::unique_ptr<bool[]> init_state ;
    // bool is_sparse;
    // uint64_t N,M;
    // Problem& p = w.p;
    // // Distribute the replicas to each worker
    // if (w.wid == 0) {
    //     std::cout << "Constructing original problem... ";
    //     p = Problem(problem_path, w.params.seed, w.params.burn_in);
    //     std::bernoulli_distribution init_gen(0.5);
        
    //     std::cout << "Distributing replicas\n";
    //     replica_workers = std::make_unique<int[]>(Nt * 2);
    //     uint64_t ind = 0;
    //     for (uint64_t i = 0; i < w.numw - 1; i++) {
    //         for (uint64_t j = 0; j < w.worker_r; j++, ind++) {
    //             replica_workers[ind] = i;
    //         }
    //     }
    //     for (uint64_t j = 0; j < w.last_r; j++, ind++) {
    //         replica_workers[ind] = w.numw - 1;
    //     }
    //     // Determine the starting beta values for each replica
    //     beta_values = std::make_unique<double[]>(Nt * 2);
    //     beta_linspace(beta_values.get(), w.params.T0, w.params.T1, Nt, w.params.inverse);
    //     std::memcpy(beta_values.get() + Nt, beta_values.get(), Nt * sizeof(double));
    //     // copy to the beta values for worker 0
    //     std::memcpy(beta_values_local.get(), beta_values.get(), w.loc_r * sizeof(double));
    //     // send local beta values to each other worker
    //     for (uint64_t i = 1; i < w.numw - 1; i++) {
    //         double* start_pointer = i * w.worker_r + beta_values.get();
    //         MPI_Send(start_pointer, w.worker_r, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    //     }
    //     // send to final
    //     double* start_pointer = (w.numw - 1) * w.worker_r + beta_values.get();
    //     MPI_Send(start_pointer, w.last_r, MPI_DOUBLE, w.numw - 1, 0, MPI_COMM_WORLD);
        

    //     init_state = std::make_unique<bool[]>(p.vcount());
    //     for (uint64_t i = 0; i < p.vcount(); i++) {
    //         init_state[i] = init_gen(p.get_rng());
    //     }

    //     is_sparse = p.is_sparse();
    //     if (is_sparse) {
    //         weights = std::make_unique<double[]>(p.ecount());
    //         src = std::make_unique<int[]>(p.ecount());
    //         dest = std::make_unique<int[]>(p.ecount());
    //     } else {
    //         weights = std::make_unique<double[]>(p.vcount() * p.vcount());
    //     }
    //     p.serialize(weights.get(), src.get(), dest.get());
    //     N = p.vcount();
    //     M = p.ecount();
    //     std::cout << "Done\n";
    // } else {
    //     MPI_Recv(beta_values_local.get(), w.loc_r, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
    // }
    // MPI_Bcast(&is_sparse, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&N, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&M, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    // if (is_sparse) {
    //     if (w.wid != 0) {
    //         weights = std::make_unique<double[]>(M);
    //         src = std::make_unique<int[]>(M);
    //         dest = std::make_unique<int[]>(M);
    //         std::cout << M << std::endl;
    //     }
        
    //     MPI_Bcast(weights.get(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //     MPI_Bcast(src.get(), M, MPI_INT, 0, MPI_COMM_WORLD);
    //     MPI_Bcast(dest.get(), M, MPI_INT, 0, MPI_COMM_WORLD);
    // } else {
    //     if (w.wid != 0)
    //         weights = std::make_unique<double[]>(N*N);
    //     MPI_Bcast(weights.get(), N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // }
    // if (w.wid != 0) {
    //     std::cout << "Constructing problem copies... ";
    //     if (is_sparse) {
    //         std::cout << N << std::endl;
    //         p = Problem(weights.get(), src.get(), dest.get(), w.params.seed, w.params.burn_in, N, M);
    //     } else {
    //         p = Problem(weights.get(), w.params.seed, w.params.burn_in, N);
    //     }
    //     std::cout << "Done\n";
    //     init_state = std::make_unique<bool[]>(p.vcount());
    // }
    // std::cout << "Variable Count" << p.vcount() << std::endl;
    // MPI_Bcast ( init_state.get(), p.vcount(), MPI_CXX_BOOL, 0, MPI_COMM_WORLD ) ;
    // StateVector init = StateVector(init_state.get(), init_state.get()+p.vcount());
    // std::vector<Replica> local_replicas = std::vector<Replica> (w.loc_r);
    // std::vector<double> dE(p.vcount());
    // double init_ene = p.calculate_dE(dE, init); 
    // uint64_t ind = 0;
    // for (Replica& r: local_replicas) {
    //     r = Replica(p, beta_values_local[ind++], init, dE, init_ene);
    // }
    // std::unique_ptr<double[]> energy_values;
    // std::unique_ptr<uint64_t[]> chain_1_exchange_order;
    // std::unique_ptr<uint64_t[]> chain_2_exchange_order;
    // std::unique_ptr<uint64_t[]> indices;
    // uint64_t exchange_index = 0;
    // if (w.wid == 0) {
    //     energy_values = std::make_unique<double[]> (Nt * 2);
    //     chain_1_exchange_order = std::make_unique<uint64_t[]> (Nt - 1);
    //     chain_2_exchange_order = std::make_unique<uint64_t[]> (Nt - 1);
    //     indices = std::make_unique<uint64_t[]> (Nt * 2);
    //     for (uint64_t i = 0; i < Nt - 1; i++) {
    //         chain_1_exchange_order[i] = i;
    //         chain_2_exchange_order[i] = i;
    //     }
    //     for (uint64_t i = 0; i < Nt * 2; i++) {
    //         indices[i] = i;
    //     }
    //     std::shuffle(chain_1_exchange_order.get(), chain_1_exchange_order.get() + Nt - 1, p.get_rng());
    //     std::shuffle(chain_2_exchange_order.get(), chain_2_exchange_order.get() + Nt - 1, p.get_rng());
    // }
    // std::unique_ptr<double[]> local_energy_values = std::make_unique<double[]> (w.loc_r);
    // uint64_t chain_1_exchange, chain_2_exchange;
    // // Main MCMC loop
    // std::cout << "Worker " << w.wid << " entering main loop "<< w.params.pts << "\n";
    // for (uint64_t i = 0; i < w.params.pts; i++) {
    //     uint64_t replica_ind = 0;
    //     for (Replica& r: local_replicas) {
    //         r.sweep(p);
    //         local_energy_values[replica_ind++] = r.get_energy();
    //     }
    //     if (w.wid == 0) {
    //         chain_1_exchange = chain_1_exchange_order[exchange_index];
    //         chain_2_exchange = Nt + chain_2_exchange_order[exchange_index];
    //         for (uint64_t worker = 1; worker < w.numw - 1; worker++) {
    //             MPI_Recv(energy_values.get() + w.worker_r * worker, w.worker_r, MPI_DOUBLE, worker, 0, MPI_COMM_WORLD, &stat);
    //         }
    //         MPI_Recv(energy_values.get() + w.worker_r * (w.numw - 1), w.last_r, MPI_DOUBLE, w.numw - 1, 0, MPI_COMM_WORLD, &stat);
    //         std::memcpy(energy_values.get(), local_energy_values.get(), w.loc_r * sizeof(double));
    //         std::vector<double> enevals (energy_values.get(), energy_values.get() + Nt*2);
    //         // std::cout << enevals << std::endl;
    //         // Perform PT comparison for chain 1
    //         double dBeta = beta_values[chain_1_exchange] - beta_values[chain_1_exchange + 1];
    //         double dE = energy_values[indices[chain_1_exchange]] - energy_values[indices[chain_1_exchange + 1]]; 
    //         double log_boltzmann_1 = dBeta * dE;
    //         double new_b1 = beta_values[chain_1_exchange+1];
    //         double new_e1 = energy_values[indices[chain_1_exchange+1]];
    //         double new_b2 = beta_values[chain_1_exchange];
    //         double new_e2 = energy_values[indices[chain_1_exchange]];
    //         if (std::log(p.rand_double()) > log_boltzmann_1) {
    //             std::cout << "Chain 1 replica exchange\n";
    //             bool swap = true;
    //             double new_b1 = beta_values[chain_1_exchange+1];
    //             double new_b2 = beta_values[chain_1_exchange];
    //             MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_1_exchange], 0, MPI_COMM_WORLD);
    //             // send the old and new beta values along with the new energy 
    //             // send info to worker with the first replica
    //             MPI_Send(&new_b1, 1, MPI_DOUBLE, replica_workers[chain_1_exchange], 0, MPI_COMM_WORLD);
    //             MPI_Send(&new_b2, 1, MPI_DOUBLE, replica_workers[chain_1_exchange], 0, MPI_COMM_WORLD);
    //             // send info to worker with the second replica
    //             if (replica_workers[chain_1_exchange+1] != replica_workers[chain_1_exchange]) {
    //                 MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_1_exchange+1], 0, MPI_COMM_WORLD);
    //                 MPI_Send(&new_b2, 1, MPI_DOUBLE, replica_workers[chain_1_exchange+1], 0, MPI_COMM_WORLD);
    //                 MPI_Send(&new_b1, 1, MPI_DOUBLE, replica_workers[chain_1_exchange+1], 0, MPI_COMM_WORLD);
    //             }
    //             std::swap(replica_workers[chain_1_exchange], replica_workers[chain_1_exchange+1]);
    //             std::swap(indices[chain_1_exchange], indices[chain_1_exchange+1]);
    //         } else {
    //             bool swap = false;
    //             MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_1_exchange], 0, MPI_COMM_WORLD);

    //             if (replica_workers[chain_1_exchange+1] != replica_workers[chain_1_exchange]) {
    //                 MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_1_exchange+1], 0, MPI_COMM_WORLD);
    //             }
    //         }
    //         for (uint64_t i = 1; i < w.numw; i++) {
    //             bool dummy = false;
    //             if (i == replica_workers[chain_1_exchange] || i == replica_workers[chain_1_exchange + 1]) {
    //                 continue;
    //             }
    //             MPI_Send(&dummy, 1, MPI_CXX_BOOL, i, 0, MPI_COMM_WORLD);
    //         }
    //         dBeta = beta_values[chain_2_exchange] - beta_values[chain_2_exchange + 1];
    //         dE = energy_values[indices[chain_2_exchange]] - energy_values[indices[chain_2_exchange + 1]]; 
    //         log_boltzmann_1 = dBeta * dE;
    //         new_b1 = beta_values[chain_2_exchange+1];
    //         new_e1 = energy_values[indices[chain_2_exchange+1]];
    //         new_b2 = beta_values[chain_2_exchange];
    //         new_e2 = energy_values[indices[chain_2_exchange]];
    //         if (std::log(p.rand_double()) > log_boltzmann_1) {
    //             bool swap = true;
    //             double new_b1 = beta_values[chain_2_exchange+1];
    //             double new_b2 = beta_values[chain_2_exchange];
    //             MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_2_exchange], 0, MPI_COMM_WORLD);
    //             // send the old and new beta values along with the new energy 
    //             // send info to worker with the first replica
    //             MPI_Send(&new_b1, 1, MPI_DOUBLE, replica_workers[chain_2_exchange], 0, MPI_COMM_WORLD);
    //             MPI_Send(&new_b2, 1, MPI_DOUBLE, replica_workers[chain_2_exchange], 0, MPI_COMM_WORLD);
    //             // send info to worker with the second replica
    //             if (replica_workers[chain_2_exchange+1] != replica_workers[chain_2_exchange]) {
    //                 MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_2_exchange+1], 0, MPI_COMM_WORLD);
    //                 MPI_Send(&new_b2, 1, MPI_DOUBLE, replica_workers[chain_2_exchange+1], 0, MPI_COMM_WORLD);
    //                 MPI_Send(&new_b1, 1, MPI_DOUBLE, replica_workers[chain_2_exchange+1], 0, MPI_COMM_WORLD);
    //             }
    //             std::swap(replica_workers[chain_2_exchange], replica_workers[chain_2_exchange+1]);
    //             std::swap(indices[chain_2_exchange], indices[chain_2_exchange+1]);
    //         } else {
    //             bool swap = false;
    //             MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_2_exchange], 0, MPI_COMM_WORLD);
    //             MPI_Send(&swap, 1, MPI_CXX_BOOL, replica_workers[chain_2_exchange+1], 0, MPI_COMM_WORLD);
    //         }
    //         for (uint64_t i = 1; i < w.numw; i++) {
    //             bool dummy = false;
    //             if (i == replica_workers[chain_2_exchange] || i == replica_workers[chain_2_exchange + 1]) {
    //                 continue;
    //             }
    //             MPI_Send(&dummy, 1, MPI_CXX_BOOL, i, 0, MPI_COMM_WORLD);
    //         }
    //         if (++exchange_index == Nt-1) {
    //             exchange_index = 0;
    //             std::shuffle(chain_1_exchange_order.get(), chain_1_exchange_order.get() + Nt-1, p.get_rng());
    //             std::shuffle(chain_2_exchange_order.get(), chain_2_exchange_order.get() + Nt-1, p.get_rng());
    //         }
    //     } else {
    //         MPI_Send(local_energy_values.get(), w.loc_r, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    //         bool exchange = false;
    //         MPI_Recv(&exchange, 1, MPI_CXX_BOOL, 0, 0, MPI_COMM_WORLD, &stat);
    //         if (exchange) {
    //             double old_beta, new_beta;
    //             MPI_Recv(&old_beta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
    //             MPI_Recv(&new_beta, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);
    //             for (uint64_t j = 0; j < w.loc_r; j++) {
    //                 if (beta_values_local[j] == old_beta) {
    //                     Replica& r = local_replicas[j];
    //                     r.set_beta(new_beta);
    //                 }
    //                 else if (beta_values_local[j] == new_beta) {
    //                     Replica& r = local_replicas[j];
    //                     r.set_beta(old_beta);
    //                 }
    //             }
                
    //         }
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
        
    // }
    // std::cout << "Worker " << w.wid << " exiting main loop "<< w.params.pts << "\n";
    // std::vector<double> local_betas (beta_values_local.get(), beta_values_local.get() + w.loc_r);
    // std::vector<double> local_ene (local_energy_values.get(), local_energy_values.get() + w.loc_r);
    // std::stringstream sstr;
    // sstr << w.wid << "\n\t" << local_betas << "\n\t" << local_ene << std::endl;
    // std::cout << sstr.str();
    // MPI_Finalize();
    // double final_ene;
    // StateVector solution;
    // chaine.solve(p, final_ene, solution, burn_in + pts);
    // uint64_t last_slash = problem_path.rfind("/"); 
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
    //     uint64_t j = 0;

    //     for (uint64_t ind = 0; ind < Nt; ind++) {
    //         double sample_size = mag_samples[ind].size();
    //         m_avg = 0.0;
    //         m2_avg = 0.0;
    //         e_avg = 0.0;
    //         e2_avg = 0.0;
    //         double beta = chain1.get_replicas()[chain1.get_indices()[ind]].get_beta();
    //         for (uint64_t i = 0; i < sample_size; i++) {
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