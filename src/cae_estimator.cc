#include "cae_estimator.hh"
#include "utils.hh"


/**
 * Utilities for spacing the temperatures to ensure a  
*/
namespace pticm {
    //_______USEFUL UTILITIES_______
    // Convert the energy value `e` to a histogram bin index 
    inline size_t to_bin(const double e, const double emin, const double estep) {
        return static_cast<size_t>(std::floor((e - emin) / estep));
    }
    // convert a bin index to an energy value (starting from emin). For log functions, emin is required to be positive
    double e_bin(const size_t index, const double estep, const double emin) {
        return index * estep + emin;
    }
    // calculate log(x1) + log(x2) (l2 is log(x2)) using more efficient trick (Newman and Barkema 1999 pg. 217 )
    inline
    double logplusexp (double x1, double l2) {
        double l1 = log(x1);
        if (l1 < l2) {
            std::swap(l1, l2);
        }
        double x = exp(l2 - l1);
        return l1 + std::log1p(x);
    }
    // calculate log(x1) + log(x2) (l1 is log(x1), l2 is log(x2)) using more efficient trick (Newman and Barkema 1999 pg. 217 )
    inline
    double logplusexp2 (double l1, double l2) {
        if (l1 < l2) {
            std::swap(l1, l2);
        }
        double x = exp(l2 - l1);
        return l1 + std::log1p(x);
    }







    // make a histogram from a provided energy log
    void make_histogram(const std::vector<double>& energy_log, size_t& nbins, double& emin, double& emax, double& estep, std::vector<uint32_t>& hist) {
        // get the min and max energy values
        emin = *std::min_element(energy_log.begin(), energy_log.end());
        emax = *std::max_element(energy_log.begin(), energy_log.end());
        // If bin number is unspecified, assume estep is then determine bin count
        if (nbins == 0) {
            assert(estep > 1e-10);
            nbins = static_cast<size_t>(std::ceil((emax - emin) / estep)) + 1;
        } else {
            // assume estep is unspecified
            estep = (emax - emin) / nbins;
        }
        // construct the histogram data structure (assume no more than 2^32 entries per bin)
        hist = std::vector<uint32_t> (nbins, 0);
        for (auto e: energy_log) {
            // get the energy index and incremement the bin
            size_t index = to_bin(e, emin, estep);
#ifndef NDEBUG
            assert(e >= emin);
            assert(e <= emax);
            assert(index < nbins);
#endif
            hist[index] ++;
        }
    }



    // TODO: Make this a configurable value. Determines number of parallel tempering steps before starting to log values for averaging
    #define BURN_IN 500
    // Conduct a preliminary PT run to construct energy histograms
    void gather_samples(PTWorker& w, std::vector<std::vector<unsigned>>& histograms,
                                     std::vector<double>& beta_values, 
                                     std::vector<double>& emin_values,
                                     double& max_emax,
                                     double& min_emin,
                                     double estep, 
                                     size_t num_samples) {
        // Ensure that only the lead process conducts prelim runs
        assert(w.wid == 0);
        // init locals and data structures
        uint64_t N = w.p.vcount();
        StateVector init_state (w.p.vcount());
        std::bernoulli_distribution bern(0.5);
        std::vector<double> energy_log;
        size_t Nt = histograms.size(); // number of temps
        std::vector<Replica> chain (Nt);
        std::vector<unsigned> swap_order (Nt-1);
        std::vector<unsigned> replica_inds (Nt);
        size_t ind = 0;
        std::vector<std::vector<double>> energylogs(Nt);
        // For each replica in the chain, fill the swap order, initialize the replica, and set its temperature
        for (auto& r: chain) {
            energylogs[ind].reserve(num_samples);
            replica_inds[ind] = ind;
            if (ind < Nt - 1)
                swap_order[ind] = ind;
            // random initial state
            for (size_t i = 0; i < N; i++) {
                init_state[i] = bern(w.p.get_rng());
            }
            std::vector<double> dE(w.p.vcount(), 0);
            // get the starting temp
            double ene_start = w.p.calculate_dE(dE, init_state);
            r = Replica (w.p, beta_values[ind++], init_state, dE, ene_start, false, BURN_IN);
        }
        // MAIN MCMC LOOP
        for (size_t i = 0; i < BURN_IN + num_samples; i++) {
            // sweep each replica
            for (size_t rind = 0; rind < Nt; rind++) {
                Replica& r = chain[replica_inds[rind]];
                r.sweep(w.p);
                // after burn in, log the energies
                if (i >= BURN_IN) {
                    energylogs[rind].push_back(r.get_energy());
                }
                ind++;
            }
            // randomly permute order to satisfy detailed balance
            std::shuffle(swap_order.begin(), swap_order.end(), w.p.get_rng());
            // attempt a swap for each neighboring pair of replicas
            for (auto index: swap_order) {
                Replica& r1 = chain[replica_inds[index]];
                Replica& r2 = chain[replica_inds[index+1]];
                double e1 = r1.get_energy();
                double e2 = r2.get_energy();
                double beta1 = beta_values[index];
                double beta2 = beta_values[index+1];
                double logboltz = (e1-e2) * (beta1-beta2);
                double rngnum = w.p.rand_double();
                if (log(rngnum) < logboltz) {
                    std::swap(replica_inds[index], replica_inds[index+1]);
                    r1.set_beta(beta2);
                    r2.set_beta(beta1);
                }
            }
        }
        // make a histogram for each temperature, get the minimum and maximum recorded temps for dos calculation
        max_emax = -__DBL_MAX__;
        for (size_t i = 0; i < Nt;i++) {
            size_t nbins = 0;
            double emax;
            make_histogram(energylogs[i], nbins, emin_values[i], emax, estep, histograms[i]);
            max_emax = std::max(emax, max_emax);
            min_emin = std::min(emin_values[i], min_emin);
        }
    }


    // calculate the partition function estimate from density of states, temperature, energy minimum and energy step
    double partition_func(std::vector<double>& logdos, double beta, double emin, double estep) {
        double logz = 1;
        bool first = true;
        for (size_t i = 0; i < logdos.size(); i++) {
            double e = e_bin(i, estep, emin);
            double logboltz = -beta * e;
            if (first) {
                first = false;
                logz = logboltz + logdos[i];
            }else {
                logz = logplusexp2(logz, logboltz + logdos[i]);
            }
        }
      
        return logz;
    } 
    // get the density of states estimate from a multi-histogram method
    void calculate_dos (const std::vector<std::vector<unsigned>>& histograms, // energy histograms at seperate temperatures
                       const std::vector<double>& partition_functions,  // the calculated partition functions (using iterative process)
                       const std::vector<double>& emin_values, // minimum values for each histogram
                       const std::vector<double>& beta_hist, // beta values used when constructing the histograms
                       double estep, // energy step
                       double min_emin, // absolute minimum energy
                       double max_emax, // absolute max energy
                       size_t n, // total number of entries recorded in each histogram
                       std::vector<double>& logdos // density of states vector (stored in log form)
                       ) {
        size_t nbins = static_cast<size_t>(std::ceil((max_emax - min_emin) / estep));
        logdos = std::vector<double> (nbins);
        double logn = std::log(n);
        double e = min_emin;
        // loop over energy values
        for (size_t i = 0; i < nbins; i++, e+=estep) {
            double lognum = 1;
            double logdenom = 1;
            bool firstnum = true;
            bool firstdenom = true;
            // for each histogram, calculate the term N[E]/(n*exp[-beta*e]/Z), where N[E] is the number of entries recorded at energy E
            for (size_t j = 0; j < histograms.size(); j++) {
                double logdenomterm = logn - partition_functions[j] - beta_hist[j] * e;
                if (firstdenom) {
                    firstdenom = false;
                    logdenom = logdenomterm;
                } else {
                    logdenom = logplusexp2(logdenomterm, logdenom);
                }
                double emin_j = emin_values[j];
                if (e < emin_j) {
                    continue;
                }
                double index = to_bin(e, emin_j, estep);
                if (index >= histograms[j].size()) {
                    continue;
                }
                if (firstnum && histograms[j][index]) {
                    firstnum = false;
                    lognum = log(histograms[j][index]);
                } else if (!firstnum && histograms[j][index]) {
                    lognum = logplusexp(histograms[j][index], lognum);
                }
            }
            if (firstnum) {
                logdos[i] = -__DBL_MAX__;
                continue;
            } else {
                logdos[i] = lognum - logdenom;
               // getchar();

            }
        }
        //getchar();
    }
    // iteratively calculate the histogram partition functions (final value of density of states also passed on)
    void get_partition_iterate (const std::vector<std::vector<unsigned>>& hist, 
                                  const std::vector<double>& beta_hist,
                                  const std::vector<double>& emin_values,
                                  std::vector<double>& logdos,
                                  std::vector<double>& partition_functions, double estep, size_t n, size_t nbins,
                                  double max_emax) {
        double max_diff;
        // initialize to all zeros
        for (size_t i = 0; i < hist.size(); i++) {
            partition_functions[i] = 0;
        }
        double emin = *std::min_element(emin_values.begin(), emin_values.end());
        do {
            max_diff = 0;
            // get the new DOS estimate
            calculate_dos(hist, partition_functions, emin_values, beta_hist, estep, emin, max_emax, n, logdos);
            // get new partition functions, determine convergence value
            for (size_t i = 0; i < hist.size(); i++) {
                double beta_i = beta_hist[i];
                double temp_z = partition_func(logdos, beta_i, emin, estep);
                max_diff += std::pow(
                    std::abs(partition_functions[i] - temp_z) / std::abs(temp_z), 2);
                partition_functions[i] = temp_z;
            }
        } while(max_diff > 1e-7);
    }
    // get the probability distribution using the density of states at a given temperature (log form)
    void get_pe_vec (const std::vector<double>& logdos, std::vector<double>& logpE, double emin, double beta, double logpartition, double estep) {
        double e = emin;
        for (size_t i = 0; i < logdos.size(); i++, e += estep) {
            logpE[i] = logdos[i] - e * beta - logpartition;
        }
    }
    // calculate the acceptance rate between neighboring replicas
    double acceptance_rate (const std::vector<double>& logdos, const std::vector<double>& p1, const std::vector<double>& p2,
                            double beta_1, double beta_2, double part_1, double part_2, double estep, double emin) {
        double r = 1;
        double p1i, p12, e1, e2;
        bool first = true;
        e1 = emin;
        #pragma GCC unroll 100
        for (size_t i = 0; i < logdos.size(); i++, e1 += estep) {
            p1i = p1[i];
            // assume beta_2 < beta_1, hence swap probability will always be 1 for e2 <= e1
            for (size_t j = 0; j <= i; j++, e2 += estep) {
                if (first) {
                    first = false;
                    r = p1i + p2[j];
                } else {
                    r = logplusexp2(r, p1i + p2[j]);
                }
            }
            e2 = e1 + estep;
            for (size_t j = i+1; j < logdos.size(); j++, e2 += estep) {
                p12  = (beta_1 - beta_2) * ( e1 - e2 );
                if (first) {
                    first = false;
                    r = p1i + p2[j] + p12;
                } else {
                    r = logplusexp2(r, p1i + p2[j] + p12);
                }
            }
        }
        // assume r is a sane value and exp(r) will not cause underflow/overflow
        return exp(r);
    }
    // Given the samples gathered in the sample stage, iteratively calculate an array of beta values such that P(\beta_i<->\beta_{i+1}) = r 
    void get_beta(PTWorker& w,
                double estep, double step_size, 
                double r, size_t num_samples){
        if (w.wid == 0) {
            size_t num_temps = w.params.Nt;
            assert(num_temps >= 1);
            size_t num_hist = num_temps;
            std::vector<std::vector<unsigned>> histograms(num_hist);
            std::vector<double> hist_betas(num_hist);
            std::vector<double> emin_values(num_hist);
            double beta_start;
            double beta_end;
            // initialize temperatures
            if (w.params.inverse) {
                beta_start = w.params.T0;
                beta_end = w.params.T1;
            } else {
                beta_start = 1 / w.params.T0;
                beta_end = 1 / w.params.T1;
            }
            if (beta_start < beta_end) {
                std::swap(beta_start, beta_end);
            }
            double beta_step = (beta_end - beta_start) / (hist_betas.size() - 1);
            std::vector<double> ene_samples;
            std::vector<unsigned> hist;
            double beta = beta_start;
            double min_emin = __DBL_MAX__;
            double max_emax = -__DBL_MAX__;
            // get the histogram samples for each temperature
            for (size_t i = 0; i < num_hist; i++, beta += beta_step)
            {
                hist_betas[i] = beta;
            }
            write_update(w, "Gathering initial samples...", 1);
            // preliminary PT runs to get data
            gather_samples(w, 
                           histograms, 
                           hist_betas, 
                           emin_values, 
                           max_emax, 
                           min_emin, 
                           estep, 
                           num_samples);


            std::vector<double> partition_functions(num_hist);
            // assume each histogram has same number of entries (sweeps - burn in)
            size_t n = std::accumulate(histograms[0].begin(), histograms[0].end(), 0, std::plus<unsigned>());
            size_t total_bins = static_cast<double>(std::ceil((max_emax - min_emin) / estep));
            std::vector<double> logdos;

            // Iteratively calculate the partition functions (stored in log form)
            write_update(w, "Determining partition functions...", 1);
            get_partition_iterate(histograms, 
                                  hist_betas, 
                                  emin_values, 
                                  logdos,
                                  partition_functions, 
                                  estep, n, total_bins, max_emax);
            // set starting beta value and get partition function
            w.beta_values[0] = beta_start;
            double partition_last =  partition_func(logdos, 
                                                    beta_start,
                                                    min_emin, 
                                                    estep);
            double step = step_size;
            std::vector<double> prev_logpE(logdos.size());
            std::vector<double> curr_logpE(logdos.size());

            write_update(w, "Selecting CAE temperatures...", 1);
            // get the probability distribution for the initial temp
            get_pe_vec(logdos, prev_logpE, min_emin, beta_start, partition_last, estep);
            size_t i;
            // add temperatures until max temp hit or replica count hit
            for (i = 1; i < num_temps; i++) {
                double beta_prev = w.beta_values[i-1];
                double& beta_curr = w.beta_values[i];
                beta_curr = beta_prev;
                double partition_curr;
                double r_curr;
                double r_last = -1;
                bool cont = false;
                do {
                    cont = false;
                    double beta_orig = beta_curr;
                    beta_curr = 1 / (1 / beta_curr + step);
                    // get partition function and probability distribution for new temp, calculate acceptance rate
                    partition_curr =  partition_func(logdos, 
                                                     beta_curr,
                                                     min_emin, 
                                                     estep); 
                    get_pe_vec(logdos, curr_logpE, min_emin, beta_curr, partition_curr, estep);
                    r_curr = acceptance_rate(logdos, prev_logpE, curr_logpE,
                                             beta_prev, beta_curr,
                                             partition_last, partition_curr, estep, min_emin);
                    // if we overshot acceptance rate, halve the step size and try again 
                    if ((r - r_curr) > 0.05) {
                        beta_curr = beta_orig;
                        step /= 2;
                        cont = true;
                        continue;
                    }
                    // if the acceptance rate is large, double the step size for the next try
                    if (r_last > 0 && std::abs(r_curr - r_last) < 5e-2) {
                        step *= 2;
                        continue;
                    }
                    if (w.params.verbosity >= 2) {
                        printf("%f->%f, %f, %f\n", beta_prev, beta_curr, r_curr, step);
                    }
                    r_last = r_curr;
                } while((r_curr > r && beta_curr > beta_end) || cont);
                if (beta_curr <= beta_end) {
                    break;
                }
                partition_last = partition_curr;
                prev_logpE = curr_logpE;
            }
        w.params.Nt = i;
        std::cout << i << " replica temperatures selected" << std::endl;
        }
    }
    
} // namespace pticm
#ifndef WHOLE
int main(int argc, char** argv) {
    pticm::PTWorker w;
    w.wid = 0;
    w.p = pticm::Problem(argv[1]);
    double step_size = -0.1;
    std::vector<double> log;
    size_t nbins = 30;
    double estep;
    estep = 1;
    std::vector<uint32_t> hist(nbins);
    w.params.T0 = 0.1;
    w.params.inverse = true;
    w.params.T1 = 10;
    w.params.Nt = 30;
    // pticm::make_histogram(log, nbins, emin, emax, estep, hist);
    pticm::get_beta(w, estep, step_size, 0.5, 500);
    // pticm::get_beta(w, beta_hist, beta_0, step_size, 0.5, num_temps, num_samples);
}
#endif