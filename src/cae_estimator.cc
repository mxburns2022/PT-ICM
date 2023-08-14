#include "cae_estimator.hh"

namespace pticm {
        inline size_t to_bin(const double e, const double emin, const double estep) {
        return static_cast<size_t>(std::floor((e - emin) / estep));
    }
    double e_bin(const size_t index, const double estep) {
        return index * estep;
    }
    double e_bin(const size_t index, const double estep, const double emin) {
        return index * estep + emin;
    }
    
    void make_histogram(const std::vector<double>& energy_log, size_t& nbins, double& emin, double& emax, double& estep, std::vector<uint32_t>& hist) {
        emin =  __DBL_MAX__;
        emax = -__DBL_MAX__;
        emin = *std::min_element(energy_log.begin(), energy_log.end());
        emax = *std::max_element(energy_log.begin(), energy_log.end());
        if (nbins == 0) {
            assert(estep > 1e-10);
            nbins = static_cast<size_t>(std::ceil((emax - emin) / estep)) + 1;
            printf("\n%f %f %ld %f\n", emax, emin, nbins, estep);fflush(stdout);
        } else {
            estep = (emax - emin) / nbins;
        }
        hist = std::vector<uint32_t> (nbins, 0);
        for (auto e: energy_log) {
            // printf("\n%f %f %f %f\n", emax, emin, e, estep);fflush(stdout);
            size_t index = to_bin(e, emin, estep);
            // std::cout << index << " " << nbins << std::endl;
            assert(e >= emin);
            assert(e <= emax);
            assert(index < nbins);
            hist[index] ++;
        }
        // for (size_t i = 0; i < nbins; i++) {
        //     printf("(%f, %d), ", e_bin(i, estep, emin), hist[i]);
        //     if (i % 8 == 0) 
        //         printf("\n");
        // }
        // printf("\n");
    }
    #define BURN_IN 200
    // construct a single temperature and record energy levels
    void gather_samples(PTWorker& w, std::vector<std::vector<unsigned>>& histograms,
                                     std::vector<double>& beta_values, 
                                     std::vector<double>& emin_values,
                                     double& max_emax,
                                     double& min_emin,
                                     double estep, 
                                     size_t num_samples) {
        assert(w.wid == 0);
        uint64_t N = w.p.vcount();
        StateVector init_state (w.p.vcount());
        std::bernoulli_distribution bern(0.5);
        std::vector<double> energy_log;
        size_t Nt = histograms.size();
        std::cout << Nt << std::endl;
        std::vector<Replica> chain (Nt);
        std::vector<unsigned> swap_order (Nt-1);
        std::vector<unsigned> replica_inds (Nt);
        size_t ind = 0;
        std::vector<std::vector<double>> energylogs(Nt);
        for (auto& r: chain) {
            energylogs[ind].reserve(num_samples);
            replica_inds[ind] = ind;
            if (ind < Nt - 1)
                swap_order[ind] = ind;
            for (size_t i = 0; i < N; i++) {
                init_state[i] = bern(w.p.get_rng());
            }
            std::vector<double> dE(w.p.vcount(), 0);
            double ene_start = w.p.calculate_dE(dE, init_state);
            r = Replica (w.p, beta_values[ind++], init_state, dE, ene_start, false, BURN_IN);
        }
        for (size_t i = 0; i < BURN_IN + num_samples; i++) {
            for (size_t rind = 0; rind < Nt; rind++) {
                Replica& r = chain[replica_inds[rind]];
                r.sweep(w.p);
                if (i >= BURN_IN) {
                    energylogs[rind].push_back(r.get_energy());
                }
                ind++;
            }

            std::shuffle(swap_order.begin(), swap_order.end(), w.p.get_rng());
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
        max_emax = -__DBL_MAX__;
        for (size_t i = 0; i < Nt;i++) {
            size_t nbins = 0;
            Replica& r = chain[replica_inds[i]];
            double emax;
            make_histogram(energylogs[i], nbins, emin_values[i], emax, estep, histograms[i]);
            std::cout << histograms[i] << std::endl;
            max_emax = std::max(emax, max_emax);
            min_emin = std::min(emin_values[i], min_emin);
        }
    };

    inline
    double logplus (double x1, double x2) {
        double l1 = log(x1);
        double l2 = log(x2);
        if (l1 < l2) {
            std::swap(l1, l2);
        }
        double x = exp(l2 - l1);
        return l1 + std::log1p(x);
    }
    inline
    double logplusexp (double x1, double l2) {
        double l1 = log(x1);
        if (l1 < l2) {
            std::swap(l1, l2);
        }
        double x = exp(l2 - l1);
        return l1 + std::log1p(x);
    }
    inline
    double logplusexp2 (double l1, double l2) {
        if (l1 < l2) {
            std::swap(l1, l2);
        }
        double x = exp(l2 - l1);
        return l1 + std::log1p(x);
    }
    // double get_partition (const std::vector<unsigned>& hist, double beta_hist, double beta_part, double estep) {
    //     double lboltzmann = (beta_part-beta_hist) * e_bin(0, estep) + log(hist[0]);
    //     double log_partition = lboltzmann;
    //     for (size_t bin = 1; bin < hist.size(); bin++) {
    //         double lboltzmann = (beta_part-beta_hist) * e_bin(bin, estep) + log(hist[bin]);
    //         log_partition = logplusexp2(log_partition, lboltzmann);
    //     }
    //     return log_partition;
    // }
    double get_partition(const std::vector<std::vector<unsigned>>& hist,
                         double emin,
                         const std::vector<double>& beta_hist,
                         const std::vector<double>& emin_values,
                         const std::vector<double>& partition_values,
                         double beta, double estep, size_t n, size_t nbins) {
        double log_partition = 1;     
        bool first = true;     
        /** The iterative formula for the partition functions is given by [1] eqn 8.31:
         * Z_k = \sum_E
         * 
         * */
        for (size_t i = 0; i < nbins; i++) {
            // get the energy of the ith bin
            double e = e_bin(i, estep, emin);
            // initialize the numerator and denominator
            double lognumerator = 1;
            double logdenominator = 1;
            bool firstnum = true;
            bool firstdenom = true;
            for (size_t j = 0; j < hist.size(); j++) {
                double logdenomterm = log(n) - partition_values[j] + (beta - beta_hist[j]) * e;
                if (firstdenom) {
                    // std::cout << j << " " << e << std::endl;
                    firstdenom = false;
                    logdenominator = logdenomterm;
                } else {
                    logdenominator = logplusexp2(logdenomterm, logdenominator);
                }
                
                double emin_j = emin_values[j];
                if (e < emin_j)
                    continue;
                size_t bins_j = hist[j].size();
                size_t index = to_bin(e, emin_j, estep);
                if (index >= bins_j) {
                    continue;
                }
                if (firstnum && hist[j][index]) {
                    firstnum = false;
                    lognumerator = log(hist[j][index]);
                } else if (!firstnum && hist[j][index]){
                    lognumerator = logplusexp(hist[j][index], lognumerator);
                }
            }
            if (firstnum)
                continue;
            if (first) {
                first = false;
                log_partition = lognumerator - logdenominator;
            } else {
                log_partition = logplusexp2(log_partition, lognumerator - logdenominator);
            }
        }
        return log_partition;
    }    

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
    void calculate_dos (const std::vector<std::vector<unsigned>>& histograms, 
                       const std::vector<double>& partition_functions, 
                       const std::vector<double>& emin_values,
                       const std::vector<double>& beta_hist,
                       double estep,
                       double min_emin,
                       double max_emax,
                       size_t n,
                       std::vector<double>& pE) {
        size_t nbins = static_cast<size_t>(std::ceil((max_emax - min_emin) / estep));
        pE = std::vector<double> (nbins);
        // loop over energy values
        double logn = std::log(n);
        double sum = 0;
        double e = min_emin;
        for (size_t i = 0; i < nbins; i++, e+=estep) {
            double lognum = 1;
            double logdenom = 1;
            bool firstnum = true;
            bool firstdenom = true;

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
                // std::cout << logdenom << " " << -partition_functions[j] << " " << -beta_hist[j] * e << " " << e << std::endl;
            }
            if (firstnum) {
                pE[i] = -__DBL_MAX__;
                continue;
            } else {
                pE[i] = lognum - logdenom;
               // getchar();

            }
        }
        //getchar();
    }
    void get_partition_iterate (const std::vector<std::vector<unsigned>>& hist, 
                                  const std::vector<double>& beta_hist,
                                  const std::vector<double>& emin_values,
                                  std::vector<double>& logdos,
                                  std::vector<double>& partition_functions, double estep, size_t n, size_t nbins,
                                  double max_emax) {
        double max_diff;
        for (size_t i = 0; i < hist.size(); i++) {
            partition_functions[i] = 0;
        }
        double emin = *std::min_element(emin_values.begin(), emin_values.end());
        do {
            max_diff = 0;
            calculate_dos(hist, partition_functions, emin_values, beta_hist, estep, emin, max_emax, n, logdos);
            for (size_t i = 0; i < hist.size(); i++) {
                double beta_i = beta_hist[i];
                double temp_z = partition_func(logdos, beta_i, emin, estep);
                max_diff += std::pow(
                    std::abs(partition_functions[i] - temp_z) / std::abs(temp_z), 2);
                partition_functions[i] = temp_z;
            }
            std::cout << partition_functions << std::endl;
            std::cout << max_diff << std::endl;
        } while(max_diff > 1e-6);
    }
   
    double get_logpe (const std::vector<unsigned>& hist, size_t ind, double beta_hist, double beta_1, double logpartition, double estep) {
        double e = e_bin(ind, estep);
        double result = exp(-(beta_1 - beta_hist) * e + log(hist[ind]) - logpartition);
        return result;
    }
    void get_pe_vec (const std::vector<double>& logdos, std::vector<double>& logpE, double emin, double beta, double logpartition, double estep) {
        size_t i = 0;
        double e = emin;
        for (size_t i = 0; i < logdos.size(); i++, e += estep) {
            logpE[i] = logdos[i] - e * beta - logpartition;
        }
    }
    double acceptance_rate (const std::vector<double>& logdos, const std::vector<double>& p1, std::vector<double>& p2,
                            double beta_1, double beta_2, double part_1, double part_2, double estep, double emin) {
        double r = 1;
        double p1i, p2i, p12, e1, e2;
        double norm1 = 0;
        double norm2 = 0;
        bool first = true;
        e1 = emin;
        
        for (size_t i = 0; i < logdos.size(); i++, e1 += estep) {
            p1i = p1[i];
            e2 = emin;
            for (size_t j = 0; j < logdos.size(); j++, e2 += estep) {
                p2[j] = logdos[j] - beta_2 * e2 - part_2;
                p12  = std::min(0.0, (beta_1 - beta_2) * ( e1 - e2 ));
                if (first) {
                    first = false;
                    r = p1i + p2[j] + p12;
                } else {
                    r = logplusexp2(r, p1i + p2[j] + p12);
                }
            }
        }
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
            double emin, emax;
            size_t nbins = 0;
            double beta = beta_start;
            double min_emin = __DBL_MAX__;
            double max_emax = -__DBL_MAX__;
            // get the histogram samples for each temperature
            for (size_t i = 0; i < num_hist; i++, beta += beta_step)
            {
                hist_betas[i] = beta;
            }
            gather_samples(w, 
                           histograms, 
                           hist_betas, 
                           emin_values, 
                           max_emax, 
                           min_emin, 
                           estep, 
                           num_samples);
            // Calculate the partition functions (stored in log form)
            std::vector<double> partition_functions(num_hist);
            size_t n = std::accumulate(histograms[0].begin(), histograms[0].end(), 0, std::plus<unsigned>());
            size_t total_bins = static_cast<double>(std::ceil((max_emax - min_emin) / estep));
            std::vector<double> logdos;
            get_partition_iterate(histograms, 
                                  hist_betas, 
                                  emin_values, 
                                  logdos,
                                  partition_functions, 
                                  estep, n, total_bins, max_emax);
            
            w.beta_values[0] = beta_start;
            double partition_last =  partition_func(logdos, 
                                                    beta_start,
                                                    min_emin, 
                                                    estep);
            // std::vector<double> uvalues = std::vector<double>(100);
            // double T0 = 1.8;
            // double T1 = 2.6;
            // double tstep = ( T1 - T0 ) / ( uvalues.size() - 1 );
            // for (size_t i = 0; i < uvalues.size(); i ++, T0+=tstep) {
            //     double partition = partition_func(logdos, 1 / T0, min_emin, estep);
            //     double e = min_emin;
            //     double u = e * exp(logdos[0] + (-(1/(T0)) * e) - partition);
            //     e += estep;
            //     for (size_t j = 1; j < logdos.size(); j++, e+=estep) {
            //         double logboltzmann = -(1/(T0)) * e;
            //         // std::cout << e << std::endl;
            //         u += e * exp(logdos[j] + logboltzmann - partition);
            //     }
            //     std::cout << u / (w.p.vcount()) << std::endl;
            // }
            std::cout << 0 <<"\t" << 1 / beta_start << " " << 0 << std::endl;
            double step = step_size;
            std::vector<double> prev_logpE(logdos.size());
            std::vector<double> curr_logpE(logdos.size());
            get_pe_vec(logdos, prev_logpE, min_emin, beta_start, partition_last, estep);
            double change_lim = 0.1;
            for (size_t i = 1; i < num_temps; i++) {
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
                    partition_curr =  partition_func(logdos, 
                                                     beta_curr,
                                                     min_emin, 
                                                     estep); 
                    r_curr = acceptance_rate(logdos, prev_logpE,curr_logpE,
                                             beta_prev, beta_curr,
                                             partition_last, partition_curr, estep, min_emin);
                    if (r > r_curr && (r - r_curr) > 0.05) {
                        beta_curr = beta_orig;
                        step /= 2;
                        cont = true;
                        continue;
                    }
                    std::cout << i <<"\t" << 1 / beta_curr << " " << r_curr << std::endl;
                    if (r_last > 0 && std::abs(r_curr - r_last) < 5e-2) {
                        step *= 2;
                        continue;
                    }
                    r_last = r_curr;
                std::cout << i <<"\t" << 1 / beta_curr << " " << r_curr << std::endl;
                } while((r_curr > r && beta_curr > beta_end) || cont);
                std::cout << i <<"\t" << 1 / beta_curr << " " << r_curr << std::endl;
                if (beta_curr <= beta_end) {
                    std::cout << beta_end << std::endl;
                    break;
                }
                partition_last = partition_curr;
                prev_logpE = curr_logpE;
            }
        }
    };
    
};// namespace pticm
#ifndef WHOLE
int main(int argc, char** argv) {
    pticm::PTWorker w;
    w.wid = 0;
    w.p = pticm::Problem(argv[1]);
    size_t num_temps = 19;
    size_t num_samples = 6000;
    double beta_0 = 1.0 / 10;
    double beta_max = 1.0 / 0.1;
    double replicas_max = 30;
    double beta_hist = 1 / 2.2;
    double step_size = -0.1;
    std::vector<double> log;
    size_t nbins = 30;
    double emin,emax,estep;
    estep = 1;
    std::vector<uint32_t> hist(nbins);
    w.params.T0 = 0.1;
    w.params.inverse = true;
    w.params.T1 = 10;
    w.params.Nt = 30;
    // pticm::make_histogram(log, nbins, emin, emax, estep, hist);
    // pticm::get_partition(hist, beta_hist, 2.0, estep);
    pticm::get_beta(w, estep, step_size, 0.5, 500);
    // std::cout << pticm::get_partition();
    // pticm::get_beta(w, beta_hist, beta_0, step_size, 0.5, num_temps, num_samples);
}
#endif