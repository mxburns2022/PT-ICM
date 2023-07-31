#include "pt_icm.hh"
#include "utils.hh"
#include <algorithm>
#include <cmath>
#include <float.h>
#define DOS_MIN 0.37
#define DOS_MAX 0.6
#define DOS_REPLICAS 10
#define DOS_PTS 500
#define NUMBINS 100


void get_histogram(ReplicaChain& rc, Problem& p, std::vector<std::vector<double>>& hist, double& emin, double& estep ) {
    size_t Nt = rc.get_replicas().size();
    hist = std::vector<std::vector<double>>(Nt);
    estep = 0.0;
    const std::vector<std::vector<double>>& enelog = rc.get_energy_log();
    double emax = -DBL_MAX;
    emin = DBL_MAX;
    for (size_t i = 0; i < Nt; i++) {
        const std::vector<double>& log_i = enelog[i];
        hist[i] = std::vector<double>(NUMBINS, 0);
        for (size_t j = 0; j < log_i.size(); j++) {
            emin = std::min(log_i[j], emin);
            emax = std::max(log_i[j], emax);
        }
    }
    estep = (emax - emin) / (NUMBINS - 1);
    for (size_t i = 0; i < Nt; i++) {
        const std::vector<double>& log_i = enelog[i];
        for (size_t j = 0; j < log_i.size(); j++) {
            double ene_val = log_i[j];
            size_t index = (ene_val - emin) / estep;
            hist[i][index] += 1;
        }
    }
}
std::vector<double> find_T(ReplicaChain& rc, 
                   Problem& p, 
                   const std::vector<std::vector<double>>& hist, 
                   double emin, 
                   double estep,
                   size_t Nt_target,
                   double rlim) {
    size_t Nt = rc.get_replicas().size();
   
    std::vector<double> r (Nt - 1);
    std::vector<double> init_beta (Nt);
    double sum_prev = std::accumulate(hist[0].begin(), hist[0].end(), 0.0);
    for (size_t i = 0; i < Nt - 1; i++) {
        r[i] = 0;
        size_t sum_next = std::accumulate(hist[i+1].begin(), hist[i+1].end(), 0.0);
        double beta_j = rc.get_replicas()[rc.get_indices()[i]].get_beta();
        double beta_k = rc.get_replicas()[rc.get_indices()[i+1]].get_beta();
        init_beta[i] = beta_j;
        for (size_t j = 0; j < NUMBINS; j++) {
            double ene_j, ene_k;
            ene_j = emin + estep * (j + 0.5);
            double p_ej = hist[i][j] / sum_prev;
            #pragma unroll
            for (size_t k = 0; k < NUMBINS; k++) {
                ene_k = emin + estep * (k + 0.5);
                double p_ek = hist[i+1][k] / sum_next;
                r[i] += p_ej * p_ek * std::min(1.0, exp((beta_j - beta_k) * (ene_j - ene_k)));
            }
        }
        sum_prev = sum_next;
    }
    init_beta[Nt-1] = rc.get_replicas()[rc.get_indices()[Nt-1]].get_beta();

    std::vector<double> beta (Nt_target);
    beta[0] = rc.get_replicas()[rc.get_indices()[0]].get_beta();
    double local_r, rdiff, betadiff;
    size_t curr_ind = 0;
    for (size_t i = 1; i < Nt_target; i++) {
        local_r = 0;
        rdiff = 0;
        beta[i] = beta[i - 1];
        while (true) {
            if (curr_ind < Nt - 1) {
                rdiff = r[curr_ind];
                betadiff = init_beta[curr_ind+1] - init_beta[curr_ind];
            } else {
                rdiff = r[Nt - 2];
                betadiff = init_beta[Nt - 1] - init_beta[Nt - 2];
            }
            if (local_r + rdiff < rlim) {
                local_r += rdiff;
                curr_ind += (curr_ind < Nt - 1);
                continue;
            } else {
                double diff = (rlim - local_r) / rdiff;
                beta[i] = beta[i-1] + betadiff * diff;
                break;
            }
        }
    }

    return beta;
}

// approximate the bin width with the average change in energy via short MCMC runs
std::vector<double> get_cae(Problem& p, size_t Nt, size_t trials) {
    size_t sweeps = p.get_burn_in() + DOS_PTS;
    ReplicaChain dos_run(p, DOS_REPLICAS, DOS_MIN, DOS_MAX, sweeps, true);
    std::cout << "Running initial temperature sweeps" << std::endl;
    for (size_t i = 0; i < sweeps; i++) 
        dos_run.do_sweeps(p);
    
    Matrix hist;
    double emin, estep;
    std::cout << "Processing Histogram" << std::endl;
    get_histogram(dos_run, p, hist, emin, estep);
    std::vector<double> result = find_T(dos_run, p, hist, emin, estep, Nt, 0.5);
    return result;
}

void get_typical_scale(const Problem& p, double maxBeta, size_t trials) {
    // ReplicaChain dos_run(p, );
}