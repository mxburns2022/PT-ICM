#ifndef CAE_EST
#define CAE_EST
#include "mpi_defs.hh"
#include "problem.hh"
namespace pticm {
    /**
     * Utilities for reconstructing constant exchange ratio temperatures for a given problem
     * Builds a histogram for the problem using a short PT run on a preset temperature range, then estimates the density of states and uses it to reconstruct a 
     *  set of temperatures satisfying r(\beta_i, \beta_{i+1}) = $\sum_{E_1, E_2}P(E_1, \beta_i)P(E_2, \beta_{i+1})P(E_1<->E_2, \beta_1<->\beta_2)$, where 
     *  P(E_1, \beta_i) = exp[-beta_i*E_1] / (\sum_E (exp[-beta_i*E_1]*P(E_1))), where \sum_E (exp[-beta_i*E_1]*P(E_1))=Z is the partition function for beta_i, and 
     *  P(E) is inferred from the histogram. 
     *  
     *  For a run consisting of temperatures [beta_1, ..., beta_N], we can construct P(E) as follows (see section 8.2 of Newman 1999 "Monte Carlo Methods in Statistical Physics"):
     *      P(E)=\sum_i(N_i(E)) / ( \sum_j(n_j * Z_j^{-1} * exp[-\beta_i E]) ), where N_i(E) is the number of observations of energy E in replica i, n_j is the number of samples taken in 
     *      replica j, and Z_j is the estimated partition function of replica j, which we iteratively solve via:
     *          Z_k=\sum_E \frac{\sum_i N_i(E)}{\sum_j n_j Z_j^{-1} exp[(\beta_k - \beta_i) E]}
     *      until a fixed point is reached. We start out with some initial estimate, then iterate until convergence
     * 
    */
    // construct a single temperature and record energy levels
    // void gather_samples(PTWorker& w, double beta, size_t num_samples);
    // Given the samples gathered in the sample stage, iteratively calculate an array of beta values such that P(\beta_i<->\beta_{i+1}) = r 
    void get_beta(PTWorker& w, 
                  double beta_hist, 
                  double beta_0, 
                  double step_size, 
                  double r, 
                  size_t num_temps, 
                  size_t num_samples);
}; // namespace pticm
#endif 