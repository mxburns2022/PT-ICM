#ifndef PROBLEM
#define PROBLEM
#include <random>
#include "definitions.hh"
namespace pticm {
    class Problem {
        private:
            bool sparse; // if the problem is sparse (in which case the matrix is an edgelist)
            Matrix formulation; // use whichever one is determined at runtime
            SparseMatrix formulation_sparse;
            uint64_t N; // number of nodes
            uint64_t M; // number of edges
            uint64_t seed;
            std::mt19937_64 rng_engine; // Mersenne Twister random engine for RNG
            std::uniform_real_distribution<double> unif;
            double characteristic_scale;
            uint64_t burn_in; // number of burn-in steps before sampling
        public:
            Problem(std::string path, uint64_t seed = 0, uint64_t _burn_in = UINT64_MAX);

            Problem(double* _formulation, 
                    uint64_t _seed, 
                    uint64_t _burn_in, 
                    uint64_t _N): N(_N),
                            seed(_seed),
                            rng_engine(_seed),
                            unif(0,1),
                            burn_in(_burn_in) {
                        sparse = false;
                        formulation = std::vector<std::vector<double>>(_N);
                        for (uint64_t i = 0; i < _N; i++) {
                            double* currptr = _formulation + i * N;
                            formulation[i] = std::vector<double> (currptr, currptr + N);
                        }
                    };
            Problem(double* weights, 
                    int* src,
                    int* dest,
                    uint64_t _seed, 
                    uint64_t _burn_in, 
                    uint64_t _N,
                    uint64_t _M): N(_N),
                                M(_M),
                                seed(_seed),
                                rng_engine(_seed),
                                unif(0,1),
                                burn_in(_burn_in) {
                sparse = true;
                formulation_sparse = SparseMatrix(N);
                for (uint64_t i = 0; i < M; i++) {
                    Edge e1 = {static_cast<uint64_t>(dest[i]), weights[i]};
                    formulation_sparse[src[i]].push_back(e1);
                    Edge e2 = {static_cast<uint64_t>(src[i]), weights[i]};
                    formulation_sparse[dest[i]].push_back(e2);
                }
            };
            Problem() {};
            double energy (const StateVector& s) const;
            double energy (bool* s) const;
            double get_char_scale () const {return characteristic_scale;};
            double maxcut (const StateVector& s) const;
            uint64_t vcount () const {return N;};
            uint64_t ecount () const {return M;};
            double rand_double () {return unif(rng_engine);}
            uint64_t get_burn_in() const {return burn_in;}
            const std::list<Edge>& get_sparse_neighbors(uint64_t ind) const {return formulation_sparse[ind];}
            const std::vector<double>& get_dense_neighbors(uint64_t ind) const {return formulation[ind];}
            bool is_sparse () const {return sparse;}
            double calculate_dE (std::vector<double>& dE, const StateVector& s) const;
            const Matrix& get_formulation() const {return formulation;}
            const SparseMatrix& get_formulation_sparse() const {return formulation_sparse;}
            std::mt19937_64& get_rng() {return rng_engine;};
            // Incremental Isakov step
            void incremental_dE(std::vector<double>& dE, const StateVector& s, uint32_t index) const;
            void serialize(double* weights, int* src, int* dest);
            Problem& operator=(const Problem& other) {
                sparse = other.sparse;
                formulation = other.formulation;
                formulation_sparse = other.formulation_sparse;
                characteristic_scale = other.characteristic_scale;
                N = other.N;
                M = other.M;
                seed = other.seed;
                rng_engine = other.rng_engine;
                sparse = other.sparse;
                unif = other.unif;
                burn_in = other.burn_in;
                return *this;
            }

    };
} //namespace pticm
#endif