#include "mpi_defs.hh"
#ifndef UTILS
#define UTILS

#define LINECHECK {\
    printf("CHECK %d\n", __LINE__);\
    fflush(stdout);\
}
namespace pticm {
    inline
    void write_update ( PTWorker& w, std::string message, uint32_t level = 0) {
        if (w.wid == 0 && level <= w.params.verbosity) {
            printf("%s\n", message.c_str());
            fflush(stdout);
        }
    }
    inline
    int sign(bool var) {
        return (static_cast<int>(var)*2 - 1);
    }
    int parse_args (int argc, 
                    char** argv,
                    std::string& problem_path,
                    Parameters& params);
    void beta_linspace(double* temps, 
                       double T0, 
                       double T1, 
                       uint64_t num_temp, 
                       bool inverse);
} 
#endif