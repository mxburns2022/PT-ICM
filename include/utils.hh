#include "pt_icm.hh"
#ifndef UTILS
#define UTILS


struct Histogram {
    double binwidth;
    double start_value;
    std::vector<uint64_t> hist;
};
std::vector<double> get_cae(Problem& p, uint64_t Nt, uint64_t trials);

#endif