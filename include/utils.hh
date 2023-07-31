#include "pt_icm.hh"
#ifndef UTILS
#define UTILS


struct Histogram {
    double binwidth;
    double start_value;
    std::vector<size_t> hist;
};
std::vector<double> get_cae(Problem& p, size_t Nt, size_t trials);

#endif