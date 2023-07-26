#include "pt_icm.hh"
#ifndef UTILS
#define UTILS


struct Histogram {
    double binwidth;
    double start_value;
    std::vector<size_t> hist;
};

void get_histogram(Problem& p, std::vector<double>& temperatures, size_t Nt);


#endif