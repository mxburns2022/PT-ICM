#include "pt_icm.hh"
#ifndef UTILS
#define UTILS


struct Histogram {
    double binwidth;
    double start_value;
    std::vector<uint64_t> hist;
};

#endif