#include "pt_icm.hh"
#ifndef UTILS
#define UTILS


struct Histogram {
    double binwidth;
    double start_value;
    std::vector<uint64_t> hist;
};
std::vector<double> get_cae(Problem& p, uint64_t Nt, uint64_t trials);
std::ostream& operator<<(std::ostream& out, std::vector<uint8_t> vec);
#endif