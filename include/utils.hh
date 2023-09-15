#include "mpi_defs.hh"
#ifndef UTILS
#define UTILS

// Declarations of useful utilites (may not be in utils.cc, might have moved)
#define LINECHECK {\
    printf("%d CHECK %d\n",w.wid, __LINE__);\
    fflush(stdout);\
}
namespace pticm {
        template <typename T>
    std::ostream& operator<<(std::ostream& out, std::vector<T> vec){
        out << "[";
        if (vec.size() > 0) {
            for (uint64_t i = 0; i < vec.size() - 1; i++) {
                out << vec.at(i) << ", ";
            }
            out << vec.back();

        }
        out << "]";
        return out;
    };
    template <typename T>
    std::ostream& operator<<(std::fstream& out, std::vector<T> vec){
        out << "\"[";
        if (vec.size() > 0) {
            for (uint64_t i = 0; i < vec.size() - 1; i++) {
                out << vec.at(i) << ", ";
            }
            out << vec.back();

        }
        out << "]\"";
        return out;
    };
    template <typename T>
    std::ostream& operator<<(std::ostream& out, std::list<T> vec){
            out << "[";
        size_t size = vec.size();
        if (size > 0) {
            size_t ind = 0;
            for (T& e: vec) {
                out << e;
                if (++ind < size) {
                    out << ", ";
                }
            }
        }
        out << "]";
        return out;
    };
    template <typename T>
    std::vector<T> operator/(std::vector<T>& left, std::vector<T>& right){
        assert(left.size() == right.size());
        std::vector<T> result = std::vector<T>(left.size());
        for (size_t i = 0; i < left.size(); i++) {
            result[i] = left[i] / right[i];
        }
        return result;
    };
    template <typename S, typename T>
    std::vector<S> operator/(std::vector<T>& left, std::vector<T>& right){
        assert(left.size() == right.size());
        std::vector<T> result = std::vector<T>(left.size());
        for (size_t i = 0; i < left.size(); i++) {
            result[i] = static_cast<S>(left[i]) / static_cast<S>(right[i]);
        }
        return result;
    };
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
    void get_beta(PTWorker& w,
                double estep, double step_size, 
                double r, size_t num_samples);
} 
#endif