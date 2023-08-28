#ifndef DEFS
#define DEFS
#include <list>
#include <vector>
#include <iostream>
#include "mpi.h"
namespace pticm {
    template <typename T>
    std::ostream& operator<<(std::ostream& out, std::vector<T> vec) {
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
    std::ostream& operator<<(std::ostream& out, std::list<T> vec) {
        out << "[";
        size_t size = vec.size();
        if (size > 0) {
            size_t ind = 0;
            for (T& e : vec) {
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
    std::vector<T> operator/(std::vector<T>& left, std::vector<T>& right) {
        assert(left.size() == right.size());
        std::vector<T> result = std::vector<T>(left.size());
        for (size_t i = 0; i < left.size(); i++) {
            result[i] = left[i] / right[i];
        }
        return result;
    };
    template <typename S, typename T>
    std::vector<S> operator/(std::vector<T>& left, std::vector<T>& right) {
        assert(left.size() == right.size());
        std::vector<T> result = std::vector<T>(left.size());
        for (size_t i = 0; i < left.size(); i++) {
            result[i] = static_cast<S>(left[i]) / static_cast<S>(right[i]);
        }
        return result;
    };
    struct Edge {
        uint64_t dest;
        double weight;
    };
    typedef std::vector<bool> StateVector; // for convenience, use boolean vectors
    std::ostream& operator<<(std::ostream& out, StateVector vec);
    std::vector<double> operator/(std::vector<uint32_t>& left, std::vector<uint32_t>& right);
    typedef std::vector<std::vector<double>> Matrix;
    typedef std::vector<std::list<Edge>> SparseMatrix;
    //the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

}
#endif