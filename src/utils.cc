#include "pt_icm.hh"
#include "utils.hh"

std::vector<double> operator/(std::vector<uint32_t>& left, std::vector<uint32_t>& right) {
    assert(left.size() == right.size());
    std::vector<double> result = std::vector<double>(left.size());
    for (size_t i = 0; i < left.size(); i++) {
        result[i] = right[i] > 0 ? static_cast<double>(left[i]) / right[i] : 0;
    }
    return result;
};
    std::ostream& operator<<(std::ostream& out, pticm::StateVector vec){
        if (vec.size() > 0) {
            for (uint64_t i = 0; i < vec.size() - 1; i++) {
                out << static_cast<uint32_t>(vec[i]);
            }
            out << static_cast<uint32_t>(vec.back());
        }
        return out;
    };
namespace pticm{
    // various utilites for printing to stdout
    std::ostream& operator<<(std::ostream& out, StateVector& vec){
        if (vec.size() > 0) {
            for (uint64_t i = 0; i < vec.size() - 1; i++) {
                out << static_cast<uint32_t>(vec[i]);
            }
            out << static_cast<uint32_t>(vec.back());
        }
        return out;
    }
    std::ostream& operator<<(std::ostream& out, std::list<Edge> vec){
        out << "[";
        size_t size = vec.size();
        if (size > 0) {
            size_t ind = 0;
            for (Edge& e: vec) {
                out << "(" << e.dest << ", " << e.weight << ")";
                if (++ind < size) {
                    out << ", ";
                }
            }
        }
        out << "]";
        return out;
    }
    std::ostream& operator<<(std::ostream& out, std::vector<uint8_t> vec){
        out << "[";
        if (vec.size() > 0) {
            for (uint64_t i = 0; i < vec.size() - 1; i++) {
                out << static_cast<int32_t>(vec[i]) << " ";
            }
            out << static_cast<int32_t>(vec.back());
        }
        out << "]";
        return out;
    }

    std::vector<double> operator/(std::vector<uint32_t>& left, std::vector<uint32_t>& right) {
        assert(left.size() == right.size());
        std::vector<double> result = std::vector<double>(left.size());
        for (size_t i = 0; i < left.size(); i++) {
            result[i] = right[i] > 0 ? static_cast<double>(left[i]) / right[i] : 0;
        }
        return result;
    }
} //namespace pt_icm
