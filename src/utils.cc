#include "pt_icm.hh"


std::ostream& operator<<(std::ostream& out, StateVector vec){
    out << "[";
    if (vec.size() > 0) {
        for (size_t i = 0; i < vec.size() - 1; i++) {
            out << static_cast<int32_t>(vec[i]) << " ";
        }
        out << static_cast<int32_t>(vec.back());
    }
    out << "]";
    return out;
};
std::ostream& operator<<(std::ostream& out, ReplicaChain& rc) {
    const std::vector<Replica> replicas = rc.get_replicas();
    if (replicas.size() == 0) {
        out << "()";
        return out;
    }
    for (size_t i = 0; i < replicas.size()-1; i++) {
        out << "(" << replicas[i].get_beta() << ": " << replicas[i].get_energy() << ", " << replicas[i].get_flips() << "), ";
    }
    out << "(" << replicas.back().get_beta() << ": " << replicas.back().get_energy()<< ", " << replicas.back().get_flips() << ")";
    return out;
}   

