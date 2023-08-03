#include "pt_icm.hh"
#include "utils.hh"

void parse_args(int argc, 
                char** argv,
                uint64_t) {

}
std::ostream& operator<<(std::ostream& out, StateVector vec){
    out << "[";
    if (vec.size() > 0) {
        for (uint64_t i = 0; i < vec.size() - 1; i++) {
            out << static_cast<int32_t>(vec[i]) << " ";
        }
        out << static_cast<int32_t>(vec.back());
    }
    out << "]";
    return out;
};
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
};
std::ostream& operator<<(std::ostream& out, ReplicaChain& rc) {
    const std::vector<Replica> replicas = rc.get_replicas();
    if (replicas.size() == 0) {
        out << "()";
        return out;
    }
    for (uint64_t i = 0; i < replicas.size()-1; i++) {
        out << "(" << replicas[i].get_beta() << ": " << replicas[i].get_energy() << ", " << replicas[i].get_flips() << "), ";
    }
    out << "(" << replicas.back().get_beta() << ": " << replicas.back().get_energy()<< ", " << replicas.back().get_flips() << ")";
    return out;
}   


// void get_histogram(Problem& p, std::vector<double>& temperatures, uint64_t Nt, uint64_t nsweeps, uint64_t bins, double r) {
//     double beta0 = 0.05;
//     double beta1 = 2;
//     ReplicaChain rc(p, Nt, beta0, beta1);
//     double min_ene = __DBL_MAX__;
//     double max_ene = -__DBL_MAX__;
//     double* energies = new double[Nt * nsweeps];
//     const std::vector<Replica>& rvec = rc.get_replicas();
//     const std::vector<uint64_t>& indvec = rc.get_indices();
//     double ene;
//     for (uint64_t sw = 0; sw < nsweeps; sw++) {
//         rc.do_sweeps(p);
//         rc.replica_exchanges(p);
//         for (uint64_t i = 0; i < Nt; i++) {
//             ene = rvec[indvec[i]].get_energy();
//             min_ene = std::min(ene, min_ene);
//             max_ene = std::max(ene, max_ene);
//             energies[i*nsweeps + sw] = ene;
//         }
//     }
//     std::vector<double> init_beta_vals(Nt);
//     double step = (max_ene - min_ene) / (bins-1);
//     double* histogram = new double[Nt * bins]{0};
//     double beta_curr = beta0;
//     for (uint64_t i = 0; i < Nt; i++) {
//         init_beta_vals[i] = rc.get_replicas()[i].get_beta();
//         for (uint64_t j = 0; j < nsweeps; j++) {
//             uint64_t index = static_cast<uint64_t>((energies[i*nsweeps + j] - min_ene) / step);
//             histogram[i * bins + index] += 1;
//         }
//     }
//     double* dos = new double[bins]{0};
//     double partition = 0;
//     for (uint64_t i = 0; i < bins; i++) {
//         double sum = 0.0;
//         for (uint64_t j = 0; j < Nt; j++) {
//             sum += histogram[j * bins + i];
//             dos[i] += histogram[j * bins + i] * histogram[j * bins + i]  / nsweeps;
        
//         partition += sum;
//         dos[i] /= sum;
//     }

//     temperatures = std::vector<double>(Nt);
//     temperatures[0] = beta0;
//     uint64_t next_ene_ind = 1;
//     std::vector<double> partition_functions(Nt);
//     for (uint64_t i = 0; i < Nt; i++) {
//         for (uint64_t j = 0; j < bins; j++) {
//             for (uint64_t j = 0; j < bins; j++) {
//                 partition_functions[i] = 
//             }
//         }

//     }
//     // delete dos;
//     delete histogram;
//     delete energies; 

// };