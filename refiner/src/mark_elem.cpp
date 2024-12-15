#include "mark_elem.hpp"

// void mark_elem_to_refine(const std::vector<double> &eta_arr,
//                          std::vector<int> &marked_elem) {
//     std::vector<int> sorted_idx(eta_arr.size());
//     std::iota(sorted_idx.begin(), sorted_idx.end(), 0);

//     std::sort(sorted_idx.begin(), sorted_idx.end(), 
//               [&eta_arr](int i1, int i2) {return eta_arr[i1] > eta_arr[i2];});

//     for (int ielem = 0; ielem < marked_elem.size(); ielem++) {
//         marked_elem[ielem] = sorted_idx[ielem];
//     }
// }

void mark_elem_to_refine(const std::vector<double> &eta_arr,
                         const std::vector<int> &load_elem,
                         std::vector<int> &marked_elem,
                         const double &threshold) {
    double eta_max = 0.0;
    for (int ielem = 0; ielem < eta_arr.size(); ielem++) {
        if (load_elem[ielem]) continue;
        eta_max = std::max(eta_max, eta_arr[ielem]);
    }

    for (int ielem = 0; ielem < eta_arr.size(); ielem++) {
        if (load_elem[ielem] || eta_arr[ielem] > threshold * eta_max) {
            marked_elem.push_back(ielem);
        }
    }
}