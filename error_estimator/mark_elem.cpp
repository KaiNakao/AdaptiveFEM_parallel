#include "mark_elem.hpp"
void mark_elem_to_refine(const std::vector<double> &eta_arr,
                         const std::vector<int> &load_elem,
                         std::vector<int> &marked_elem,
                         const double &threshold) {
    double eta_max = 0.0, eta_max_local = 0.0;
    for (int ielem = 0; ielem < eta_arr.size(); ielem++) {
        if (load_elem[ielem])  {
            std::cout << "load elem: " << ielem << std::endl;
            continue;
        }
        eta_max = std::max(eta_max, eta_arr[ielem]);
    }
    eta_max_local = eta_max;
    MPI_Allreduce(&eta_max_local, &eta_max, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    marked_elem.resize(eta_arr.size());
    for (int ielem = 0; ielem < eta_arr.size(); ielem++) {
        if (load_elem[ielem] || eta_arr[ielem] > threshold * eta_max) {
            marked_elem.at(ielem) = 1;
        }
    }
}