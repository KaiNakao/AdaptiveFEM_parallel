#pragma once
#include <vector>
#include <iostream>

#include "mpi.h"
namespace nsp_error_estimator {
void mark_elem_to_refine(const std::vector<double> &eta_arr,
                         const std::vector<int> &load_elem,
                         std::vector<int> &marked_elem,
                         const double &threshold);
}
