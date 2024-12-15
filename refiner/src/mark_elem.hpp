#pragma once
#include <numeric>
#include <vector>
#include <algorithm>

// void mark_elem_to_refine(const std::vector<double> &eta_arr,
//                          std::vector<int> &marked_elem);
void mark_elem_to_refine(const std::vector<double> &eta_arr,
                         const std::vector<int> &load_elem,
                         std::vector<int> &marked_elem,
                         const double &threshold);