#pragma once
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "mpi.h"
void extend_mesh(const int &myid,
                 const std::map<int, std::vector<int>> &neighbor_map,
                 std::vector<std::vector<int>> &cny,
                 std::vector<std::vector<double>> &coor,
                 std::vector<int> &matid_arr,
                 std::map<std::set<int>, std::set<int>> &face_to_elems,
                 std::map<int, std::vector<int>> &import_index,
                 std::map<int, std::vector<int>> &export_index);
