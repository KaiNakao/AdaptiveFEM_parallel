#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <mpi.h>

void find_nodes_on_obs_points(const std::vector<std::vector<double>> &xnode_obs,
                              const std::vector<std::vector<int>> &cny,
                              const std::vector<std::vector<double>> &coor, 
                              const int &nnode, const int &nelem, const double &ds,
                              std::vector<int> &node_id_obs);
