#pragma once
#include <iostream>
#include <vector>
#include <map>
#include "hdf_lib.h"

void write_shape(const int &myid, const int &nnode, const int &nelem, const int &nmaterial);

void write_mesh(const int &myid, 
                const std::vector<std::vector<int>> &cny,
                const std::vector<std::vector<double>> &coor,
                const std::vector<int> &matid_arr);

void write_neighbor(const int &myid, const std::map<int, std::vector<int>> &neighbor_nodes);
