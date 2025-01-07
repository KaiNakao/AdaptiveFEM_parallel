#pragma once
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "mpi.h"

void read_shape(const int &myid, int &nelem, int &nnode, int &nmaterial,
                int &nload);

void read_mesh(const int &myid, const int &nelem, const int &nnode,
               std::vector<std::vector<int>> &cny,
               std::vector<std::vector<double>> &coor,
               std::vector<int> &madid_arr);

void read_load_elem(const int &myid, std::vector<int> &load_elem, int nelem);

void read_material(const int &myid, const int &nmaterial,
                   std::vector<std::vector<double>> &material);

void read_neighbor(const int &myid,
                   std::map<int, std::vector<int>> &neighbor_map);

void read_displacement(const int &myid, const int &iload, const int &nnode,
                       std::vector<std::vector<double>> &displacement);

void read_obs_point(const int &myid,
                    std::vector<std::vector<double>> &obs_point);