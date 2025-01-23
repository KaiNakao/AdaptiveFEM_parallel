#pragma once
#include <mpi.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void calc_projection(const int &myid,
                     const std::vector<std::vector<double>> &coor,
                     std::vector<std::vector<double>> &obs_point);

void read_loads(int &nload, std::vector<std::vector<double>> &load_arr);

void read_surf_cny(const int &myid, std::vector<std::vector<int>> &surf_cny);