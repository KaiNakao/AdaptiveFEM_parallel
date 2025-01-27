#pragma once
#include <mpi.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void write_nodal_force(const int &myid, const std::vector<int> &node_id_obs,
                       const std::vector<double> &fvec);