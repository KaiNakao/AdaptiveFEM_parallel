#pragma once
#include <iostream>

#include "mpi.h"
#include "reader.hpp"
#include "writer.hpp"

#ifdef SURF
#include "surf.hpp"
#endif

class Obs_Displacement {
public:
    Obs_Displacement(int myid, int numprocs);
    void prep();
    void exec(int iload, const std::vector<double> &uvec);
    int myid, numprocs;
private:
    
    std::vector<std::vector<double>> coor;
    std::vector<std::vector<int>> cny;
    int nfound = 0;
    std::vector<int> iobs_arr;
    std::vector<int> ielem_arr;
    std::vector<std::vector<double>> unode;
    std::vector<std::vector<double>> result;
    std::vector<std::vector<double>> r_arr;
    std::vector<double> nvec, u, u_arr_buf;

    std::vector<double> u_arr_buf_total;
    std::vector<int> displs, nfound_arr;
    std::vector<int> iobs_arr_total;
    int nfound_total;
};
