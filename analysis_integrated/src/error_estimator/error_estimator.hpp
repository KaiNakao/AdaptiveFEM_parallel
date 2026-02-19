#pragma once
#include <iostream>
#include <map>
#include <numeric>
#include <set>

#include "extend_mesh.hpp"
#include "mark_elem.hpp"
#include "mpi.h"
#include "posterior_error.hpp"
#include "reader.hpp"
#include "writer.hpp"

class Error_Estimator {
public:
    Error_Estimator(int myid, int numprocs);
    void prep();
    void exec(int iload, const std::vector<double> &uvec);
    void post();
    int myid, numprocs;
    std::vector<int> load_elem;
private:
    std::vector<std::vector<double>> coor;
    std::vector<std::vector<int>> cny;
    std::vector<int> matid_arr;
    std::vector<std::vector<double>> material;
    std::map<std::set<int>, std::set<int>> face_to_elems;
    int nelem, nnode, nmaterial, nload;
    std::map<int, std::vector<int>> import_index, export_index;
    int nelem_org;
    int nnode_org;
    std::vector<std::set<int>> face_arr;
    std::vector<double> eta_arr;
};