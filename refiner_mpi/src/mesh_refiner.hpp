#pragma once

// Dependency
#include "reader.hpp"
#include "writer.hpp"   
#include "refinement_scheme.hpp"

// STL
#include <map>
#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <mpi.h>

class Refiner {
public: // public member function
    Refiner();
    Refiner(const int myid, const int numprocs);
    ~Refiner();

    void executeRefinement();
    void readData();
    void writeData();

private: // private member function
    int m_myid;
    int m_numprocs;
    int m_nmaterial;
    // mesh info
    std::vector<std::vector<int>> m_connectivity;      // size: (number of elements, 4)
    std::vector<std::vector<double>> m_coordinates;    // size: (number of nodes, 3)
    std::set<int> m_marked_elems_id;                   // size: number of marked elements
    std::vector<int> m_matid_arr;                      // size: number of elements
    std::map<int, std::vector<int>> m_neighbor_nodes; // size: number of nodes, list of neighbor nodes

    // refined mesh info
    std::vector<std::vector<int>> m_new_connectivity;      // size: (number of refined elements, 4)
    std::vector<std::vector<double>> m_new_coordinates;    // size: (number of refined nodes, 3)
    std::vector<int> m_new_matid_arr;                      // size: number of refined elements
    std::vector<int> m_original;                          // size: number of refined elements, flag for original element
    std::map<int, std::vector<int>> m_new_neighbor_nodes;
};
