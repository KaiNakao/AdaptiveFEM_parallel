#pragma once

// Dependency
#include "reader.hpp"
#include "refinement_scheme.hpp"

// STL
#include <map>
#include <set>
#include <string>
#include <vector>

class Refiner {
public: // public member function
    Refiner();
    Refiner(const std::string &data_dir);
    ~Refiner();

    void executeRefinement();
    void readData();
    void writeNewMesh();

private: // private member function
    // mesh info
    std::string m_data_dir;
    std::vector<std::vector<int>> m_connectivity;      // size: (number of elements, 4)
    std::vector<std::vector<double>> m_coordinates;    // size: (number of nodes, 3)
    std::set<int> m_marked_elems_id;                   // size: number of marked elements
    std::vector<int> m_matid_arr;                      // size: number of elements

    // refined mesh info
    std::vector<std::vector<int>> m_new_connectivity;      // size: (number of refined elements, 4)
    std::vector<std::vector<double>> m_new_coordinates;    // size: (number of refined nodes, 3)
    std::vector<int> m_new_matid_arr;                      // size: number of refined elements
    std::vector<int> m_original;                          // size: number of refined elements, flag for original element
};
