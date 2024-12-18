#pragma once

// Dependency
#include "geo_util.hpp"
#include "reader.hpp"
#include "refinement_scheme.hpp"

// STL
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>

class Refiner
{
public: // public member function
    Refiner();
    Refiner(const std::string &data_dir);
    ~Refiner();

    void executeRefinement();

public: // public member variable



private: // private member function
    // Switching Refinement Method 
    int switchScheme(int _elem_id);

    // moving node tools

    // Other routines

private: // private member variable
    std::string m_data_dir;
    // mesh info
    std::vector<std::vector<int>> m_connectivity;      // size: (number of elements, 4)
    std::vector<std::vector<double>> m_coordinates;    // size: (number of nodes, 3)
    std::vector<std::vector<int>> m_adj_elements;      // size: (number of elements, 4)
    std::vector<double> m_error;                       // size: number of elements
    std::vector<int> m_marked_elems_id;                // size: number of MARKED elements
    std::vector<int> m_matid_arr;                      // size: number of elements
    std::map<std::set<int>, std::vector<int>> m_face_to_elems; // face to elements

    // general info
    std::set<int> m_elem_refine;                       // element ids to be refined (node addition) <- not changed from initial data
    std::set<int> m_elem_smooth;                       // element ids to be smoothed

    // Mesh refinement tool
    std::vector<int> m_deleted_elems;
    std::unordered_map<int, std::vector<int>>  m_elements;
    std::set<int> m_candidate_elements;
    std::set<int> m_new_elements;

    /// TODO: implementation of suitable data structure
    std::vector<std::vector<double>> m_points;
};