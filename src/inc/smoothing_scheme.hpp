#pragma once

//Dependecy
#include "geo_util.hpp"

// STL
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>

class SmoothingScheme
{
public: // public member function
    SmoothingScheme(const std::set<int> &elem_smooth,
                    std::vector<std::vector<int>> &connectivity,
                    std::vector<std::vector<double>> &coordinates);
    ~SmoothingScheme();

    void executeSmoothing(std::vector<std::vector<int>> &new_conn, 
                          std::vector<std::vector<double>> &new_coor);

private: // private member function

    // moving node tools
    void fetchNodes();
    std::vector<double> calculateMovingDirection(int _node_id);
    void search_adj_nodes();
    void moveNodes();
    double findMaxAspectRatio();
    bool compareMaxAspectRatio();

    // Other routines

    // Update mesh data
    void updateMeshData();

private: // private member variable
    // mesh info
    std::vector<std::vector<int>> m_connectivity;      // size: (number of elements, 4)
    std::vector<std::vector<double>> m_coordinates;    // size: (number of nodes, 3)
    std::vector<std::vector<int>> m_adj_elements;      // size: (number of elements, 4)
    std::vector<int> m_elem_smooth;                    // element ids to be smoothed

    // smoothing tool
    std::set<int> m_node_smooth;
    std::map<int, std::set<int>> m_adj_nodes;
};