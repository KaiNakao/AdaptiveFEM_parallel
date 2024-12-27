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
                    std::vector<std::vector<double>> &coordinates,
                    std::vector<int> &matid_arr,
                    std::map<std::set<int>, std::vector<int>> &face_to_elems,
                    std::map<int, std::vector<std::set<int>>> &node_to_faces,
                    std::vector<std::vector<int>> &adj_elems);
    ~SmoothingScheme();

    void executeSmoothing(std::vector<std::vector<int>> &new_conn, 
                          std::vector<std::vector<double>> &new_coor);

private: // private member function

    // moving node tools
    void fetchNodes();
    std::vector<double> calculateMovingDirection(int _node_id);
    void search_adj_nodes();
    void createMovableTabel();
    void moveNodes();
    double findMaxAspectRatio();
    double findEuclidianNorm();

    // Other routines

    // Update mesh data

private: // private member variable
    // mesh info
    std::vector<std::vector<int>> m_connectivity;      // size: (number of elements, 4)
    std::vector<std::vector<double>> m_coordinates;    // size: (number of nodes, 3)
    std::vector<std::vector<double>> m_coordinates_prev;
    std::vector<std::vector<int>> m_adj_elements;      // size: (number of elements, 4)
    std::vector<int> m_elem_smooth;                    // element ids to be smoothed
    std::unordered_map<int, int> m_matid_map;
    std::map<std::set<int>, std::vector<int>> m_face_to_elems;
    std::map<int, std::vector<std::set<int>>> m_node_to_faces;

    // smoothing tool
    std::set<int> m_node_smooth;
    std::map<int, bool> m_node_movable;
    std::map<int, std::set<int>> m_adj_nodes;
};