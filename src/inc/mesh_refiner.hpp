#pragma once

// Dependency
#include "geo_util.hpp"
#include "reader.hpp"
#include "smoothing_scheme.hpp"
#include "refinement_scheme.hpp"

// STL
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>

class Refiner
{
public: // public member function
    Refiner();
    Refiner(const std::string &data_dir);
    ~Refiner();

    void executeRefinement();

    inline std::vector<std::vector<int>> getConnectivity() { return m_connectivity; }
    inline std::vector<std::vector<double>> getCoordinates() { return m_coordinates; }
    inline std::vector<int> getScheme() { return m_elem_to_scheme; }

public: // public member variable

private: // private member function
    // Switching Refinement Method 
    int switchScheme(int _elem_id);

    // Schemes
    void elementRefine(int _elem_id);
    void nodeSmoothing();

    // Mesh refinement tools
    void addPoints(int _elem_id);
    bool checkFlippable(int _elem_id);
    void performFlip(int _elem_id);
    void split14(int _tetra_id, std::vector<double> &_pt);
    void flip23(int _tetra_id);
    void replaceAdjElem(int _elem_id, int _old_id, int _new_id);
    int findVertexPos(int _target_tetra, int _known_tetra, int _vert_pos);

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
    std::string m_data_dir;
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

    // smoothing tool
    std::set<int> m_node_smooth;
    std::map<int, std::set<int>> m_adj_nodes;
    std::vector<int> m_elem_to_scheme;

    /// TODO: implementation of suitable data structure
    std::vector<std::vector<double>> m_points;
};