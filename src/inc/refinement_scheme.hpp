#pragma once
#include <vector>
#include <unordered_map>
#include <map>
#include <set>

class Refinement_scheme {
public: // public member function
    Refinement_scheme(const std::set<int> &elem_refine,
                      std::vector<std::vector<int>> &connectivity,
                      std::vector<std::vector<double>> &coordinates,
                      std::vector<int> &matid_arr,
                      std::map<std::set<int>, std::vector<int>> &face_to_elems);
    ~Refinement_scheme();
    void executeRefinement(std::vector<std::vector<int>> &new_conn, 
                           std::vector<std::vector<double>> &new_coor);

private: // private member function
    int m_max_elem_id;
    std::vector<std::vector<double>> m_points;
    std::unordered_map<int, std::vector<int>> m_elems_in_mesh;
    std::unordered_map<int, int> m_matid_map;
    std::unordered_map<int, std::vector<double>> m_coor_map;
    std::set<int> m_candidate_elems;
    std::map<std::set<int>, std::vector<int>> m_face_to_elems;

    void split14(int _tetra_id, std::vector<double> &_pt, std::set<int> &elems_for_flip);
    void performFlip(int _elem_id);
    bool checkFlippable(int _o_tetra_id);
    void flip23(int _o_tetra_id);
    int findVertexPos(int _target_tetra, int _known_tetra, int _vert_pos);
    bool checkConcave(const std::vector<std::vector<double>> &_tri,
                      const std::vector<double> &_pt);
};