#pragma once
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include "entity.hpp"
#include "geo_util.hpp"

class Refinement_scheme {
public: // public member function
    Refinement_scheme(const std::set<int> &elem_refine,
                      const std::vector<std::vector<int>> &connectivity,
                      const std::vector<std::vector<double>> &coordinates,
                      const std::vector<int> &matid_arr);
    ~Refinement_scheme();
    void executeRefinement_bisect(std::vector<std::vector<int>> &new_conn,
                                  std::vector<std::vector<double>> &new_coor,
                                  std::vector<int> &new_matid_arr,
                                  std::vector<int> &original);
    const std::map<int, Edge> &get_newnode_edge() const { return m_newnode_edge; }

private: // private member function
    std::unordered_map<int, int> m_matid_map;
    std::unordered_map<int, std::vector<double>> m_coor_map;
    const std::set<int> &m_elem_refine;
    const std::vector<std::vector<int>> &m_connectivity;

    void initialize_tets(std::set<Tetrahedron> &t, std::set<Tetrahedron> &s, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets);
    void local_refine(std::set<Tetrahedron> &t, std::set<Tetrahedron> &s, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets);
    void bisect_tets(std::set<Tetrahedron> &t, std::set<Tetrahedron> &s,
                     std::map<Edge, int> &edges_cut, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets);
    void bisect_tet(std::set<Tetrahedron> &t, const Tetrahedron &tet,
                    std::map<Edge, int> &edges_cut, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets);
    int find_tet_type(const Tetrahedron &tet);
    void refine_to_conformity(std::set<Tetrahedron> &t, std::set<Tetrahedron> &s,
                              std::map<Edge, int> &edges_cut, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets);
    void collect_tet_edges(const Tetrahedron &tet, std::vector<Edge> &edges);
    std::map<int, Edge> m_newnode_edge;
    void output_new_mesh(const std::set<Tetrahedron> &t, std::vector<std::vector<int>> &new_conn,
                         std::vector<std::vector<double>> &new_coor,
                         std::vector<int> &new_matid_arr, std::vector<int> &original);

    void tmp_output(const std::set<Tetrahedron> &s);
};
