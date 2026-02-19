#pragma once
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <mpi.h>
#include "entity.hpp"

class Refinement_scheme {
public: // public member function
    Refinement_scheme(const int myid, const int numprocs,
                      const std::set<int> &elem_refine,
                      const std::vector<std::vector<int>> &connectivity,
                      const std::vector<std::vector<double>> &coordinates,
                      const std::vector<int> &matid_arr,
                      const std::map<int, std::vector<int>> &neighbor_nodes);
    ~Refinement_scheme();
    void executeRefinement_bisect(std::vector<std::vector<int>> &new_conn,
                                  std::vector<std::vector<double>> &new_coor,
                                  std::vector<int> &new_matid_arr,
                                  std::vector<int> &original,
                                  std::map<int, std::vector<int>> &new_neighbor_nodes);

private: // private member function
    int m_myid;
    int m_numprocs;
    std::unordered_map<int, int> m_matid_map;
    std::unordered_map<int, std::vector<double>> m_coor_map;
    const std::set<int> &m_elem_refine;
    const std::vector<std::vector<int>> &m_connectivity;
    const std::map<int, std::vector<int>> &m_neighbor_nodes;

    void initialize_tets(std::set<Tetrahedron> &t, std::set<Tetrahedron> &s, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets);
    void local_refine(std::set<Tetrahedron> &t, std::set<Tetrahedron> &s, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets,
                      std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges);
    void bisect_tets(std::set<Tetrahedron> &t, std::set<Tetrahedron> &s,
                     std::map<Edge, int> &edges_cut, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets,
                     std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges);
    void bisect_tet(std::set<Tetrahedron> &t, const Tetrahedron &tet,
                    std::map<Edge, int> &edges_cut, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets);
    int find_tet_type(const Tetrahedron &tet);
    void refine_to_conformity(std::set<Tetrahedron> &t, std::set<Tetrahedron> &s,
                              std::map<Edge, int> &edges_cut, std::map<Edge, std::set<Tetrahedron>> &edge_to_tets,
                              std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges);
    void collect_tet_edges(const Tetrahedron &tet, std::vector<Edge> &edges);
    void output_new_mesh(const std::set<Tetrahedron> &t, std::vector<std::vector<int>> &new_conn,
                         std::vector<std::vector<double>> &new_coor,
                         std::vector<int> &new_matid_arr, std::vector<int> &original);
    void build_neighbor_edges(const std::set<Tetrahedron> &t,
                               std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges);
    void check_neighbor_edges(const std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges, 
                              const std::map<Edge, std::set<Tetrahedron>> &edge_to_tets);   
    void sync_edges_cut(std::map<Edge, int> &edges_cut,
                        std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges);
    void tmp_output(std::set<Tetrahedron> &s);
    void build_neighbor_nodes(const std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges,
                              std::map<int, std::vector<int>> &neighbor_nodes);
    void check_neighbor_nodes(const std::map<int, std::vector<int>> &neighbor_nodes);
};
