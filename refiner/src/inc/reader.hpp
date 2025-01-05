#pragma once
#include <string>
void read_shape(const std::string &data_dir, 
                int &nelem, int &nnode_linear, int &nnode_quad, int &nmaterial, int &nelem_marked);
                
void read_mesh(const std::string &data_dir, 
               const int &nelem, const int &nnode, 
               std::vector<std::vector<int>> &cny,
               std::vector<std::vector<double>> &coor,
               std::vector<int> &matid_arr);

void read_displacement(const std::string &data_dir, 
                       const int &nnode, 
                       std::vector<std::vector<double>> &displacement);

void read_material(const std::string &data_dir, 
                   std::vector<std::vector<double>> &material);

void read_load_elem(const std::string &data_dir, 
                    const int &nelem, 
                    std::vector<int> &load_elem);

void read_marked_elem(const std::string &data_dir,
                      const int &nelem_marked,
                      std::vector<int> &marked_elem);

void read_refinement_edge(const std::string &data_dir,
                          const int &nelem,
                          std::vector<std::vector<int>> &refinement_edge);

void read_marked_edge(const std::string &data_dir,
                      const int &nelem,
                      std::vector<std::vector<std::vector<int>>> &marked_edge);

void read_tet_flag(const std::string &data_dir,
                   const int &nelem,
                   std::vector<bool> &tet_flag);