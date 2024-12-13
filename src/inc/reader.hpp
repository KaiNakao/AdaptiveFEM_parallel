#pragma once
#include <string>
#include <vector>

void read_shape(const std::string& data_dir, 
                int& nelem, int& nnode_linear, int& nnode_quad, int& nelem_marked);
                
void read_mesh(const std::string &data_dir, 
               const int &nelem, const int &nnode, 
               std::vector<std::vector<int>> &cny,
               std::vector<std::vector<double>> &coor);

void search_adj_element(const std::vector<std::vector<int>> &cny, 
                        const std::vector<std::vector<double>> &coor, 
                        std::vector<std::vector<int>> &adj_elem);

void read_displacement(const std::string &data_dir, 
                       const int &nnode, 
                       std::vector<std::vector<double>> &displacement);

void read_eta(const std::string &data_dir,
              const int &nelem,
              std::vector<double> &eta);

void read_marked_elem(const std::string &data_dir,
                      const int &nelem_marked,
                      std::vector<int> &marked_elem);