#pragma once
#include <string>
void read_shape(const std::string& data_dir, 
                int& nelem, int& nnode_linear, int& nnode_quad);
                
void read_mesh(const std::string &data_dir, 
               const int &nelem, const int &nnode, 
               std::vector<std::vector<int>> &cny,
               std::vector<std::vector<double>> &coor);
