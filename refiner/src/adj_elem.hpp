#pragma once
#include <set>
#include <map>
#include <vector>

void search_adj_element(const std::vector<std::vector<int>> &cny, 
                        const std::vector<std::vector<double>> &coor, 
                        std::vector<std::vector<int>> &adj_elem,
                        std::map<std::set<int>, std::vector<int>> &face_to_elems);
