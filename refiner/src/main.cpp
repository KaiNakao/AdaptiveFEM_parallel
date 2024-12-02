#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "reader.hpp"

int main() {
    std::string data_dir("/data6/itou/AFEM/analysis_result/");

    int nelem, nnode_linear, nnode_quad;
    read_shape(data_dir, nelem, nnode_linear, nnode_quad);

    std::vector<std::vector<int>> cny(nelem, std::vector<int>(10));
    std::vector<std::vector<double>> coor(nnode_quad, std::vector<double>(3));

    read_mesh(data_dir, nelem, nnode_quad, cny, coor);
}
