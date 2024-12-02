#include <string>
#include <fstream>
#include <iostream>

#include "reader.hpp"

int main() {
    std::string data_dir("/data6/itou/AFEM/analysis_result/");

    int nelem, nnode_linear, nnode_quad;
    read_shape(data_dir, nelem, nnode_linear, nnode_quad);

}