#include <string>
#include <fstream>
#include <iostream>
#include "reader.hpp"

void read_shape(const std::string &data_dir, 
                int &nelem, int &nnode_linear, int &nnode_quad) {

    std::ifstream ifs(data_dir + "shape.dat");
    if (!ifs) {
        std::cerr << "Error: cannot open file shape.dat" << std::endl;
        return;
    }

    std::string buf;
    ifs >> buf; // read "nelem"
    ifs >> nelem;
    ifs >> buf; // read "nnode_linear"
    ifs >> nnode_linear;
    ifs >> buf; // read "nnode_quad"
    ifs >> nnode_quad;

    std::cout << "nelem: " << nelem << std::endl;
    std::cout << "nnode_linear: " << nnode_linear << std::endl;
    std::cout << "nnode_quad: " << nnode_quad << std::endl;
}