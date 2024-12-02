#include <string>
#include <fstream>
#include <iostream>
#include <vector>

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

void read_mesh(const std::string &data_dir, 
               const int &nelem, const int &nnode, 
               std::vector<std::vector<int>> &cny,
               std::vector<std::vector<double>> &coor) {
    std::vector<int> buf_cny(11 * nelem);
    std::vector<double> buf_coor(3 * nnode);

    FILE *fp;
    // Read cny
    if ((fp = fopen((data_dir + "cny_quad.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file mesh.dat" << std::endl;
        return;
    }
    fread(&buf_cny[0], sizeof(int), 11 * nelem, fp);
    fclose(fp);
    for (int ielem = 0; ielem < nelem; ielem++) {
        for (int inode = 0; inode < 10; inode++) {
            cny[ielem][inode] = buf_cny[ielem * 11 + inode];
        }
    }

    // Read coor
    if ((fp = fopen((data_dir + "coor_quad.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file mesh.dat" << std::endl;
        return;
    }
    fread(&buf_coor[0], sizeof(double), 3 * nnode, fp);
    fclose(fp);
    for (int inode = 0; inode < nnode; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            coor[inode][idim] = buf_coor[inode * 3 + idim];
        }
    }
}
