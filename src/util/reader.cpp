#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "../inc/reader.hpp"

void read_shape(const std::string &data_dir, 
                int &nelem, int &nnode_linear, int &nnode_quad, int &nelem_marked) {

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
    ifs >> buf; // read "nmaterials"
    ifs >> buf;
    ifs >> buf; // read "nelem_marked"
    ifs >> nelem_marked;

    std::cout << "nelem: " << nelem << std::endl;
    std::cout << "nnode_linear: " << nnode_linear << std::endl;
    std::cout << "nnode_quad: " << nnode_quad << std::endl;
    std::cout << "nelem_marked: " << nelem_marked << std::endl;
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
            cny[ielem][inode] = buf_cny[ielem * 11 + inode] - 1;
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

void read_displacement(const std::string &data_dir, 
                       const int &nnode, 
                       std::vector<std::vector<double>> &displacement) {
    std::vector<double> buf_displacement(3 * nnode);

    FILE *fp;
    // Read displacement
    if ((fp = fopen((data_dir + "displacement.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file displacement.bin" << std::endl;
        return;
    }
    fread(&buf_displacement[0], sizeof(double), 3 * nnode, fp);
    fclose(fp);
    for (int inode = 0; inode < nnode; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            displacement[inode][idim] = buf_displacement[inode * 3 + idim];
        }
    }
}

void read_eta(const std::string &data_dir,
              const int &nelem,
              std::vector<double> &eta) {
    std::vector<double> buf_eta(nelem);

    FILE *fp;
    // Read eta
    if ((fp = fopen((data_dir + "eta_quad.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file eta_quad.bin" << std::endl;
        return;
    }
    fread(&buf_eta[0], sizeof(double), nelem, fp);
    fclose(fp);
    for (int ielem = 0; ielem < nelem; ielem++) {
        eta[ielem] = buf_eta[ielem];
    }            
}

void read_marked_elem(const std::string &data_dir,
                      const int &nelem_marked,
                      std::vector<int> &marked_elem) {
    std::vector<int> buf_eta(nelem_marked);

    FILE *fp;
    // Read eta
    if ((fp = fopen((data_dir + "marked_elem_quad.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file maerked_elem_quad.bin" << std::endl;
        return;
    }
    fread(&buf_eta[0], sizeof(int), nelem_marked, fp);
    fclose(fp);
    for (int ielem = 0; ielem < nelem_marked; ielem++) {
        marked_elem[ielem] = buf_eta[ielem];
    }            
}