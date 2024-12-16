#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "reader.hpp"

void read_shape(const std::string &data_dir, 
                int &nelem, int &nnode_linear, int &nnode_quad, int &nmaterial, int &nelem_marked) {

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
    ifs >> buf; // read "nmaterial"
    ifs >> nmaterial;
    ifs >> buf; // read "nelem_marked"
    ifs >> nelem_marked;

    std::cout << "nelem: " << nelem << std::endl;
    std::cout << "nnode_linear: " << nnode_linear << std::endl;
    std::cout << "nnode_quad: " << nnode_quad << std::endl;
    std::cout << "nmaterial: " << nmaterial << std::endl;
    std::cout << "nelem_marked: " << nelem_marked << std::endl;
}

void read_mesh(const std::string &data_dir, 
               const int &nelem, const int &nnode, 
               std::vector<std::vector<int>> &cny,
               std::vector<std::vector<double>> &coor,
               std::vector<int> &matid_arr) {
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
        matid_arr[ielem] = buf_cny[ielem * 11 + 10] - 1;
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
    if ((fp = fopen((data_dir + "displacement_quad.bin").c_str(), "r")) == NULL) {
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

void read_material(const std::string &data_dir, 
                   std::vector<std::vector<double>> &material) {
    std::vector<double> buf_material(3 * material.size());
    FILE *fp;
    if ((fp = fopen((data_dir + "material.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file material.bin" << std::endl;
        return;
    }
    fread(&buf_material[0], sizeof(double), 3 * material.size(), fp);
    fclose(fp);

    for (int imaterial = 0; imaterial < material.size(); imaterial++) {
        double vp = buf_material[imaterial * 3];
        double vs = buf_material[imaterial * 3 + 1];
        double rho = buf_material[imaterial * 3 + 2];
        double lam = rho * (vp * vp - 2.0 * vs * vs);
        double mu = rho * vs * vs;
        material[imaterial][0] = lam;
        material[imaterial][1] = mu;
    }

    for (int imaterial = 0; imaterial < material.size(); imaterial++) {
        std::cout << "material: " << imaterial << " ";
        for (int idim = 0; idim < 2; idim++) {
            std::cout << material[imaterial][idim] << " ";
        }
        std::cout << std::endl;
    }   
}

void read_load_elem(const std::string &data_dir, 
                    const int &nelem, 
                    std::vector<int> &load_elem) {
    FILE *fp;
    if ((fp = fopen((data_dir + "load_elem.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file load_elem.bin" << std::endl;
        return;
    }
    fread(&load_elem[0], sizeof(int), nelem, fp);
    fclose(fp);
}

void read_marked_elem(const std::string &data_dir,
                      const int &nelem_marked,
                      std::vector<int> &marked_elem) {
    FILE *fp;
    if ((fp = fopen((data_dir + "marked_elem_quad.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file marked_elem_quad.bin" << std::endl;
        return;
    }
    fread(&marked_elem[0], sizeof(int), nelem_marked, fp);
    fclose(fp);
}