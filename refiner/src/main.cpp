#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "reader.hpp"
#include "adj_elem.hpp"
#include "posterior_error.hpp"

int main() {
    std::string data_dir("/data6/itou/AFEM/data/analysis_result_iburi/");

    int nelem, nnode_linear, nnode_quad;
    read_shape(data_dir, nelem, nnode_linear, nnode_quad);

    std::vector<std::vector<int>> cny(nelem, std::vector<int>(10));
    std::vector<std::vector<double>> coor(nnode_quad, std::vector<double>(3));
    std::vector<std::vector<double>> displacement(nnode_quad, std::vector<double>(3));
    read_mesh(data_dir, nelem, nnode_quad, cny, coor);
    read_displacement(data_dir, nnode_quad, displacement);

    // print coordinates of element 0
    for (int inode = 0; inode < 10; inode++) {
        std::cout << "inode: " << inode << " ";
        for (int idim = 0; idim < 3; idim++) {
            std::cout << coor[cny[0][inode]][idim] << " ";
        }
        std::cout << std::endl;
    }

    std::vector<std::vector<int>> adj_elem(nelem, std::vector<int>(4));
    std::map<std::set<int>, std::vector<int>> face_to_elems;
    search_adj_element(cny, coor, adj_elem, face_to_elems);

    posterior_error_estimation(cny, coor, displacement, face_to_elems);
}
