#include "../inc/adj_elem.hpp"
#include <map>
#include <set>
#include <iostream>

void search_adj_element(const std::vector<std::vector<int>> &cny, 
                        const std::vector<std::vector<double>> &coor, 
                        std::vector<std::vector<int>> &adj_elem) {
    
    std::map<std::set<int>, std::vector<int>> face_to_elems;
    for (int ielem = 0; ielem < cny.size(); ielem++) {
        const std::vector<int> &elem = cny[ielem];

        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(elem[(iface + inode + 1) % 4]);
            }

            face_to_elems[face].push_back(ielem);
        }
    }

    for (int ielem = 0; ielem < cny.size(); ielem++) {
        const std::vector<int> &elem = cny[ielem];
        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(elem[(iface + inode + 1 ) % 4]);
            }
            if (face_to_elems[face].size() == 2) {
                for (int i = 0; i < 2; i++) {
                    if (face_to_elems[face][i] != ielem) {
                        adj_elem[ielem][iface] = face_to_elems[face][i];
                    }
                }
            } 
            else {
                adj_elem[ielem][iface] = -1;
            }
        }
    }

    // for (int ielem = 0; ielem < cny.size(); ielem++) {
    //     std::cout << "ielem: " << ielem << std::endl;
    //     for (int iface = 0; iface < 4; iface++) {
    //         std::cout << "iface: " << iface << " " << adj_elem[ielem][iface] << std::endl;
    //     }
    // }

}