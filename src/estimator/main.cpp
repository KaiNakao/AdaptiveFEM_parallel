// #include <string>
// #include <fstream>
// #include <iostream>
// #include <vector>

// #include "reader.hpp"
// #include "adj_elem.hpp"
// #include "posterior_error.hpp"
// #include "mark_elem.hpp"

// int main() {
//     std::string data_dir("/data6/itou/AFEM/data/analysis_result_iburi_large/");

//     int nelem, nnode_linear, nnode_quad, nmaterial;
//     read_shape(data_dir, nelem, nnode_linear, nnode_quad, nmaterial);

//     std::vector<std::vector<int>> cny(nelem, std::vector<int>(10));
//     std::vector<int> matid_arr(nelem);
//     std::vector<std::vector<double>> coor(nnode_quad, std::vector<double>(3));
//     std::vector<std::vector<double>> displacement(nnode_quad, std::vector<double>(3));
//     std::vector<int> load_elem(nelem);
//     read_mesh(data_dir, nelem, nnode_quad, cny, coor, matid_arr);
//     read_displacement(data_dir, nnode_quad, displacement);
//     read_load_elem(data_dir, nelem, load_elem);

//     for (int ielem = 0; ielem < nelem; ielem++) {
//         if (load_elem[ielem] == 1) {
//             std::cout << "load_elem: " << ielem << std::endl;
//         }
//     }

//     for (int inode = 0; inode < 10; inode++) {
//         std::cout << "inode: " << inode << " ";
//         for (int idim = 0; idim < 3; idim++) {
//             std::cout << std::scientific << displacement[inode][idim] << " ";
//         }
//         std::cout << std::endl;
//     }

//     std::vector<std::vector<double>> material(nmaterial, std::vector<double>(2));
//     read_material(data_dir, material);

//     // print coordinates of element 0
//     for (int inode = 0; inode < 10; inode++) {
//         std::cout << "inode: " << inode << " ";
//         for (int idim = 0; idim < 3; idim++) {
//             std::cout << coor[cny[0][inode]][idim] << " ";
//         }
//         std::cout << std::endl;
//     }

//     std::vector<std::vector<int>> adj_elem(nelem, std::vector<int>(4));
//     std::map<std::set<int>, std::vector<int>> face_to_elems;
//     search_adj_element(cny, coor, adj_elem, face_to_elems);

//     std::vector<double> eta_arr;
//     posterior_error_estimation(cny, coor, displacement, material, matid_arr, face_to_elems, eta_arr);

//     std::vector<int> marked_elem;
//     // mark_elem_to_refine(eta_arr, marked_elem);
//     mark_elem_to_refine(eta_arr, load_elem, marked_elem, 0.1);
//     int nelem_marked = marked_elem.size();
//     for (int ielem = 0; ielem < nelem_marked; ielem++) {
//         std::cout << "marked_elem: " << marked_elem[ielem] << " " << eta_arr[marked_elem[ielem]] << std::endl;
//     }

//     FILE *fp;
//     if ((fp = fopen((data_dir + "eta_quad.bin").c_str(), "w")) == NULL) {
//         std::cerr << "Error: cannot open file eta_quad.bin" << std::endl;
//         return 1;
//     }
//     fwrite(&eta_arr[0], sizeof(double), nelem, fp);
//     fclose(fp);

//     if ((fp = fopen((data_dir + "marked_elem_quad.bin").c_str(), "w")) == NULL) {
//         std::cerr << "Error: cannot open file marked_elem_quad.bin" << std::endl;
//         return 1;
//     }
//     fwrite(&marked_elem[0], sizeof(int), nelem_marked, fp);
//     fclose(fp);

//     std::ofstream ofs(data_dir + "shape.dat", std::ios::app);
//     if (!ofs) {
//         std::cerr << "Error: cannot open file shape.dat" << std::endl;
//         return 1;
//     }
//     ofs << "nelem_marked" << std::endl;
//     ofs << nelem_marked << std::endl;

// }
