#include "posterior_error.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

#include "adj_elem.hpp"
#include "mark_elem.hpp"
#include "reader.hpp"

ErrorEstimator::ErrorEstimator(const std::string &data_dir) {
    int nelem, nnode_linear, nnode_quad, nmaterial, dummy;
    read_shape(data_dir, nelem, nnode_linear, nnode_quad, nmaterial, dummy);
    int nload;
    read_nload(nload);

    std::vector<std::vector<int>> cny(nelem, std::vector<int>(10));
    std::vector<int> matid_arr(nelem);
    std::vector<std::vector<double>> coor(nnode_quad, std::vector<double>(3));
    std::vector<std::vector<double>> displacement(nnode_quad,
                                                  std::vector<double>(3));
    std::vector<int> load_elem(nelem);
    read_mesh(data_dir, nelem, nnode_quad, cny, coor, matid_arr);
    read_load_elem(data_dir, nelem, load_elem);

    for (int ielem = 0; ielem < nelem; ielem++) {
        if (load_elem[ielem] == 1) {
            std::cout << "load_elem: " << ielem << std::endl;
        }
    }

    std::vector<std::vector<double>> material(nmaterial,
                                              std::vector<double>(2));
    read_material(data_dir, material);

    std::vector<std::vector<int>> adj_elem(nelem, std::vector<int>(4));
    std::map<std::set<int>, std::vector<int>> face_to_elems;
    std::map<int, std::vector<std::set<int>>> node_to_faces;
    search_adj_element(cny, coor, adj_elem, face_to_elems, node_to_faces);

    std::vector<double> eta_arr;
    posterior_error_estimation(cny, coor, displacement, material, matid_arr,
                               face_to_elems, nload, data_dir, eta_arr);

    std::vector<int> marked_elem;
    // mark_elem_to_refine(eta_arr, marked_elem);
    mark_elem_to_refine(eta_arr, load_elem, marked_elem, 0.1);
    int nelem_marked = marked_elem.size();
    for (int ielem = 0; ielem < nelem_marked; ielem++) {
        std::cout << "marked_elem: " << marked_elem[ielem] << " "
                  << eta_arr[marked_elem[ielem]] << std::endl;
    }

    FILE *fp;
    if ((fp = fopen((data_dir + "eta_quad.bin").c_str(), "w")) == NULL) {
        std::cerr << "Error: cannot open file eta_quad.bin" << std::endl;
        std::exit(1);
    }
    fwrite(&eta_arr[0], sizeof(double), nelem, fp);
    fclose(fp);

    if ((fp = fopen((data_dir + "marked_elem_quad.bin").c_str(), "w")) ==
        NULL) {
        std::cerr << "Error: cannot open file marked_elem_quad.bin"
                  << std::endl;
        std::exit(1);
    }
    for (int ielem = 0; ielem < nelem_marked; ielem++) {
        marked_elem[ielem] += 1;
    }
    fwrite(&marked_elem[0], sizeof(int), nelem_marked, fp);
    fclose(fp);

    std::ofstream ofs(data_dir + "shape.dat");
    if (!ofs) {
        std::cerr << "Error: cannot open file shape.dat" << std::endl;
        std::exit(1);
    }
    ofs << "nelem" << std::endl;
    ofs << nelem << std::endl;
    ofs << "nnode_linear" << std::endl;
    ofs << nnode_linear << std::endl;
    ofs << "nnode_quad" << std::endl;
    ofs << nnode_quad << std::endl;
    ofs << "nmaterial" << std::endl;
    ofs << nmaterial << std::endl;
    ofs << "nelem_marked" << std::endl;
    ofs << nelem_marked << std::endl;
    ofs.close();
}

ErrorEstimator::~ErrorEstimator() {}

void ErrorEstimator::calc_rk(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &displacement,
    const std::vector<std::vector<double>> &material,
    const std::vector<int> &matid_arr, std::vector<double> &rk_arr) {
    rk_arr.resize(cny.size());

    // dndrdr[i][j][k] = dn_i/(dr_j dr_k)
    std::vector<std::vector<std::vector<double>>> dndrdr = calc_dndrdr();

    for (int ielem = 0; ielem < cny.size(); ielem++) {
        std::vector<int> elem_id = cny[ielem];
        std::vector<std::vector<double>> xnode(10, std::vector<double>(3));
        std::vector<double> uvec(30);
        for (int inode = 0; inode < 10; inode++) {
            for (int idim = 0; idim < 3; idim++) {
                xnode[inode][idim] = coor[elem_id[inode]][idim];
                uvec[inode * 3 + idim] = displacement[elem_id[inode]][idim];
            }
        }

        int matid = matid_arr[ielem];
        double lam = material[matid][0];
        double mu = material[matid][1];
        std::vector<std::vector<double>> dmat(6, std::vector<double>(6));
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                dmat[i][j] = 0.0;
            }
        }
        dmat[0][0] = 2.0 * mu + lam;
        dmat[0][1] = lam;
        dmat[0][2] = lam;
        dmat[1][0] = lam;
        dmat[1][1] = 2.0 * mu + lam;
        dmat[1][2] = lam;
        dmat[2][0] = lam;
        dmat[2][1] = lam;
        dmat[2][2] = 2.0 * mu + lam;
        dmat[3][3] = mu;
        dmat[4][4] = mu;
        dmat[5][5] = mu;

        // dxdr[i][j] = dx_i/dr_j
        std::vector<std::vector<double>> dxdr = calc_dxdr(xnode);
        // drdx[i][j] = dr_i/dx_j
        std::vector<std::vector<double>> drdx = calc_drdx(dxdr);
        // volume of the element
        double vol = calc_volume(dxdr);
        // dndxdx[i][j][k] = dn_i/(dx_j dx_k)
        std::vector<std::vector<std::vector<double>>> dndxdx =
            calc_dndxdx(dndrdr, drdx);

        // internal residual
        std::vector<double> div_sigma(3);
        for (int idim = 0; idim < 3; idim++) {
            div_sigma[idim] = 0.0;
        }

        for (int idim = 0; idim < 3; idim++) {
            std::vector<std::vector<double>> bmat_grad(6,
                                                       std::vector<double>(30));
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 30; j++) {
                    bmat_grad[i][j] = 0.0;
                }
            }
            for (int inode = 0; inode < 10; inode++) {
                bmat_grad[0][inode * 3] = dndxdx[inode][0][idim];
                bmat_grad[1][inode * 3 + 1] = dndxdx[inode][1][idim];
                bmat_grad[2][inode * 3 + 2] = dndxdx[inode][2][idim];
                bmat_grad[3][inode * 3] = dndxdx[inode][1][idim];
                bmat_grad[3][inode * 3 + 1] = dndxdx[inode][0][idim];
                bmat_grad[4][inode * 3 + 1] = dndxdx[inode][2][idim];
                bmat_grad[4][inode * 3 + 2] = dndxdx[inode][1][idim];
                bmat_grad[5][inode * 3] = dndxdx[inode][2][idim];
                bmat_grad[5][inode * 3 + 2] = dndxdx[inode][0][idim];
            }
            std::vector<double> sigma_grad(6);
            for (int i = 0; i < 6; i++) {
                double tmp = 0.0;
                for (int j = 0; j < 6; j++) {
                    for (int k = 0; k < 30; k++) {
                        tmp += dmat[i][j] * bmat_grad[j][k] * uvec[k];
                    }
                }
                sigma_grad[i] = tmp;
            }
            if (idim == 0) {
                div_sigma[0] += sigma_grad[0];  // sigma_xx
                div_sigma[1] += sigma_grad[3];  // sigma_xy
                div_sigma[2] += sigma_grad[5];  // sigma_xz

            } else if (idim == 1) {
                div_sigma[0] += sigma_grad[3];  // sigma_xy
                div_sigma[1] += sigma_grad[1];  // sigma_yy
                div_sigma[2] += sigma_grad[4];  // sigma_yz
            } else if (idim == 2) {
                div_sigma[0] += sigma_grad[5];  // sigma_xz
                div_sigma[1] += sigma_grad[4];  // sigma_yz
                div_sigma[2] += sigma_grad[2];  // sigma_zz
            }

            // 断層ではmoment tensorによる力
        }

        // ノルムを要素内積分
        double rk = 0.0;
        for (int idim = 0; idim < 3; idim++) {
            rk += pow(div_sigma[idim], 2);
        }
        rk *= vol;

        rk_arr[ielem] = rk;
    }
}

void ErrorEstimator::calc_re(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &displacement,
    const std::vector<std::vector<double>> &material,
    const std::vector<int> &matid_arr,
    const std::map<std::set<int>, std::vector<int>> &face_to_elems,
    std::map<std::set<int>, double> &re_map) {
    std::map<std::vector<int>, std::vector<std::vector<double>>>
        gauss_point_dict = set_gauss_point_dict();

    for (const auto &p : face_to_elems) {
        std::set<int> face = p.first;

        std::vector<int> elems = p.second;
        double re, beta;
        if (elems.size() == 1) {
            // 境界
            re = 0.;
        } else {
            // stress jump
            re = calc_stress_jump(cny, coor, displacement, material, matid_arr,
                                  face, elems, gauss_point_dict);
        }
        re_map[face] = re;
    }
}

double ErrorEstimator::calc_stress_jump(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &displacement,
    const std::vector<std::vector<double>> &material,
    const std::vector<int> &matid_arr, const std::set<int> &face,
    const std::vector<int> &elems,
    const std::map<std::vector<int>, std::vector<std::vector<double>>>
        &gauss_point_dict) {
    std::vector<std::vector<std::vector<double>>> gauss_tractions(
        2, std::vector<std::vector<double>>(3, std::vector<double>(3)));

    std::vector<double> normal_vec = calc_normal_vec(coor, face);
    double area = calc_area(coor, face);

    for (int ielem = 0; ielem < 2; ielem++) {
        int elem_id = elems[ielem];
        std::vector<int> inode_face_arr;
        for (int inode : face) {
            for (int jnode = 0; jnode < 4; jnode++) {
                if (cny[elem_id][jnode] == inode) {
                    inode_face_arr.push_back(jnode);
                }
            }
        }

        // {r1, r2, r3, weight}
        std::vector<std::vector<double>> gauss_points =
            gauss_point_dict.at(inode_face_arr);

        std::vector<std::vector<double>> xnode(10, std::vector<double>(3));
        std::vector<double> uvec(30);
        for (int inode = 0; inode < 10; inode++) {
            int node_id = cny[elem_id][inode];
            for (int idim = 0; idim < 3; idim++) {
                xnode[inode][idim] = coor[node_id][idim];
                uvec[inode * 3 + idim] = displacement[node_id][idim];
            }
        }

        int matid = matid_arr[elem_id];
        double lam = material[matid][0];
        double mu = material[matid][1];
        std::vector<std::vector<double>> dmat(6, std::vector<double>(6));
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                dmat[i][j] = 0.0;
            }
        }
        dmat[0][0] = 2.0 * mu + lam;
        dmat[0][1] = lam;
        dmat[0][2] = lam;
        dmat[1][0] = lam;
        dmat[1][1] = 2.0 * mu + lam;
        dmat[1][2] = lam;
        dmat[2][0] = lam;
        dmat[2][1] = lam;
        dmat[2][2] = 2.0 * mu + lam;
        dmat[3][3] = mu;
        dmat[4][4] = mu;
        dmat[5][5] = mu;

        // dxdr[i][j] = dx_i/dr_j
        std::vector<std::vector<double>> dxdr = calc_dxdr(xnode);
        // drdx[i][j] = dr_i/dx_j
        std::vector<std::vector<double>> drdx = calc_drdx(dxdr);

        for (int ipoint = 0; ipoint < 3; ipoint++) {
            double r1 = gauss_points[ipoint][0];
            double r2 = gauss_points[ipoint][1];
            double r3 = gauss_points[ipoint][2];
            double weight = gauss_points[ipoint][3];
            double l0 = 1.0 - r1 - r2 - r3, l1 = r1, l2 = r2, l3 = r3;

            std::vector<std::vector<double>> dndr = calc_dndr(l0, l1, l2, l3);
            std::vector<std::vector<double>> dndx = calc_dndx(drdx, dndr);

            std::vector<std::vector<double>> bmat(6, std::vector<double>(30));
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 30; j++) {
                    bmat[i][j] = 0.0;
                }
            }
            for (int inode = 0; inode < 10; inode++) {
                bmat[0][inode * 3] = dndx[inode][0];
                bmat[1][inode * 3 + 1] = dndx[inode][1];
                bmat[2][inode * 3 + 2] = dndx[inode][2];
                bmat[3][inode * 3] = dndx[inode][1];
                bmat[3][inode * 3 + 1] = dndx[inode][0];
                bmat[4][inode * 3 + 1] = dndx[inode][2];
                bmat[4][inode * 3 + 2] = dndx[inode][1];
                bmat[5][inode * 3] = dndx[inode][2];
                bmat[5][inode * 3 + 2] = dndx[inode][0];
            }

            std::vector<double> sigma(6);
            for (int i = 0; i < 6; i++) {
                double tmp = 0.0;
                for (int j = 0; j < 6; j++) {
                    for (int k = 0; k < 30; k++) {
                        tmp += dmat[i][j] * bmat[j][k] * uvec[k];
                    }
                }
                sigma[i] = tmp;
            }

            // traction
            gauss_tractions[ielem][ipoint][0] =
                (sigma[0] * normal_vec[0] + sigma[3] * normal_vec[1] +
                 sigma[5] * normal_vec[2]) *
                weight;
            gauss_tractions[ielem][ipoint][1] =
                (sigma[3] * normal_vec[0] + sigma[1] * normal_vec[1] +
                 sigma[4] * normal_vec[2]) *
                weight;
            gauss_tractions[ielem][ipoint][2] =
                (sigma[5] * normal_vec[0] + sigma[4] * normal_vec[1] +
                 sigma[2] * normal_vec[2]) *
                weight;
        }
    }

    double re = 0.0;
    for (int ipoint = 0; ipoint < 3; ipoint++) {
        for (int idim = 0; idim < 3; idim++) {
            re += pow(gauss_tractions[0][ipoint][idim] -
                          gauss_tractions[1][ipoint][idim],
                      2);
        }
    }
    re *= area;
    return re;
}

std::vector<std::vector<double>> ErrorEstimator::calc_dndr(double l0, double l1,
                                                           double l2,
                                                           double l3) {
    std::vector<std::vector<double>> dndr(10, std::vector<double>(3));
    dndr[0] = {1.0 - 4.0 * l0, 1.0 - 4.0 * l0, 1.0 - 4.0 * l0};
    dndr[1] = {-1.0 + 4.0 * l0, 0.0, 0.0};
    dndr[2] = {0.0, -1.0 + 4.0 * l2, 0.0};
    dndr[3] = {0.0, 0.0, -1.0 + 4.0 * l3};
    dndr[4] = {-4.0 * l0 + 4.0 * l0, -4.0 * l0, -4.0 * l0};
    dndr[5] = {4.0 * l2, 4.0 * l0, 0.0};
    dndr[6] = {-4.0 * l2, -4.0 * l2 + 4.0 * l0, -4.0 * l2};
    dndr[7] = {-4.0 * l3, -4.0 * l3, -4.0 * l3 + 4.0 * l0};
    dndr[8] = {4.0 * l3, 0.0, 4.0 * l0};
    dndr[9] = {0.0, 4.0 * l3, 4.0 * l2};
    return dndr;
}

std::vector<std::vector<double>> ErrorEstimator::calc_dndx(
    const std::vector<std::vector<double>> &drdx,
    const std::vector<std::vector<double>> &dndr) {
    std::vector<std::vector<double>> dndx(10, std::vector<double>(3));
    for (int inode = 0; inode < 10; inode++) {
        for (int ix = 0; ix < 3; ix++) {
            dndx[inode][ix] = 0.0;
            for (int ir = 0; ir < 3; ir++) {
                dndx[inode][ix] += dndr[inode][ir] * drdx[ir][ix];
            }
        }
    }
    return dndx;
}

std::map<std::vector<int>, std::vector<std::vector<double>>>
ErrorEstimator::set_gauss_point_dict() {
    std::map<std::vector<int>, std::vector<std::vector<double>>>
        gauss_point_dict;
    std::map<int, std::vector<double>> node_to_point;
    // plane r3 = 0
    node_to_point[0] = {0.5, 0.5, 0.0, 1.0 / 3.0};
    node_to_point[1] = {0.0, 0.5, 0.0, 1.0 / 3.0};
    node_to_point[2] = {0.5, 0.0, 0.0, 1.0 / 3.0};
    gauss_point_dict[{0, 1, 2}] = {node_to_point[0], node_to_point[1],
                                   node_to_point[2]};
    gauss_point_dict[{0, 2, 1}] = {node_to_point[0], node_to_point[2],
                                   node_to_point[1]};
    gauss_point_dict[{1, 0, 2}] = {node_to_point[1], node_to_point[0],
                                   node_to_point[2]};
    gauss_point_dict[{1, 2, 0}] = {node_to_point[1], node_to_point[2],
                                   node_to_point[0]};
    gauss_point_dict[{2, 0, 1}] = {node_to_point[2], node_to_point[0],
                                   node_to_point[1]};
    gauss_point_dict[{2, 1, 0}] = {node_to_point[2], node_to_point[1],
                                   node_to_point[0]};
    node_to_point.clear();

    // plane r1 + r2 + r3 = 1
    node_to_point[1] = {0.0, 0.5, 0.5, 1.0 / 3.0};
    node_to_point[2] = {0.5, 0.0, 0.5, 1.0 / 3.0};
    node_to_point[3] = {0.5, 0.5, 0.0, 1.0 / 3.0};
    gauss_point_dict[{1, 2, 3}] = {node_to_point[1], node_to_point[2],
                                   node_to_point[3]};
    gauss_point_dict[{1, 3, 2}] = {node_to_point[1], node_to_point[3],
                                   node_to_point[2]};
    gauss_point_dict[{2, 1, 3}] = {node_to_point[2], node_to_point[1],
                                   node_to_point[3]};
    gauss_point_dict[{2, 3, 1}] = {node_to_point[2], node_to_point[3],
                                   node_to_point[1]};
    gauss_point_dict[{3, 1, 2}] = {node_to_point[3], node_to_point[1],
                                   node_to_point[2]};
    gauss_point_dict[{3, 2, 1}] = {node_to_point[3], node_to_point[2],
                                   node_to_point[1]};
    node_to_point.clear();

    // plane r1 = 0
    node_to_point[0] = {0.0, 0.5, 0.5, 1.0 / 3.0};
    node_to_point[2] = {0.0, 0.0, 0.5, 1.0 / 3.0};
    node_to_point[3] = {0.0, 0.5, 0.0, 1.0 / 3.0};
    gauss_point_dict[{0, 2, 3}] = {node_to_point[0], node_to_point[2],
                                   node_to_point[3]};
    gauss_point_dict[{0, 3, 2}] = {node_to_point[0], node_to_point[3],
                                   node_to_point[2]};
    gauss_point_dict[{2, 0, 3}] = {node_to_point[2], node_to_point[0],
                                   node_to_point[3]};
    gauss_point_dict[{2, 3, 0}] = {node_to_point[2], node_to_point[3],
                                   node_to_point[0]};
    gauss_point_dict[{3, 0, 2}] = {node_to_point[3], node_to_point[0],
                                   node_to_point[2]};
    gauss_point_dict[{3, 2, 0}] = {node_to_point[3], node_to_point[2],
                                   node_to_point[0]};
    node_to_point.clear();

    // plane r2 = 0
    node_to_point[0] = {0.5, 0.0, 0.5, 1.0 / 3.0};
    node_to_point[1] = {0.0, 0.0, 0.5, 1.0 / 3.0};
    node_to_point[3] = {0.5, 0.0, 0.0, 1.0 / 3.0};
    gauss_point_dict[{0, 1, 3}] = {node_to_point[0], node_to_point[1],
                                   node_to_point[3]};
    gauss_point_dict[{0, 3, 1}] = {node_to_point[0], node_to_point[3],
                                   node_to_point[1]};
    gauss_point_dict[{1, 0, 3}] = {node_to_point[1], node_to_point[0],
                                   node_to_point[3]};
    gauss_point_dict[{1, 3, 0}] = {node_to_point[1], node_to_point[3],
                                   node_to_point[0]};
    gauss_point_dict[{3, 0, 1}] = {node_to_point[3], node_to_point[0],
                                   node_to_point[1]};
    gauss_point_dict[{3, 1, 0}] = {node_to_point[3], node_to_point[1],
                                   node_to_point[0]};
    node_to_point.clear();

    return gauss_point_dict;
}

double ErrorEstimator::calc_hk(const std::vector<std::vector<double>> &xnode) {
    // length of each edge
    std::vector<double> edge_length(6);
    for (int iedge = 0; iedge < 6; iedge++) {
        int inode1, inode2;
        if (iedge == 0) {
            inode1 = 0;
            inode2 = 1;
        } else if (iedge == 1) {
            inode1 = 0;
            inode2 = 2;
        } else if (iedge == 2) {
            inode1 = 0;
            inode2 = 3;
        } else if (iedge == 3) {
            inode1 = 1;
            inode2 = 2;
        } else if (iedge == 4) {
            inode1 = 1;
            inode2 = 3;
        } else if (iedge == 5) {
            inode1 = 2;
            inode2 = 3;
        }
        for (int idim = 0; idim < 3; idim++) {
            edge_length[iedge] +=
                pow(xnode[inode1][idim] - xnode[inode2][idim], 2);
        }
        edge_length[iedge] = sqrt(edge_length[iedge]);
    }
    // longest edge
    double hk = *std::max_element(edge_length.begin(), edge_length.end());
    return hk;
}

std::vector<double> ErrorEstimator::calc_normal_vec(
    const std::vector<std::vector<double>> &coor, const std::set<int> &face) {
    std::vector<int> node_id_arr(face.begin(), face.end());
    std::vector<double> normal_vec(3);
    std::vector<double> vec1(3), vec2(3);
    for (int idim = 0; idim < 3; idim++) {
        vec1[idim] = coor[node_id_arr[1]][idim] - coor[node_id_arr[0]][idim];
        vec2[idim] = coor[node_id_arr[2]][idim] - coor[node_id_arr[0]][idim];
    }
    normal_vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    normal_vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    normal_vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    // normalize
    double norm = sqrt(pow(normal_vec[0], 2) + pow(normal_vec[1], 2) +
                       pow(normal_vec[2], 2));
    for (int idim = 0; idim < 3; idim++) {
        normal_vec[idim] /= norm;
    }
    return normal_vec;
}

double ErrorEstimator::calc_area(const std::vector<std::vector<double>> &coor,
                                 const std::set<int> &face) {
    std::vector<int> node_id_arr(face.begin(), face.end());
    std::vector<double> vec1(3), vec2(3);
    for (int idim = 0; idim < 3; idim++) {
        vec1[idim] = coor[node_id_arr[1]][idim] - coor[node_id_arr[0]][idim];
        vec2[idim] = coor[node_id_arr[2]][idim] - coor[node_id_arr[0]][idim];
    }
    std::vector<double> cross_product(3);
    cross_product[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    cross_product[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    cross_product[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    double area = sqrt(pow(cross_product[0], 2) + pow(cross_product[1], 2) +
                       pow(cross_product[2], 2)) /
                  2.0;
    return area;
}

double ErrorEstimator::calc_he(const std::vector<std::vector<double>> &coor,
                               const std::set<int> &face) {
    std::vector<int> node_id_arr(face.begin(), face.end());
    // length of each edge
    std::vector<double> edge_length(3);
    for (int iedge = 0; iedge < 3; iedge++) {
        int inode1, inode2;
        if (iedge == 0) {
            inode1 = 0;
            inode2 = 1;
        } else if (iedge == 1) {
            inode1 = 0;
            inode2 = 2;
        } else if (iedge == 2) {
            inode1 = 1;
            inode2 = 2;
        }
        for (int idim = 0; idim < 3; idim++) {
            edge_length[iedge] += pow(coor[node_id_arr[inode1]][idim] -
                                          coor[node_id_arr[inode2]][idim],
                                      2);
        }
        edge_length[iedge] = sqrt(edge_length[iedge]);
    }

    // longest edge
    double he = *std::max_element(edge_length.begin(), edge_length.end());
    return he;
}

std::vector<std::vector<std::vector<double>>> ErrorEstimator::calc_dndrdr() {
    std::vector<std::vector<std::vector<double>>> dndrdr(
        10, std::vector<std::vector<double>>(3, std::vector<double>(3)));

    dndrdr[0][0][0] = 4.0;
    dndrdr[0][0][1] = 4.0;
    dndrdr[0][0][2] = 4.0;
    dndrdr[0][1][0] = 4.0;
    dndrdr[0][1][1] = 4.0;
    dndrdr[0][1][2] = 4.0;
    dndrdr[0][2][0] = 4.0;
    dndrdr[0][2][1] = 4.0;
    dndrdr[0][2][2] = 4.0;

    dndrdr[1][0][0] = 4.0;
    dndrdr[1][0][1] = 0.0;
    dndrdr[1][0][2] = 0.0;
    dndrdr[1][1][0] = 0.0;
    dndrdr[1][1][1] = 0.0;
    dndrdr[1][1][2] = 0.0;
    dndrdr[1][2][0] = 0.0;
    dndrdr[1][2][1] = 0.0;
    dndrdr[1][2][2] = 0.0;

    dndrdr[2][0][0] = 0.0;
    dndrdr[2][0][1] = 0.0;
    dndrdr[2][0][2] = 0.0;
    dndrdr[2][1][0] = 0.0;
    dndrdr[2][1][1] = 4.0;
    dndrdr[2][1][2] = 0.0;
    dndrdr[2][2][0] = 0.0;
    dndrdr[2][2][1] = 0.0;
    dndrdr[2][2][2] = 0.0;

    dndrdr[3][0][0] = 0.0;
    dndrdr[3][0][1] = 0.0;
    dndrdr[3][0][2] = 0.0;
    dndrdr[3][1][0] = 0.0;
    dndrdr[3][1][1] = 0.0;
    dndrdr[3][1][2] = 0.0;
    dndrdr[3][2][0] = 0.0;
    dndrdr[3][2][1] = 0.0;
    dndrdr[3][2][2] = 4.0;

    dndrdr[4][0][0] = -8.0;
    dndrdr[4][0][1] = -4.0;
    dndrdr[4][0][2] = -4.0;
    dndrdr[4][1][0] = -4.0;
    dndrdr[4][1][1] = 0.0;
    dndrdr[4][1][2] = 4.0;
    dndrdr[4][2][0] = -4.0;
    dndrdr[4][2][1] = 0.0;
    dndrdr[4][2][2] = 0.0;

    dndrdr[5][0][0] = 0.0;
    dndrdr[5][0][1] = 4.0;
    dndrdr[5][0][2] = 0.0;
    dndrdr[5][1][0] = 4.0;
    dndrdr[5][1][1] = 0.0;
    dndrdr[5][1][2] = 0.0;
    dndrdr[5][2][0] = 0.0;
    dndrdr[5][2][1] = 0.0;
    dndrdr[5][2][2] = 0.0;

    dndrdr[6][0][0] = 0.0;
    dndrdr[6][0][1] = -4.0;
    dndrdr[6][0][2] = 0.0;
    dndrdr[6][1][0] = -4.0;
    dndrdr[6][1][1] = -8.0;
    dndrdr[6][1][2] = -4.0;
    dndrdr[6][2][0] = 0.0;
    dndrdr[6][2][1] = -4.0;
    dndrdr[6][2][2] = 0.0;

    dndrdr[7][0][0] = 0.0;
    dndrdr[7][0][1] = 0.0;
    dndrdr[7][0][2] = -4.0;
    dndrdr[7][1][0] = 0.0;
    dndrdr[7][1][1] = 0.0;
    dndrdr[7][1][2] = -4.0;
    dndrdr[7][2][0] = -4.0;
    dndrdr[7][2][1] = -4.0;
    dndrdr[7][2][2] = -8.0;

    dndrdr[8][0][0] = 0.0;
    dndrdr[8][0][1] = 0.0;
    dndrdr[8][0][2] = 4.0;
    dndrdr[8][1][0] = 0.0;
    dndrdr[8][1][1] = 0.0;
    dndrdr[8][1][2] = 0.0;
    dndrdr[8][2][0] = 4.0;
    dndrdr[8][2][1] = 0.0;
    dndrdr[8][2][2] = 0.0;

    dndrdr[9][0][0] = 0.0;
    dndrdr[9][0][1] = 0.0;
    dndrdr[9][0][2] = 0.0;
    dndrdr[9][1][0] = 0.0;
    dndrdr[9][1][1] = 0.0;
    dndrdr[9][1][2] = 4.0;
    dndrdr[9][2][0] = 0.0;
    dndrdr[9][2][1] = 4.0;
    dndrdr[9][2][2] = 0.0;

    return dndrdr;
}

std::vector<std::vector<double>> ErrorEstimator::calc_dxdr(
    const std::vector<std::vector<double>> &xnode) {
    std::vector<std::vector<double>> dxdr(3, std::vector<double>(3));
    double x1 = xnode[0][0], y1 = xnode[0][1], z1 = xnode[0][2];
    double x2 = xnode[1][0], y2 = xnode[1][1], z2 = xnode[1][2];
    double x3 = xnode[2][0], y3 = xnode[2][1], z3 = xnode[2][2];
    double x4 = xnode[3][0], y4 = xnode[3][1], z4 = xnode[3][2];
    dxdr[0][0] = -x1 + x2;
    dxdr[0][1] = -x1 + x3;
    dxdr[0][2] = -x1 + x4;
    dxdr[1][0] = -y1 + y2;
    dxdr[1][1] = -y1 + y3;
    dxdr[1][2] = -y1 + y4;
    dxdr[2][0] = -z1 + z2;
    dxdr[2][1] = -z1 + z3;
    dxdr[2][2] = -z1 + z4;
    return dxdr;
}

std::vector<std::vector<std::vector<double>>> ErrorEstimator::calc_dndxdx(
    const std::vector<std::vector<std::vector<double>>> &dndrdr,
    const std::vector<std::vector<double>> &drdx) {
    std::vector<std::vector<std::vector<double>>> dndxdx(
        10, std::vector<std::vector<double>>(3, std::vector<double>(3)));
    for (int inode = 0; inode < 10; inode++) {
        for (int ix = 0; ix < 3; ix++) {
            for (int jx = 0; jx < 3; jx++) {
                double tmp = 0.0;
                for (int ir = 0; ir < 3; ir++) {
                    for (int jr = 0; jr < 3; jr++) {
                        tmp +=
                            dndrdr[inode][ir][jr] * drdx[ir][ix] * drdx[jr][jx];
                    }
                }
                dndxdx[inode][ix][jx] = tmp;
            }
        }
    }
    return dndxdx;
}

double ErrorEstimator::calc_volume(std::vector<std::vector<double>> &dxdr) {
    double detj = dxdr[0][0] * dxdr[1][1] * dxdr[2][2] +
                  dxdr[1][0] * dxdr[2][1] * dxdr[0][2] +
                  dxdr[2][0] * dxdr[0][1] * dxdr[1][2] -
                  dxdr[0][0] * dxdr[2][1] * dxdr[1][2] -
                  dxdr[1][0] * dxdr[0][1] * dxdr[2][2] -
                  dxdr[2][0] * dxdr[1][1] * dxdr[0][2];
    double vol = detj / 6.0;
    return vol;
}

std::vector<std::vector<double>> ErrorEstimator::calc_drdx(
    const std::vector<std::vector<double>> &dxdr) {
    std::vector<std::vector<double>> drdx(3, std::vector<double>(3));
    double detj = dxdr[0][0] * dxdr[1][1] * dxdr[2][2] +
                  dxdr[1][0] * dxdr[2][1] * dxdr[0][2] +
                  dxdr[2][0] * dxdr[0][1] * dxdr[1][2] -
                  dxdr[0][0] * dxdr[2][1] * dxdr[1][2] -
                  dxdr[1][0] * dxdr[0][1] * dxdr[2][2] -
                  dxdr[2][0] * dxdr[1][1] * dxdr[0][2];
    drdx[0][0] = (dxdr[1][1] * dxdr[2][2] - dxdr[1][2] * dxdr[2][1]) / detj;
    drdx[0][1] = (dxdr[0][2] * dxdr[2][1] - dxdr[0][1] * dxdr[2][2]) / detj;
    drdx[0][2] = (dxdr[0][1] * dxdr[1][2] - dxdr[0][2] * dxdr[1][1]) / detj;
    drdx[1][0] = (dxdr[1][2] * dxdr[2][0] - dxdr[1][0] * dxdr[2][2]) / detj;
    drdx[1][1] = (dxdr[0][0] * dxdr[2][2] - dxdr[0][2] * dxdr[2][0]) / detj;
    drdx[1][2] = (dxdr[0][2] * dxdr[1][0] - dxdr[0][0] * dxdr[1][2]) / detj;
    drdx[2][0] = (dxdr[1][0] * dxdr[2][1] - dxdr[1][1] * dxdr[2][0]) / detj;
    drdx[2][1] = (dxdr[0][1] * dxdr[2][0] - dxdr[0][0] * dxdr[2][1]) / detj;
    drdx[2][2] = (dxdr[0][0] * dxdr[1][1] - dxdr[0][1] * dxdr[1][0]) / detj;
    return drdx;
}

void ErrorEstimator::posterior_error_estimation(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    std::vector<std::vector<double>> &displacement,
    const std::vector<std::vector<double>> &material,
    const std::vector<int> &matid_arr,
    const std::map<std::set<int>, std::vector<int>> &face_to_elems,
    const int nload, const std::string &data_dir,
    std::vector<double> &eta_arr) {
    eta_arr.resize(cny.size());

    for (int iload = 0; iload < nload; iload++) {
        std::cout << "iload: " << iload << std::endl;
        read_displacement(data_dir, coor.size(), displacement, iload);
        std::vector<double> rk_arr;
        std::map<std::set<int>, double> re_map;
        calc_rk(cny, coor, displacement, material, matid_arr, rk_arr);
        calc_re(cny, coor, displacement, material, matid_arr, face_to_elems,
                re_map);

        for (int ielem = 0; ielem < cny.size(); ielem++) {
            double eta = 0.0;

            std::vector<int> elem = cny[ielem];
            std::vector<std::vector<double>> xnode(10, std::vector<double>(3));
            for (int inode = 0; inode < 10; inode++) {
                for (int idim = 0; idim < 3; idim++) {
                    xnode[inode][idim] = coor[elem[inode]][idim];
                }
            }

            double rk = rk_arr[ielem];
            double hk = calc_hk(xnode);

            eta += pow(hk, 2) * rk;

            for (int iface = 0; iface < 4; iface++) {
                std::set<int> face;
                for (int inode = 0; inode < 3; inode++) {
                    face.insert(elem[(iface + inode) % 4]);
                }
                double beta = 1.0 / face_to_elems.at(face).size();
                double re = re_map.at(face);
                double he = calc_he(coor, face);

                eta += he * re;
            }

            eta = sqrt(eta);
            eta_arr[ielem] = std::max(eta_arr[ielem], eta);
        }
    }
}
