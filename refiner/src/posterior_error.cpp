#include <cmath>
#include <algorithm>
#include <iostream>

#include "posterior_error.hpp"

void poterior_error_estimation(const std::vector<std::vector<int>> &cny, 
                               const std::vector<std::vector<double>> &coor, 
                               const std::vector<std::vector<double>> &displacement,
                               const std::vector<std::vector<int>> &adj_elem) {
    // 実際の物性への対応は未実装
    double lam = 1.0, mu = 1.0;
    std::vector<std::vector<double>> dmat(6, std::vector<double>(6));
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            dmat[i][j] = 0.0;
        }
    }
    dmat[0][0] = 2.0 * mu + lam; dmat[0][1] = lam; dmat[0][2] = lam;
    dmat[1][0] = lam; dmat[1][1] = 2.0 * mu + lam; dmat[1][2] = lam;
    dmat[2][0] = lam; dmat[2][1] = lam; dmat[2][2] = 2.0 * mu + lam;
    dmat[3][3] = mu; dmat[4][4] = mu; dmat[5][5] = mu;

    // dndrdr[i][j][k] = dn_i/(dr_j dr_k)
    std::vector<std::vector<std::vector<double>>> dndrdr = calc_dndrdr();

    for (int ielem = 1; ielem < cny.size(); ielem++) {
        double eta = 0.0;
        std::vector<int> elem_id = cny[ielem];
        std::vector<std::vector<double>> xnode(10, std::vector<double>(3));
        std::vector<double> uvec(30);
        for (int inode = 0; inode < 10; inode++) {
            for (int idim = 0; idim < 3; idim++) {
                xnode[inode][idim] = coor[elem_id[inode]][idim];
                uvec[inode * 3 + idim] = displacement[elem_id[inode]][idim];
            }
        }

        // dxdr[i][j] = dx_i/dr_j
        std::vector<std::vector<double>> dxdr = calc_dxdr(xnode);
        // drdx[i][j] = dr_i/dx_j
        std::vector<std::vector<double>> drdx = calc_drdx(dxdr);
        // volume of the element
        double vol = calc_volume(dxdr);
        // dndxdx[i][j][k] = dn_i/(dx_j dx_k)
        std::vector<std::vector<std::vector<double>>> dndxdx = calc_dndxdx(dndrdr, drdx);

        // length of the longest edge
        double hk = calc_hk(xnode);

        // internal residual
        std::vector<double> rk(3);
        for (int idim = 0; idim < 3; idim++) {
            rk[idim] = 0.0;
        }

        for (int idim = 0; idim < 3; idim++) {
            std::vector<std::vector<double>> bmat_grad(6, std::vector<double>(30));
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
                rk[0] += sigma_grad[0]; // sigma_xx
                rk[1] += sigma_grad[3]; // sigma_xy
                rk[2] += sigma_grad[5]; // sigma_xz
            }
            else if (idim == 1) {
                rk[0] += sigma_grad[3]; // sigma_xy
                rk[1] += sigma_grad[1]; // sigma_yy
                rk[2] += sigma_grad[4]; // sigma_yz
            }
            else if (idim == 2) {
                rk[0] += sigma_grad[5]; // sigma_xz
                rk[1] += sigma_grad[4]; // sigma_yz
                rk[2] += sigma_grad[2]; // sigma_zz
            }
        }
        // 断層ではmoment tensorによる力をrkに加える

        // 定数を要素内積分
        for (int idim = 0; idim < 3; idim++) {
            rk[idim] *= vol;
        }

        eta += pow(hk, 0) * (pow(rk[0], 2) + pow(rk[1], 2) + pow(rk[2], 2));
        std::cout << "element: " << ielem << " eta: " << eta << std::endl;
    }
}

double calc_hk(const std::vector<std::vector<double>> &xnode) {

    // length of each edge
    std::vector<double> edge_length(6);
    for (int iedge = 0; iedge < 6; iedge++) {
        int inode1, inode2;
        if (iedge == 0) {
            inode1 = 0; inode2 = 1;
        }
        else if (iedge == 1) {
            inode1 = 0; inode2 = 2;
        }
        else if (iedge == 2) {
            inode1 = 0; inode2 = 3;
        }
        else if (iedge == 3) {
            inode1 = 1; inode2 = 2;
        }
        else if (iedge == 4) {
            inode1 = 1; inode2 = 3;
        }
        else if (iedge == 5) {
            inode1 = 2; inode2 = 3;
        }
        for (int idim = 0; idim < 3; idim++) {
            edge_length[iedge] += pow(xnode[inode1][idim] - xnode[inode2][idim], 2);
        }
        edge_length[iedge] = sqrt(edge_length[iedge]);
    }
    // longest edge
    double hk = *std::max_element(edge_length.begin(), edge_length.end());
    return hk;
}

std::vector<double> calc_normal_vec(const std::vector<std::vector<double>> &xnode) {
    std::vector<double> normal_vec(3);
    std::vector<double> vec1(3), vec2(3);
    for (int idim = 0; idim < 3; idim++) {
        vec1[idim] = xnode[1][idim] - xnode[0][idim];
        vec2[idim] = xnode[2][idim] - xnode[0][idim];
    }
    normal_vec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    normal_vec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    normal_vec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
    return normal_vec;
}

std::vector<std::vector<std::vector<double>>> calc_dndrdr() {
    std::vector<std::vector<std::vector<double>>> 
    dndrdr(10, std::vector<std::vector<double>>(3, std::vector<double>(3)));

    dndrdr[0][0][0] = 4.0; dndrdr[0][0][1] = 4.0; dndrdr[0][0][2] = 4.0;
    dndrdr[0][1][0] = 4.0; dndrdr[0][1][1] = 4.0; dndrdr[0][1][2] = 4.0;
    dndrdr[0][2][0] = 4.0; dndrdr[0][2][1] = 4.0; dndrdr[0][2][2] = 4.0;

    dndrdr[1][0][0] = 4.0; dndrdr[1][0][1] = 0.0; dndrdr[1][0][2] = 0.0;
    dndrdr[1][1][0] = 0.0; dndrdr[1][1][1] = 0.0; dndrdr[1][1][2] = 0.0;
    dndrdr[1][2][0] = 0.0; dndrdr[1][2][1] = 0.0; dndrdr[1][2][2] = 0.0;

    dndrdr[2][0][0] = 0.0; dndrdr[2][0][1] = 0.0; dndrdr[2][0][2] = 0.0;
    dndrdr[2][1][0] = 0.0; dndrdr[2][1][1] = 4.0; dndrdr[2][1][2] = 0.0;
    dndrdr[2][2][0] = 0.0; dndrdr[2][2][1] = 0.0; dndrdr[2][2][2] = 0.0;

    dndrdr[3][0][0] = 0.0; dndrdr[3][0][1] = 0.0; dndrdr[3][0][2] = 0.0;
    dndrdr[3][1][0] = 0.0; dndrdr[3][1][1] = 0.0; dndrdr[3][1][2] = 0.0;
    dndrdr[3][2][0] = 0.0; dndrdr[3][2][1] = 0.0; dndrdr[3][2][2] = 4.0;

    dndrdr[4][0][0] = -8.0; dndrdr[4][0][1] = -4.0; dndrdr[4][0][2] = -4.0;
    dndrdr[4][1][0] = -4.0; dndrdr[4][1][1] = 0.0; dndrdr[4][1][2] = 4.0;
    dndrdr[4][2][0] = -4.0; dndrdr[4][2][1] = 0.0; dndrdr[4][2][2] = 0.0;

    dndrdr[5][0][0] = 0.0; dndrdr[5][0][1] = 4.0; dndrdr[5][0][2] = 0.0;
    dndrdr[5][1][0] = 4.0; dndrdr[5][1][1] = 0.0; dndrdr[5][1][2] = 0.0;
    dndrdr[5][2][0] = 0.0; dndrdr[5][2][1] = 0.0; dndrdr[5][2][2] = 0.0;

    dndrdr[6][0][0] = 0.0; dndrdr[6][0][1] = -4.0; dndrdr[6][0][2] = 0.0;
    dndrdr[6][1][0] = -4.0; dndrdr[6][1][1] = -8.0; dndrdr[6][1][2] = -4.0;
    dndrdr[6][2][0] = 0.0; dndrdr[6][2][1] = -4.0; dndrdr[6][2][2] = 0.0;

    dndrdr[7][0][0] = 0.0; dndrdr[7][0][1] = 0.0; dndrdr[7][0][2] = -4.0;
    dndrdr[7][1][0] = 0.0; dndrdr[7][1][1] = 0.0; dndrdr[7][1][2] = -4.0;
    dndrdr[7][2][0] = -4.0; dndrdr[7][2][1] = -4.0; dndrdr[7][2][2] = -8.0;

    dndrdr[8][0][0] = 0.0; dndrdr[8][0][1] = 0.0; dndrdr[8][0][2] = 4.0;
    dndrdr[8][1][0] = 0.0; dndrdr[8][1][1] = 0.0; dndrdr[8][1][2] = 0.0;
    dndrdr[8][2][0] = 4.0; dndrdr[8][2][1] = 0.0; dndrdr[8][2][2] = 0.0;

    dndrdr[9][0][0] = 0.0; dndrdr[9][0][1] = 0.0; dndrdr[9][0][2] = 0.0;
    dndrdr[9][1][0] = 0.0; dndrdr[9][1][1] = 0.0; dndrdr[9][1][2] = 4.0;
    dndrdr[9][2][0] = 0.0; dndrdr[9][2][1] = 4.0; dndrdr[9][2][2] = 0.0;

    return dndrdr;
}

std::vector<std::vector<double>> calc_dxdr(const std::vector<std::vector<double>> &xnode) {
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

std::vector<std::vector<std::vector<double>>> calc_dndxdx(
    const std::vector<std::vector<std::vector<double>>> &dndrdr, 
    const std::vector<std::vector<double>> &drdx) {
    std::vector<std::vector<std::vector<double>>> dndxdx(10, std::vector<std::vector<double>>(3, std::vector<double>(3)));
    for (int inode = 0; inode < 10; inode++) {
        for (int ix = 0; ix < 3; ix++) {
            for (int jx = 0; jx < 3; jx++) {
                double tmp = 0.0;
                for (int ir = 0; ir < 3; ir++) {
                    for (int jr = 0; jr < 3; jr++) {
                        tmp += dndrdr[inode][ir][jr] * drdx[ir][ix] * drdx[jr][jx];
                    }
                }
                dndxdx[inode][ix][jx] = tmp;
            }
        }
    }
    return dndxdx;
}

double calc_volume(std::vector<std::vector<double>> &dxdr) {
    double detj = dxdr[0][0] * dxdr[1][1] * dxdr[2][2] + dxdr[1][0] * dxdr[2][1] * dxdr[0][2] + dxdr[2][0] * dxdr[0][1] * dxdr[1][2]
               - dxdr[0][0] * dxdr[2][1] * dxdr[1][2] - dxdr[1][0] * dxdr[0][1] * dxdr[2][2] - dxdr[2][0] * dxdr[1][1] * dxdr[0][2];
    double vol = detj / 6.0;
    return vol;
}

std::vector<std::vector<double>> calc_drdx(
    const std::vector<std::vector<double>> &dxdr) {
    std::vector<std::vector<double>> drdx(3, std::vector<double>(3));
    double detj = dxdr[0][0] * dxdr[1][1] * dxdr[2][2] + dxdr[1][0] * dxdr[2][1] * dxdr[0][2] + dxdr[2][0] * dxdr[0][1] * dxdr[1][2]
               - dxdr[0][0] * dxdr[2][1] * dxdr[1][2] - dxdr[1][0] * dxdr[0][1] * dxdr[2][2] - dxdr[2][0] * dxdr[1][1] * dxdr[0][2];
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