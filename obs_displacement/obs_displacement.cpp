
#include <iostream>

#include "mpi.h"
#include "reader.hpp"
#include "writer.hpp"

int main(int argc, char *argv[]) {
    int myid, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int nelem, nnode, nmaterial, nload;
    read_shape(myid, nelem, nnode, nmaterial, nload);

    std::vector<std::vector<int>> cny;
    std::vector<int> matid_arr;
    std::vector<std::vector<double>> coor;
    read_mesh(myid, nelem, nnode, cny, coor, matid_arr);

    double xmin = coor[0][0];
    double xmax = coor[0][0];
    double ymin = coor[0][1];
    double ymax = coor[0][1];
    double zmin = coor[0][2];
    double zmax = coor[0][2];
    for (int inode = 1; inode < nnode; inode++) {
        xmin = std::min(xmin, coor[inode][0]);
        xmax = std::max(xmax, coor[inode][0]);
        ymin = std::min(ymin, coor[inode][1]);
        ymax = std::max(ymax, coor[inode][1]);
        zmin = std::min(zmin, coor[inode][2]);
        zmax = std::max(zmax, coor[inode][2]);
    }

    std::vector<std::vector<double>> obs_point;
    read_obs_point(myid, obs_point);

    std::vector<int> found_arr(obs_point.size(), 0);
    std::vector<int> iobs_arr;
    std::vector<int> ielem_arr;
    std::vector<std::vector<double>> r_arr;

    std::vector<double> dx(3), dr(3);
    std::vector<std::vector<double>> dxdr(3, std::vector<double>(3));
    std::vector<std::vector<double>> drdx(3, std::vector<double>(3));
    std::vector<std::vector<double>> xnode(10, std::vector<double>(3));
    for (int iobs = 0; iobs < obs_point.size(); iobs++) {
        double xobs = obs_point[iobs][0];
        double yobs = obs_point[iobs][1];
        double zobs = obs_point[iobs][2];
        if (xobs < xmin || xobs > xmax || yobs < ymin || yobs > ymax ||
            zobs < zmin || zobs > zmax) {
            continue;
        }

        for (int ielem = 0; ielem < nelem; ielem++) {
            double xmin_elem, xmax_elem, ymin_elem, ymax_elem, zmin_elem,
                zmax_elem;
            double detj;
            for (int inode = 0; inode < 10; inode++) {
                int node_id = cny[ielem][inode];
                xnode[inode] = coor[node_id];
            }
            xmin_elem = xnode[0][0];
            xmax_elem = xnode[0][0];
            ymin_elem = xnode[0][1];
            ymax_elem = xnode[0][1];
            zmin_elem = xnode[0][2];
            zmax_elem = xnode[0][2];
            for (int inode = 1; inode < 10; inode++) {
                xmin_elem = std::min(xmin_elem, xnode[inode][0]);
                xmax_elem = std::max(xmax_elem, xnode[inode][0]);
                ymin_elem = std::min(ymin_elem, xnode[inode][1]);
                ymax_elem = std::max(ymax_elem, xnode[inode][1]);
                zmin_elem = std::min(zmin_elem, xnode[inode][2]);
                zmax_elem = std::max(zmax_elem, xnode[inode][2]);
            }
            if (xobs < xmin_elem || xobs > xmax_elem || yobs < ymin_elem ||
                yobs > ymax_elem || zobs < zmin_elem || zobs > zmax_elem) {
                continue;
            }
            for (int ix = 0; ix < 3; ix++) {
                for (int ir = 0; ir < 3; ir++) {
                    dxdr[ix][ir] = xnode[ir + 1][ix] - xnode[0][ix];
                }
            }
            detj = dxdr[0][0] *
                       (dxdr[1][1] * dxdr[2][2] - dxdr[1][2] * dxdr[2][1]) -
                   dxdr[0][1] *
                       (dxdr[1][0] * dxdr[2][2] - dxdr[1][2] * dxdr[2][0]) +
                   dxdr[0][2] *
                       (dxdr[1][0] * dxdr[2][1] - dxdr[1][1] * dxdr[2][0]);
            std::vector<std::vector<double>> drdx(3, std::vector<double>(3));
            drdx[0][0] =
                (dxdr[1][1] * dxdr[2][2] - dxdr[1][2] * dxdr[2][1]) / detj;
            drdx[0][1] =
                (dxdr[0][2] * dxdr[2][1] - dxdr[0][1] * dxdr[2][2]) / detj;
            drdx[0][2] =
                (dxdr[0][1] * dxdr[1][2] - dxdr[0][2] * dxdr[1][1]) / detj;
            drdx[1][0] =
                (dxdr[1][2] * dxdr[2][0] - dxdr[1][0] * dxdr[2][2]) / detj;
            drdx[1][1] =
                (dxdr[0][0] * dxdr[2][2] - dxdr[0][2] * dxdr[2][0]) / detj;
            drdx[1][2] =
                (dxdr[0][2] * dxdr[1][0] - dxdr[0][0] * dxdr[1][2]) / detj;
            drdx[2][0] =
                (dxdr[1][0] * dxdr[2][1] - dxdr[1][1] * dxdr[2][0]) / detj;
            drdx[2][1] =
                (dxdr[0][1] * dxdr[2][0] - dxdr[0][0] * dxdr[2][1]) / detj;
            drdx[2][2] =
                (dxdr[0][0] * dxdr[1][1] - dxdr[0][1] * dxdr[1][0]) / detj;

            dx[0] = xobs - xnode[0][0];
            dx[1] = yobs - xnode[0][1];
            dx[2] = zobs - xnode[0][2];

            double r1 =
                drdx[0][0] * dx[0] + drdx[0][1] * dx[1] + drdx[0][2] * dx[2];
            double r2 =
                drdx[1][0] * dx[0] + drdx[1][1] * dx[1] + drdx[1][2] * dx[2];
            double r3 =
                drdx[2][0] * dx[0] + drdx[2][1] * dx[1] + drdx[2][2] * dx[2];

            if (r1 < -1e-8 || r2 < -1e-8 || r3 < -1e-8 ||
                r1 + r2 + r3 > 1.0 + 1e-8) {
                continue;
            }

            found_arr[iobs] = 1;
            iobs_arr.push_back(iobs);
            ielem_arr.push_back(ielem);
            r_arr.push_back({r1, r2, r3});
            break;
        }
    }

    int nfound = 0;
    for (int iobs = 0; iobs < found_arr.size(); iobs++) {
        nfound += found_arr[iobs];
    }
    std::vector<int> nfound_arr, displs;
    if (myid == 0) {
        nfound_arr.resize(numprocs);
    }
    MPI_Gather(&nfound, 1, MPI_INT, &nfound_arr[0], 1, MPI_INT, 0,
               MPI_COMM_WORLD);

    int nfound_total = 0;
    std::vector<int> iobs_arr_total;
    std::vector<double> u_arr_buf_total;
    if (myid == 0) {
        displs.resize(numprocs);
        displs[0] = 0;
        for (int iproc = 0; iproc < numprocs; iproc++) {
            nfound_total += nfound_arr[iproc];
            if (iproc == 0) {
                continue;
            }
            displs[iproc] = displs[iproc - 1] + nfound_arr[iproc - 1];
        }
        iobs_arr_total.resize(nfound_total);
        u_arr_buf_total.resize(3 * nfound_total);
    }

    MPI_Gatherv(&iobs_arr[0], nfound, MPI_INT, &iobs_arr_total[0],
                &nfound_arr[0], &displs[0], MPI_INT, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        for (int iproc = 0; iproc < numprocs; iproc++) {
            nfound_arr[iproc] *= 3;
            displs[iproc] *= 3;
        }
    }

    std::vector<double> nvec(10), u(3), u_arr_buf(3 * nfound);
    std::vector<std::vector<double>> unode(10, std::vector<double>(3)), result;
    if (myid == 0) {
        result.resize(obs_point.size(), std::vector<double>(3));
    }

    for (int iload = 0; iload < nload; iload++) {
        if (myid == 0) {
            std::cout << "iload: " << iload << std::endl;
        }
        std::vector<std::vector<double>> displacement(coor.size(),
                                                      std::vector<double>(3));
        read_displacement(myid, iload, nnode, displacement);

        for (int ipoint = 0; ipoint < nfound; ipoint++) {
            int iobs = iobs_arr[ipoint];
            int ielem = ielem_arr[ipoint];
            for (int inode = 0; inode < 10; inode++) {
                int node_id = cny[ielem][inode];
                unode[inode] = displacement[node_id];
            }
            double r1 = r_arr[ipoint][0];
            double r2 = r_arr[ipoint][1];
            double r3 = r_arr[ipoint][2];
            double l0 = 1.0 - r1 - r2 - r3;
            double l1 = r1;
            double l2 = r2;
            double l3 = r3;
            nvec[0] = l0 * (2.0 * l0 - 1.0);
            nvec[1] = l1 * (2.0 * l1 - 1.0);
            nvec[2] = l2 * (2.0 * l2 - 1.0);
            nvec[3] = l3 * (2.0 * l3 - 1.0);
            nvec[4] = 4.0 * l0 * l1;
            nvec[5] = 4.0 * l1 * l2;
            nvec[6] = 4.0 * l2 * l0;
            nvec[7] = 4.0 * l3 * l0;
            nvec[8] = 4.0 * l3 * l1;
            nvec[9] = 4.0 * l3 * l2;

            for (int idim = 0; idim < 3; idim++) {
                u[idim] = 0.0;
            }
            for (int inode = 0; inode < 10; inode++) {
                for (int idim = 0; idim < 3; idim++) {
                    u[idim] += unode[inode][idim] * nvec[inode];
                }
            }
            for (int idim = 0; idim < 3; idim++) {
                u_arr_buf[ipoint * 3 + idim] = u[idim];
            }
        }
        MPI_Gatherv(&u_arr_buf[0], 3 * nfound, MPI_DOUBLE, &u_arr_buf_total[0],
                    &nfound_arr[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (myid == 0) {
            for (int ipoint = 0; ipoint < nfound_total; ipoint++) {
                int iobs = iobs_arr_total[ipoint];
                for (int idim = 0; idim < 3; idim++) {
                    result[iobs][idim] = u_arr_buf_total[ipoint * 3 + idim];
                }
            }
            write_displacement_obs(result, iload);
        }
    }
}
