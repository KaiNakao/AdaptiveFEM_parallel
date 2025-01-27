#include <mpi.h>

#include <iostream>
#include <vector>
#include <cmath>

#include "reader.hpp"
#include "sub_f.hpp"
#include "find_nodes.hpp"
#include "writer.hpp"

int main(int argc, char **argv) {
    int myid, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int nobs;
    std::vector<double> xnode_obs_buf(3 * 8);
    std::vector<std::vector<double>> xnode_obs(8, std::vector<double>(3));
    if (myid == 0) {
        read_nobs_xnode_(nobs, xnode_obs_buf[0]);
    }
    MPI_Bcast(&nobs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xnode_obs_buf[0], 3 * 8, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int inode = 0; inode < 8; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            xnode_obs[inode][idim] = xnode_obs_buf[inode * 3 + idim];
        }
    }

    std::vector<double> xc(3), mvec(6);
    if (myid == 0) {
        read_centroid_(xc[0], mvec[0]);
    }
    MPI_Bcast(&xc[0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mvec[0], 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::vector<double> fvec(24);
    calc_fvec_(fvec[0], xnode_obs_buf[0], xc[0], mvec[0]);

    int nelem, nnode, nmaterial, nload;
    std::vector<std::vector<int>> cny;
    std::vector<int> matid_arr;
    std::vector<std::vector<double>> coor;
    read_shape(myid, nelem, nnode, nmaterial, nload);
    read_mesh(myid, nelem, nnode, cny, coor, matid_arr);

    double ds;
    if (myid == 0) {
        read_ds(ds);
        std::cout << "ds = " << ds << std::endl;
    }
    MPI_Bcast(&ds, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::vector<int> node_id_obs(8);
    find_nodes_on_obs_points(xnode_obs, cny, coor, nnode, nelem, ds, node_id_obs);

    write_nodal_force(myid, node_id_obs, fvec);

    MPI_Finalize();
}