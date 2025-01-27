#include "find_nodes.hpp"

void find_nodes_on_obs_points(const std::vector<std::vector<double>> &xnode_obs,
                              const std::vector<std::vector<int>> &cny,
                              const std::vector<std::vector<double>> &coor, 
                              const int &nnode, const int &nelem, const double &ds,
                              std::vector<int> &node_id_obs) {


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

    std::vector<std::vector<double>> xnode(10, std::vector<double>(3));
    for (int ipoint = 0; ipoint < 8; ipoint++) {
        // std::cout << "obs point " << ipoint << std::endl;
        int found_id = -1;
        bool found_local = false, found_global = false;
        double xobs = xnode_obs[ipoint][0];
        double yobs = xnode_obs[ipoint][1];
        double zobs = xnode_obs[ipoint][2];
        for (int ielem = 0; ielem < nelem; ielem++) {
            double xmin_elem, xmax_elem, ymin_elem, ymax_elem, zmin_elem,
                zmax_elem;
            for (int inode = 0; inode < 4; inode++) {
                int node_id = cny[ielem][inode];
                xnode[inode] = coor[node_id];
            }
            
            for (int inode = 0; inode < 4; inode++) {
                double dist = sqrt(pow(xnode[inode][0] - xobs, 2) +
                                   pow(xnode[inode][1] - yobs, 2) +
                                   pow(xnode[inode][2] - zobs, 2));
                if (dist < ds*1e-6) {
                    std::cout << "obs point " << ipoint << " is on node "
                                << cny[ielem][inode] << std::endl;
                    found_id = cny[ielem][inode];
                    found_local = true;
                    break;
                }
            }
            if (found_local) {
                break;
            }
        }
        MPI_Allreduce(&found_local, &found_global, 1, MPI_C_BOOL, MPI_LOR,
                      MPI_COMM_WORLD); 
        if (!found_global) {
            std::cout << "obs point " << ipoint << " is not on any node"
                      << std::endl;
            MPI_Finalize();
            std::exit(1);
        }
        if (found_global) {
            node_id_obs[ipoint] = found_id;
        }
        else {
            node_id_obs[ipoint] = -1;
        }
    }

}