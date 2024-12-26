#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "mpi.h"
extern "C" {
#include "hdf_lib.c"
}

bool is_on_surf(const std::set<int> &face,
                const std::map<int, std::vector<int>> &neighbor_map_arr,
                const std::map<int, std::set<std::set<int>>> &neighbor_face) {
    // std::cout << "key in neighbor_face: ";
    // for (const auto &p : neighbor_face)
    //     std::cout << p.first << " ";
    // std::cout << std::endl;
    for (const auto &p : neighbor_map_arr) {
        const int rank = p.first;
        if (neighbor_face.find(rank) == neighbor_face.end()) {
            continue;
        }
        const std::vector<int> &nodes = p.second;
        const std::set<std::set<int>> &neighbor_face_tmp =
            neighbor_face.at(rank);
        int cnt = 0;
        int pos_arr[3];
        for (int node_id : face) {
            auto pos = std::find(nodes.begin(), nodes.end(), node_id);
            if (pos != nodes.end()) {
                pos_arr[cnt] = std::distance(nodes.begin(), pos);
                cnt++;
            }
        }
        if (cnt < 3) continue;
        std::set<int> face_idx{pos_arr[0], pos_arr[1], pos_arr[2]};
        std::cout << "face_idx: ";
        for (int idx : face_idx) {
            std::cout << idx << " ";
        }

        std::cout << "(rank: " << rank << ")" << std::endl;
        if (neighbor_face_tmp.find(face_idx) != neighbor_face_tmp.end()) {
            return false;
        }
    }
    return true;
}

bool is_on_domain_boundary(const std::set<int> &face,
                           const std::vector<std::vector<double>> &coor,
                           const double &xmin, const double &xmax,
                           const double &ymin, const double &ymax,
                           const double &zmin) {
    std::vector<double> gpoint(3);
    for (int idim = 0; idim < 3; idim++) {
        gpoint[idim] = 0.0;
        for (int node_id : face) {
            gpoint[idim] += coor[node_id][idim];
        }
        gpoint[idim] /= 3.0;
    }
    if (gpoint[0] == xmin || gpoint[0] == xmax || gpoint[1] == ymin ||
        gpoint[1] == ymax || gpoint[2] == zmin) {
        return true;
    }
    return false;
}

int main(int argc, char **argv) {
    int myid, numprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    long fid = 0, count, size, start;
    // read mdata
    std::ostringstream ss;
    ss << std::setw(6) << std::setfill('0') << myid;
    std::string filename = "mdata/" + ss.str() + ".data.h5";
    std::cout << "filename: " << filename << std::endl;

    hdf_open_file_(filename.c_str(), &fid);
    size = 3;
    int settings[size];
    hdf_read_int_array_(&fid, "/setting", settings, &size, &count);
    int nnode = settings[0];
    int nelem = settings[1];
    std::cout << "nelem: " << nelem << std::endl;

    size = 11 * nelem;
    std::vector<std::vector<int>> cny(nelem, std::vector<int>(10));
    std::vector<int> material(nelem);
    int conn_buf[size];
    hdf_read_int_array_(&fid, "/conn", conn_buf, &size, &count);
    for (int ielem = 0; ielem < nelem; ielem++) {
        for (int inode = 0; inode < 10; inode++) {
            cny[ielem][inode] = conn_buf[ielem * 11 + inode] - 1;
        }
        material[ielem] = conn_buf[ielem * 11 + 10] - 1;
    }

    size = 3 * nnode;
    std::vector<std::vector<double>> coor(nnode, std::vector<double>(3));
    double coor_buf[size];
    hdf_read_double_array_(&fid, "/coor", coor_buf, &size, &count);
    for (int inode = 0; inode < nnode; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            coor[inode][idim] = coor_buf[inode * 3 + idim];
        }
    }

    start = 0;
    size = 1;
    std::vector<int> mpinode_buf;
    mpinode_buf.resize(size);
    hdf_read_int_array_part_(&fid, "/MPInode", &start, &size,
                             &(mpinode_buf[0]));
    int nneighbor = mpinode_buf[0];
    std::vector<int> neighbor_size(nneighbor);
    std::vector<int> neighbor_rank(nneighbor);
    start += size;
    size = 2 * nneighbor;
    mpinode_buf.resize(size);
    hdf_read_int_array_part_(&fid, "/MPInode", &start, &size,
                             &(mpinode_buf[0]));
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        neighbor_rank[ineighbor] = mpinode_buf[2 * ineighbor];
        neighbor_size[ineighbor] = mpinode_buf[2 * ineighbor + 1];
    }
    start += size;
    size = std::accumulate(neighbor_size.begin(), neighbor_size.end(), 0);
    mpinode_buf.resize(size);
    hdf_read_int_array_part_(&fid, "/MPInode", &start, &size,
                             &(mpinode_buf[0]));
    std::map<int, std::set<int>> neighbor_map_map;
    std::map<int, std::vector<int>> neighbor_map_arr;
    int pt = 0;
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        for (int i = 0; i < neighbor_size[ineighbor]; i++) {
            neighbor_map_arr[neighbor_rank[ineighbor]].push_back(
                mpinode_buf[pt + i] - 1);
        }
        pt += neighbor_size[ineighbor];
    }
    for (const auto &p : neighbor_map_arr) {
        neighbor_map_map[p.first] =
            std::set<int>(p.second.begin(), p.second.end());
    }
    hdf_close_file_(&fid);

    std::cout << "number of neighbor: " << nneighbor << std::endl;
    std::cout << "neighbor rank: " << std::endl;
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        int rank = neighbor_rank[ineighbor];
        std::cout << rank << " ";
    }
    std::cout << std::endl;

    // read model domain
    double xmin, xmax, ymin, ymax, zmin;
    if (myid == 0) {
        filename = "data/modeldomain.dat";
        std::ifstream ifs(filename);
        std::string buf;
        getline(ifs, buf);  // xmin-plane
        getline(ifs, buf);
        xmin = std::stod(buf);
        getline(ifs, buf);  // xmax-plane
        getline(ifs, buf);
        xmax = std::stod(buf);
        getline(ifs, buf);  // ymin-plane
        getline(ifs, buf);
        ymin = std::stod(buf);
        getline(ifs, buf);  // ymax-plane
        getline(ifs, buf);
        ymax = std::stod(buf);
        getline(ifs, buf);  // zmin-plane
        getline(ifs, buf);
        zmin = std::stod(buf);
    }
    MPI_Bcast(&xmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ymin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ymax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&zmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::map<std::set<int>, std::vector<int>> face_to_elem;
    for (int ielem = 0; ielem < nelem; ielem++) {
        std::vector<int> nodes = cny[ielem];
        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(nodes[(iface + inode) % 4]);
            }

            face_to_elem[face].push_back(ielem);
        }
    }

    std::map<int, std::set<std::set<int>>> boundary_face;
    for (const auto &p : face_to_elem) {
        const std::set<int> &face = p.first;
        for (const auto &q : neighbor_map_arr) {
            const int rank = q.first;
            const std::vector<int> &nodes = q.second;
            int cnt = 0;
            int pos_arr[3];
            for (int node_id : face) {
                auto pos = std::find(nodes.begin(), nodes.end(), node_id);
                if (pos != nodes.end()) {
                    pos_arr[cnt] = std::distance(nodes.begin(), pos);
                    cnt++;
                }
            }
            if (cnt == 3) {
                boundary_face[rank].insert(
                    std::set<int>{pos_arr[0], pos_arr[1], pos_arr[2]});
            }
        }
    }

    // if (myid == 10) {
    //     std::cout << "boundary_face: " << std::endl;
    //     neighbor_map[10]
    // }

    MPI_Request request_arr[nneighbor][2];
    MPI_Status status_arr[nneighbor][2];
    std::map<int, std::set<std::set<int>>> neighbor_face;
    std::vector<int> nface_send(nneighbor);
    std::vector<int> nface_recv(nneighbor);
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        int rank = neighbor_rank[ineighbor];
        nface_send[ineighbor] = boundary_face[rank].size();
    }
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        int rank = neighbor_rank[ineighbor];
        MPI_Isend(&nface_send[ineighbor], 1, MPI_INT, rank, 0, MPI_COMM_WORLD,
                  &request_arr[ineighbor][0]);
        MPI_Irecv(&nface_recv[ineighbor], 1, MPI_INT, rank, 0, MPI_COMM_WORLD,
                  &request_arr[ineighbor][1]);
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);

    // for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
    //     int rank = neighbor_rank[ineighbor];
    //     std::cout << myid << " -> " << rank << ": " << nface_send[ineighbor]
    //               << std::endl;
    //     std::cout << myid << " <- " << rank << ": " << nface_recv[ineighbor]
    //               << std::endl;
    // }

    std::vector<std::vector<int>> buf_send(nneighbor);
    std::vector<std::vector<int>> buf_recv(nneighbor);
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        int rank = neighbor_rank[ineighbor];
        buf_send[ineighbor].resize(3 * nface_send[ineighbor]);
        int cnt = 0;
        for (const auto &face : boundary_face[rank]) {
            for (int node_id : face) {
                buf_send[ineighbor][cnt] = node_id;
                cnt++;
            }
        }
        buf_recv[ineighbor].resize(3 * nface_recv[ineighbor]);
    }
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        int rank = neighbor_rank[ineighbor];
        MPI_Isend(&buf_send[ineighbor][0], buf_send[ineighbor].size(), MPI_INT,
                  rank, 0, MPI_COMM_WORLD, &request_arr[ineighbor][0]);
        MPI_Irecv(&buf_recv[ineighbor][0], buf_recv[ineighbor].size(), MPI_INT,
                  rank, 0, MPI_COMM_WORLD, &request_arr[ineighbor][1]);
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        int rank = neighbor_rank[ineighbor];
        for (int iface = 0; iface < nface_recv[ineighbor]; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(buf_recv[ineighbor][3 * iface + inode]);
            }
            neighbor_face[rank].insert(face);
        }
    }

    std::cout << "comunication done" << std::endl;

    // for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
    //     int rank = neighbor_rank[ineighbor];
    //     std::cout << myid << " <- " << rank << ": "
    //               << neighbor_face[rank].size() << std::endl;
    //     for (const auto &face : neighbor_face[rank]) {
    //         for (int node_id : face) {
    //             std::cout << node_id << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    std::set<int> surf_nodes;
    for (auto p : face_to_elem) {
        std::set<int> face = p.first;
        // if (face == std::set<int>{1011, 1103, 1004}) {
        if (face == std::set<int>{905, 911, 1005}) {
            std::cout << "face: ";
            for (int node_id : face) {
                std::cout << node_id << " ";
            }
            std::cout << std::endl;
            std::cout << "nelem for face: " << p.second.size() << std::endl;
            for (int elem_id : p.second) {
                std::cout << elem_id << " ";
            }
            std::cout << std::endl;
            std::cout << "is_on_domain_boundary: "
                      << is_on_domain_boundary(face, coor, xmin, xmax, ymin,
                                               ymax, zmin)
                      << std::endl;
            std::cout << "is_on_surf: "
                      << is_on_surf(face, neighbor_map_arr, neighbor_face)
                      << std::endl;
        }
        std::vector<int> elems = p.second;
        if (elems.size() == 1) {
            if (is_on_domain_boundary(face, coor, xmin, xmax, ymin, ymax,
                                      zmin)) {
                continue;
            }
            // if (face == std::set<int>{2189, 1992, 2175}) {
            //     std::cout << "face: ";
            //     for (int node_id : face) {
            //         std::cout << node_id << " ";
            //     }
            //     std::cout << std::endl;
            //     std::cout << "here" << std::endl;
            // }

            if (is_on_surf(face, neighbor_map_arr, neighbor_face)) {
                for (int node_id : face) {
                    // if (node_id == 1005) {
                    //     std::cout << "node 1005" << std::endl;
                    //     for (int node_id : face) {
                    //         std::cout << node_id << " ";
                    //     }
                    //     std::cout << std::endl;
                    // }
                    surf_nodes.insert(node_id);
                }
            }
        }
    }

    // std::cout << "surf_nodes: " << surf_nodes.size() << std::endl;
    // for (int inode : surf_nodes) {
    //     std::cout << inode << std::endl;
    //     std::cout << coor[inode][0] << " " << coor[inode][1] << " "
    //               << coor[inode][2] << std::endl;
    // }
    // std::cout << std::endl;

    std::cout << "1011: " << surf_nodes.count(1011) << std::endl;
    std::cout << "1103: " << surf_nodes.count(1103) << std::endl;
    std::cout << "1004: " << surf_nodes.count(1004) << std::endl;

    // connectivity of surface nodes
    std::vector<std::vector<int>> surf_cny;
    for (int ielem = 0; ielem < nelem; ielem++) {
        std::vector<int> nodes;
        std::vector<std::vector<double>> xnode;
        int cnt = 0;
        for (int inode = 0; inode < 4; inode++) {
            if (surf_nodes.find(cny[ielem][inode]) != surf_nodes.end()) {
                nodes.push_back(cny[ielem][inode]);
                xnode.push_back(coor[cny[ielem][inode]]);
                cnt++;
            }
        }

        if (ielem == 4181) {
            std::cout << "elem 4181" << std::endl;
            std::cout << "cnt: " << cnt << std::endl;
            std::cout << "nodes: ";
            for (int node_id : cny[ielem]) {
                std::cout << node_id << " ";
            }
            std::cout << std::endl;
            std::cout << "xnode: ";
            for (int node_id : cny[ielem]) {
                for (int idim = 0; idim < 3; idim++) {
                    std::cout << coor[node_id][idim] << " ";
                }
                std::cout << std::endl;
            }
        }

        // surface element
        if (cnt != 3) continue;
        std::vector<double> v1(3), v2(3), cross(3);
        for (int idim = 0; idim < 3; idim++) {
            v1[idim] = xnode[1][idim] - xnode[0][idim];
            v2[idim] = xnode[2][idim] - xnode[0][idim];
        }
        cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
        cross[1] = v1[2] * v2[0] - v1[0] * v2[2];
        cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
        if (cross[2] < 0) {
            std::cout << "cross: ";
            for (int idim = 0; idim < 3; idim++) {
                std::cout << cross[idim] << " ";
            }
            std::swap(nodes[1], nodes[2]);
            for (int idim = 0; idim < 3; idim++) {
                std::swap(xnode[1][idim], xnode[2][idim]);
            }
        }

        nodes.push_back(ielem);
        surf_cny.push_back(nodes);
    }

    filename = "surf_mesh/" + ss.str() + "_nelem.dat";
    std::ofstream ofs(filename);
    ofs << surf_cny.size() << std::endl;
    ofs.close();

    filename = "surf_mesh/" + ss.str() + "_cny.bin";
    std::ofstream ofs_bin(filename, std::ios::binary);
    for (int ielem = 0; ielem < surf_cny.size(); ielem++) {
        for (int inode = 0; inode < 3; inode++) {
            surf_cny[ielem][inode] += 1;
            ofs_bin.write((char *)&surf_cny[ielem][inode], sizeof(int));
        }
        surf_cny[ielem][3] += 1;
        ofs_bin.write((char *)&surf_cny[ielem][3], sizeof(int));
    }
    ofs_bin.close();

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
}
