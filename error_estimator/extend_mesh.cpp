#include "extend_mesh.hpp"

void extend_mesh(const int &myid,
                 const std::map<int, std::vector<int>> &neighbor_map,
                 std::vector<std::vector<int>> &cny,
                 std::vector<std::vector<double>> &coor,
                 std::vector<int> &matid_arr,
                 std::map<std::set<int>, std::set<int>> &face_to_elems,
                 std::map<int, std::vector<int>> &import_index,
                 std::map<int, std::vector<int>> &export_index) {
    int nelem_org = cny.size();
    int nnode_org = coor.size();
    int nneighbor = neighbor_map.size();
    std::map<int, std::vector<int>> send_buf, recv_buf;
    MPI_Request request_arr[nneighbor][2];
    MPI_Status status_arr[nneighbor][2];

    // send/recv neighbor_map
    for (const auto &p : neighbor_map) {
        const int &rank = p.first;
        const std::vector<int> &nodes = p.second;
        int size = nodes.size();
        send_buf[rank] = nodes;
        recv_buf[rank].resize(size);
    }
    int ineighbor = 0;
    for (const auto &p : neighbor_map) {
        const int &rank = p.first;
        const std::vector<int> &nodes = p.second;
        int size = nodes.size();
        int tag = 0;
        MPI_Isend(&(send_buf[rank][0]), size, MPI_INT, rank, tag,
                  MPI_COMM_WORLD, &request_arr[ineighbor][0]);
        MPI_Irecv(&(recv_buf[rank][0]), size, MPI_INT, rank, tag,
                  MPI_COMM_WORLD, &request_arr[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);

    // create neighbor_node_map
    std::map<int, std::map<int, int>> neighbor_node_map;
    for (const auto &p : neighbor_map) {
        const int &rank = p.first;
        const std::vector<int> &nodes = p.second;
        neighbor_node_map[rank] = std::map<int, int>();
        for (int inode = 0; inode < nodes.size(); inode++) {
            neighbor_node_map.at(rank)[recv_buf.at(rank).at(inode)] =
                nodes.at(inode);
        }
    }

    // create face_to_elems
    face_to_elems.clear();
    for (int ielem = 0; ielem < cny.size(); ielem++) {
        const std::vector<int> &elem = cny[ielem];
        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(elem[(iface + inode) % 4]);
            }
            face_to_elems[face].insert(ielem);
        }
    }

    // find boundary faces
    std::map<int, std::set<std::set<int>>> boundary_face_tmp1;
    // initialize
    for (const auto &p : neighbor_map) {
        const int &rank = p.first;
        boundary_face_tmp1[rank] = std::set<std::set<int>>();
    }
    for (const auto &p : face_to_elems) {
        const std::set<int> &face = p.first;
        const std::set<int> &elems = p.second;
        if (elems.size() == 2) continue;
        for (const auto &q : neighbor_map) {
            const int &rank = q.first;
            const std::vector<int> &nodes = q.second;
            int cnt = 0;
            for (const int &node : face) {
                auto pos = std::find(nodes.begin(), nodes.end(), node);
                if (pos != nodes.end()) cnt++;
            }
            if (cnt == 3) {
                boundary_face_tmp1[rank].insert(face);
                break;
            }
        }
    }

    // send/recv boundary_face_tmp1
    send_buf.clear();
    recv_buf.clear();
    ineighbor = 0;
    if (boundary_face_tmp1.size() != nneighbor) {
        std::cout << "Error: boundary_face_tmp1.size() != nneighbor"
                  << std::endl;
        std::cout << "boundary_face_tmp1.size(): " << boundary_face_tmp1.size()
                  << std::endl;
        std::cout << "nneighbor: " << nneighbor << std::endl;
    }
    for (const auto &p : boundary_face_tmp1) {
        const int &rank = p.first;
        const std::set<std::set<int>> &faces = p.second;
        send_buf[rank] = {(int)faces.size() * 3};
        recv_buf[rank].resize(1);
        int tag = 0;
        MPI_Isend(&(send_buf[rank][0]), 1, MPI_INT, rank, tag, MPI_COMM_WORLD,
                  &request_arr[ineighbor][0]);
        MPI_Irecv(&(recv_buf[rank][0]), 1, MPI_INT, rank, tag, MPI_COMM_WORLD,
                  &request_arr[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "here" << std::endl;

    for (const auto &p : boundary_face_tmp1) {
        const int &rank = p.first;
        const std::set<std::set<int>> &faces = p.second;
        send_buf[rank].resize(faces.size() * 3);
        int cnt = 0;
        for (const std::set<int> &face : faces) {
            for (const int &node_id : face) {
                send_buf[rank][cnt] = node_id;
                cnt++;
            }
        }
        recv_buf[rank].resize(recv_buf[rank][0]);
    }
    ineighbor = 0;
    for (const auto &p : boundary_face_tmp1) {
        const int &rank = p.first;
        int tag = 0;
        MPI_Isend(&(send_buf[rank][0]), send_buf[rank].size(), MPI_INT, rank,
                  tag, MPI_COMM_WORLD, &request_arr[ineighbor][0]);
        MPI_Irecv(&(recv_buf[rank][0]), recv_buf[rank].size(), MPI_INT, rank,
                  tag, MPI_COMM_WORLD, &request_arr[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);

    std::map<int, std::set<std::set<int>>> boundary_face_tmp2;
    for (const auto &p : boundary_face_tmp1) {
        const int &rank = p.first;
        for (int iface = 0; iface < recv_buf[rank].size() / 3; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                int node_id = recv_buf[rank][3 * iface + inode];
                node_id = neighbor_node_map[rank][node_id];
                face.insert(node_id);
            }
            boundary_face_tmp2[rank].insert(face);
        }
    }

    std::map<int, std::set<std::set<int>>> boundary_face;
    for (const auto &p : boundary_face_tmp1) {
        const int &rank = p.first;
        const std::set<std::set<int>> &faces = p.second;
        for (const std::set<int> &face : faces) {
            if (boundary_face_tmp2[rank].count(face) == 1) {
                boundary_face[rank].insert(face);
            }
        }
    }

    // collect boundary nodes/cny
    std::map<int, std::set<int>> boundary_nodes;
    std::map<int, std::set<std::vector<int>>> boundary_cny;
    for (const auto &p : boundary_face) {
        const int &rank = p.first;
        const std::set<std::set<int>> &faces = p.second;
        for (const std::set<int> &face : faces) {
            int elem_id = *face_to_elems[face].begin();
            std::vector<int> nodes = cny[elem_id];
            nodes.push_back(matid_arr[elem_id]);
            for (const int &node_id : nodes) {
                boundary_nodes[rank].insert(node_id);
            }
            boundary_cny[rank].insert(nodes);
        }
    }

    // send/recv boundary_nodes
    send_buf.clear();
    recv_buf.clear();
    std::map<int, std::vector<double>> send_buf_d;
    std::map<int, std::vector<double>> recv_buf_d;
    ineighbor = 0;
    for (const auto &p : boundary_nodes) {
        const int &rank = p.first;
        const std::set<int> &nodes = p.second;
        send_buf[rank] = {(int)nodes.size()};
        recv_buf[rank].resize(1);
        int tag = 0;
        MPI_Isend(&(send_buf[rank][0]), 1, MPI_INT, rank, tag, MPI_COMM_WORLD,
                  &request_arr[ineighbor][0]);
        MPI_Irecv(&(recv_buf[rank][0]), 1, MPI_INT, rank, tag, MPI_COMM_WORLD,
                  &request_arr[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);

    export_index.clear();
    import_index.clear();
    for (const auto &p : boundary_nodes) {
        const int &rank = p.first;
        const std::set<int> &nodes = p.second;
        send_buf[rank].resize(nodes.size());
        send_buf_d[rank].resize(nodes.size() * 3);
        int cnt = 0;
        for (const int &node_id : nodes) {
            send_buf[rank][cnt] = node_id;
            for (int idim = 0; idim < 3; idim++) {
                send_buf_d[rank][3 * cnt + idim] = coor.at(node_id).at(idim);
            }
            cnt++;
        }
        export_index[rank] = send_buf[rank];
        int nrecv = recv_buf[rank][0];
        recv_buf[rank].resize(nrecv);
        recv_buf_d[rank].resize(nrecv * 3);
    }

    ineighbor = 0;
    for (const auto &p : boundary_nodes) {
        const int &rank = p.first;
        int tag = 0;
        MPI_Isend(&(send_buf[rank][0]), send_buf[rank].size(), MPI_INT, rank,
                  tag, MPI_COMM_WORLD, &request_arr[ineighbor][0]);
        MPI_Irecv(&(recv_buf[rank][0]), recv_buf[rank].size(), MPI_INT, rank,
                  tag, MPI_COMM_WORLD, &request_arr[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);

    ineighbor = 0;
    for (const auto &p : boundary_nodes) {
        const int &rank = p.first;
        int tag = 0;
        MPI_Isend(&(send_buf_d[rank][0]), send_buf_d[rank].size(), MPI_DOUBLE,
                  rank, tag, MPI_COMM_WORLD, &request_arr[ineighbor][0]);
        MPI_Irecv(&(recv_buf_d[rank][0]), recv_buf_d[rank].size(), MPI_DOUBLE,
                  rank, tag, MPI_COMM_WORLD, &request_arr[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);

    for (const auto &p : boundary_nodes) {
        const int &rank = p.first;
        for (int inode = 0; inode < recv_buf[rank].size(); inode++) {
            int node_id = recv_buf[rank][inode];
            if (neighbor_node_map[rank].count(node_id) == 0) {
                neighbor_node_map[rank][node_id] = coor.size();
                coor.push_back({recv_buf_d[rank].at(3 * inode),
                                recv_buf_d[rank].at(3 * inode + 1),
                                recv_buf_d[rank].at(3 * inode + 2)});
            }
            import_index[rank].push_back(neighbor_node_map[rank][node_id]);
        }
    }

    // send/recv boundary_cny
    send_buf.clear();
    recv_buf.clear();
    ineighbor = 0;
    for (const auto &p : boundary_cny) {
        const int &rank = p.first;
        const std::set<std::vector<int>> &b_cny = p.second;
        send_buf[rank] = {(int)b_cny.size() * 11};
        recv_buf[rank].resize(1);
        int tag = 0;
        MPI_Isend(&(send_buf[rank][0]), 1, MPI_INT, rank, tag, MPI_COMM_WORLD,
                  &request_arr[ineighbor][0]);
        MPI_Irecv(&(recv_buf[rank][0]), 1, MPI_INT, rank, tag, MPI_COMM_WORLD,
                  &request_arr[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);

    for (const auto &p : boundary_cny) {
        const int &rank = p.first;
        const std::set<std::vector<int>> &b_cny = p.second;
        send_buf[rank].resize(b_cny.size() * 11);
        int cnt = 0;
        for (const std::vector<int> &nodes : b_cny) {
            for (const int &node_id : nodes) {
                send_buf[rank][cnt] = node_id;
                cnt++;
            }
        }
        recv_buf[rank].resize(recv_buf[rank][0]);
    }

    ineighbor = 0;
    for (const auto &p : boundary_cny) {
        const int &rank = p.first;
        int tag = 0;
        MPI_Isend(&(send_buf[rank][0]), send_buf[rank].size(), MPI_INT, rank,
                  tag, MPI_COMM_WORLD, &request_arr[ineighbor][0]);
        MPI_Irecv(&(recv_buf[rank][0]), recv_buf[rank].size(), MPI_INT, rank,
                  tag, MPI_COMM_WORLD, &request_arr[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &request_arr[0][0], &status_arr[0][0]);

    for (const auto &p : boundary_cny) {
        const int &rank = p.first;
        for (int icny = 0; icny < recv_buf[rank].size() / 11; icny++) {
            std::vector<int> nodes(10);
            for (int inode = 0; inode < 10; inode++) {
                int node_id = recv_buf[rank][11 * icny + inode];
                if (neighbor_node_map[rank].count(node_id) == 0) {
                    std::cout << "Error: node_id not found" << std::endl;
                    std::cout << "myid: " << myid << " rank: " << rank
                              << " node_id: " << node_id << std::endl;
                }
                nodes[inode] = neighbor_node_map[rank][node_id];
            }
            int mat_id = recv_buf[rank][11 * icny + 10];
            cny.push_back(nodes);
            matid_arr.push_back(mat_id);
        }
    }

    for (int ielem = nelem_org; ielem < cny.size(); ielem++) {
        const std::vector<int> &elem = cny[ielem];
        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(elem[(iface + inode) % 4]);
            }
            if (face_to_elems.count(face) == 0) continue;
            face_to_elems[face].insert(ielem);
        }
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    // for (int iproc = 0; iproc < 16; iproc++) {
    //     if (myid == iproc) {
    //         for (const auto &p : face_to_elems) {
    //             const std::set<int> &face = p.first;
    //             const std::set<int> &elems = p.second;
    //             if (elems.size() == 2) continue;
    //             std::vector<double> center(3);
    //             for (const int &node_id : face) {
    //                 for (int idim = 0; idim < 3; idim++) {
    //                     center[idim] += coor[node_id][idim];
    //                 }
    //             }
    //             for (int idim = 0; idim < 3; idim++) {
    //                 center[idim] /= 3;
    //             }
    //             if (fabs(center[0]) < 10.0 ||
    //                 fabs(center[0] - 105000.0) < 10.0 ||
    //                 fabs(center[1]) < 10.0 ||
    //                 fabs(center[1] - 110000.0) < 10.0 || fabs(center[2])
    //                 < 10.0) continue;
    //             std::cout << "center: " << center[0] << " " << center[1] << "
    //             "
    //                       << center[2] << std::endl;
    //         }
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }
}
