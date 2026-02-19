#include "refinement_scheme.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <tuple>

#include "geo_util.hpp"


Refinement_scheme::Refinement_scheme(
    const int myid, const int numprocs,
    const std::set<int> &elem_refine,
    const std::vector<std::vector<int>> &connectivity,
    const std::vector<std::vector<double>> &coordinates,
    const std::vector<int> &matid_arr,
    const std::map<int, std::vector<int>> &neighbor_nodes)
    : m_elem_refine(elem_refine), m_connectivity(connectivity),
      m_neighbor_nodes(neighbor_nodes){
    for (int ielem = 0; ielem < static_cast<int>(connectivity.size()); ielem++) {
        m_matid_map[ielem] = matid_arr[ielem];
    }
    for (int inode = 0; inode < static_cast<int>(coordinates.size()); inode++) {
        m_coor_map[inode] = coordinates[inode];
    }
    m_myid = myid;
    m_numprocs = numprocs;
}

Refinement_scheme::~Refinement_scheme() {}

void Refinement_scheme::executeRefinement_bisect(
    std::vector<std::vector<int>> &new_conn,
    std::vector<std::vector<double>> &new_coor, std::vector<int> &new_matid_arr,
    std::vector<int> &original,
    std::map<int, std::vector<int>> &new_neighbor_nodes) {

    std::set<Tetrahedron> t, s;
    std::map<Edge, std::set<Tetrahedron>> edge_to_tets;
    std::map<int, std::vector<std::tuple<std::vector<int>, int>>> neighbor_edges; // neighbor_id -> vector of (edge, cut flag)
    initialize_tets(t, s, edge_to_tets);
    MPI_Barrier(MPI_COMM_WORLD);
    build_neighbor_edges(t, neighbor_edges);
    check_neighbor_edges(neighbor_edges, edge_to_tets);
    local_refine(t, s, edge_to_tets, neighbor_edges);
    build_neighbor_nodes(neighbor_edges, new_neighbor_nodes);
    check_neighbor_nodes(new_neighbor_nodes);
    output_new_mesh(t, new_conn, new_coor, new_matid_arr, original);
}

void Refinement_scheme::initialize_tets(std::set<Tetrahedron> &t, 
                                        std::set<Tetrahedron> &s,
                                        std::map<Edge, std::set<Tetrahedron>> &edge_to_tets) {
    // initialize t, s
    for (int ielem = 0; ielem < static_cast<int>(m_connectivity.size()); ielem++) {
        int elem_id = ielem;
        std::vector<int> tet_nodes = m_connectivity[ielem];
        Tetrahedron tet;
        // tet.nodes
        for (int node_id : tet_nodes) {
            tet.nodes.insert(node_id);
        }
        // tet.faces
        for (int iface = 0; iface < 4; iface++) {
            Face face;
            // face.nodes
            std::vector<int> face_nodes;
            for (int inode = 0; inode < 3; inode++) {
                face.nodes.insert(tet_nodes.at((iface + inode + 1) % 4));
                face_nodes.push_back(tet_nodes.at((iface + inode + 1) % 4));
            }
            // face.marked_edge
            std::tuple<double, std::vector<double>, int, int> longest;
            longest = std::make_tuple(0., std::vector<double>(3, 0.), -1, -1);
            for (int iedge = 0; iedge < 3; iedge++) {
                std::vector<int> edge_arr;
                edge_arr.push_back(face_nodes.at(iedge));
                edge_arr.push_back(face_nodes.at((iedge + 1) % 3));
                std::sort(edge_arr.begin(), edge_arr.end());
                int idx1 = edge_arr.at(0);
                int idx2 = edge_arr.at(1);
                std::vector<double> p1 = m_coor_map.at(idx1),
                                    p2 = m_coor_map.at(idx2);
                std::vector<double> edge_center(3, 0.);
                for (int idim = 0; idim < 3; idim++) {
                    edge_center[idim] += p1[idim];
                    edge_center[idim] += p2[idim];
                    edge_center[idim] /= 2.0;
                }
                std::tuple<double, std::vector<double>, int, int> length =
                    std::make_tuple(findLength(p1, p2), edge_center, idx1, idx2);
                if (length > longest) {
                    longest = length;
                }
            }
            face.marked_edge.nodes.insert(std::get<2>(longest)); // mark longest edge on face
            face.marked_edge.nodes.insert(std::get<3>(longest)); // mark longest edge on face
            tet.faces.push_back(face);
        }
        // tet.refinement_edge
        std::tuple<double, std::vector<double>, int, int> longest = std::make_tuple(0., std::vector<double>(3, 0.), -1, -1);
        for (const auto &face : tet.faces) {
            std::vector<int> edge_nodes;
            for (int inode : face.marked_edge.nodes) {
                edge_nodes.push_back(inode);
            }
            std::sort(edge_nodes.begin(), edge_nodes.end());
            int idx1 = edge_nodes.at(0), idx2 = edge_nodes.at(1);
            std::vector<double> p1 = m_coor_map.at(idx1),
                                p2 = m_coor_map.at(idx2);
            std::vector<double> edge_center(3, 0.);
            for (int idim = 0; idim < 3; idim++) {
                edge_center[idim] += p1[idim];
                edge_center[idim] += p2[idim];
                edge_center[idim] /= 2.0;
            }
            std::tuple<double, std::vector<double>, int, int> length = std::make_tuple(findLength(p1, p2), edge_center, idx1, idx2);
            if (length > longest) {
                longest = length;
            }
        }
        tet.refinement_edge.nodes.insert(std::get<2>(longest)); // refinement edge: longest edge in tet
        tet.refinement_edge.nodes.insert(std::get<3>(longest));
        // tet.flag
        tet.flag = false;
        // tet.material
        tet.material = m_matid_map.at(elem_id);
        // tet.original
        tet.original = true;
        t.insert(tet);
        if (m_elem_refine.count(elem_id)) {
            s.insert(tet);
        }
    }

    for (const auto &tet : t) {
        std::vector<Edge> edges;
        collect_tet_edges(tet, edges);
        for (const auto &edge : edges) {
            edge_to_tets[edge].insert(tet);
        }
    }
}

// t: all tets
// s: refinemnt tets
void Refinement_scheme::local_refine(std::set<Tetrahedron> &t,
                                     std::set<Tetrahedron> &s,
                                     std::map<Edge, std::set<Tetrahedron>> &edge_to_tets,
                                     std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges) {
    std::map<Edge, int> edges_cut;
    bisect_tets(t, s, edges_cut, edge_to_tets, neighbor_edges);
    refine_to_conformity(t, s, edges_cut, edge_to_tets, neighbor_edges);
    return;
}

void Refinement_scheme::bisect_tets(std::set<Tetrahedron> &t,
                                    std::set<Tetrahedron> &s,
                                    std::map<Edge, int> &edges_cut,
                                    std::map<Edge, std::set<Tetrahedron>> &edge_to_tets,
                                    std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges) {
    for (const Tetrahedron &tet : s) {
        bisect_tet(t, tet, edges_cut, edge_to_tets);
        // remove
        t.erase(tet);
    }
    s.clear();

    sync_edges_cut(edges_cut, neighbor_edges);
    check_neighbor_edges(neighbor_edges, edge_to_tets);
}

int Refinement_scheme::find_tet_type(const Tetrahedron &tet) {
    int parent_type = -1;
    // find tet type
    int cnt = 0;
    std::set<int> nodes_all;
    for (const auto &face : tet.faces) {
        for (int node_id : face.marked_edge.nodes) {
            if (tet.refinement_edge.nodes.count(node_id)) {
                cnt += 1;
            }
            nodes_all.insert(node_id);
        }
    }
    if (cnt == 6) {
        if (nodes_all.size() == 3) {
            if (tet.flag) {
                parent_type = 0;  // Pftype
            } else {
                parent_type = 1;  // Putype
            }
        } else {
            parent_type = 2;  // Atype
        }
    } else if (cnt == 5) {
        parent_type = 3;  // Mtype
    } else if (cnt == 4) {
        parent_type = 4;  // Otype
    }
    if (parent_type == -1) {
        std::cerr << "ERROR: tet type not found" << std::endl;
        std::cerr << "cnt: " << cnt << std::endl;
        std::cerr << "nodes_all.size(): " << nodes_all.size() << std::endl;
        std::exit(1);
    }
    return parent_type;
}

void Refinement_scheme::bisect_tet(std::set<Tetrahedron> &t,
                                    const Tetrahedron &tet,
                                    std::map<Edge, int> &edges_cut,
                                    std::map<Edge, std::set<Tetrahedron>> &edge_to_tets) {
    int parent_type = find_tet_type(tet);
    std::vector<int> nodes_refine_edge, nodes_opposite_edge;
    for (int node_id : tet.nodes) {
        if (tet.refinement_edge.nodes.count(node_id)) {
            nodes_refine_edge.push_back(node_id);
        } else {
            nodes_opposite_edge.push_back(node_id);
        }
    }
    int anode = nodes_refine_edge.at(0);
    int bnode = nodes_refine_edge.at(1);
    int cnode = nodes_opposite_edge.at(0);
    int dnode = nodes_opposite_edge.at(1);

    // midpoint
    std::vector<double> mid_point(3);
    for (int node_id : tet.refinement_edge.nodes) {
        for (int idim = 0; idim < 3; idim++) {
            mid_point[idim] += m_coor_map[node_id][idim];
        }
    }
    for (int idim = 0; idim < 3; idim++) {
        mid_point[idim] /= 2.0;
    }
    int enode;
    if (edges_cut.count(tet.refinement_edge)) {
        enode = edges_cut.at(tet.refinement_edge);
    } else {
        enode = m_coor_map.size();
        edges_cut[tet.refinement_edge] = enode;
        m_coor_map[enode] = mid_point;
    }
    // ///////////////
    // std::vector<double> tet_center(3, 0.);
    // for (int node_id : tet.nodes) {
    //     for (int idim = 0; idim < 3; idim++) {
    //         tet_center[idim] += m_coor_map[node_id][idim];
    //     }
    // }
    // for (int idim = 0; idim < 3; idim++) {
    //     tet_center[idim] /= 4.0;
    // }
    // std::cout << "tet center: ";
    // for (int idim = 0; idim < 3; idim++) {
    //     std::cout << tet_center[idim] << " ";
    // }
    // std::vector<double> edge_center(3, 0.);
    // for (int node_id : tet.refinement_edge.nodes) {
    //     for (int idim = 0; idim < 3; idim++) {
    //         edge_center[idim] += m_coor_map[node_id][idim];
    //     }
    // }
    // for (int idim = 0; idim < 3; idim++) {
    //     edge_center[idim] /= 2.0;
    // }
    // std::cout << "refinement edge center: ";
    // for (int idim = 0; idim < 3; idim++) {
    //     std::cout << edge_center[idim] << " ";
    // }
    // std::cout << std::endl;
    // ///////////////

    std::vector<Edge> parent_edges;
    collect_tet_edges(tet, parent_edges);
    for (const auto &edge : parent_edges) {
        if (!edge_to_tets[edge].count(tet)) {
            std::cerr << "ERROR: edge_to_tets does not contain tet" << std::endl;
            std::exit(1);
        }
        edge_to_tets[edge].erase(tet);
    }

    // new tetrahedra
    Tetrahedron t1, t2;
    // tet.nodes
    t1.nodes = {anode, enode, cnode, dnode};
    t2.nodes = {enode, bnode, cnode, dnode};
    // tet.faces
    // inherited face
    Face face11, face21, face_org;
    face11.nodes = {anode, dnode, cnode};
    for (const auto &face : tet.faces) {
        if (face.nodes == face11.nodes) {
            face_org = face;
            break;
        }
    }
    face11.marked_edge = face_org.marked_edge;
    t1.faces.push_back(face11);
    t1.refinement_edge = face11.marked_edge;

    face21.nodes = {bnode, cnode, dnode};
    for (const auto &face : tet.faces) {
        if (face.nodes == face21.nodes) {
            face_org = face;
            break;
        }
    }
    face21.marked_edge = face_org.marked_edge;
    t2.faces.push_back(face21);
    t2.refinement_edge = face21.marked_edge;

    // cut face
    Face face12, face13, face22, face23;
    face12.nodes = {anode, enode, cnode};
    face13.nodes = {anode, enode, dnode};
    face12.marked_edge.nodes = {anode, cnode};
    face13.marked_edge.nodes = {anode, dnode};
    t1.faces.push_back(face12);
    t1.faces.push_back(face13);

    face22.nodes = {enode, bnode, cnode};
    face23.nodes = {enode, bnode, dnode};
    face22.marked_edge.nodes = {bnode, cnode};
    face23.marked_edge.nodes = {bnode, dnode};
    t2.faces.push_back(face22);
    t2.faces.push_back(face23);

    // new face
    Face face14, face24;
    face14.nodes = {enode, cnode, dnode};
    // Pftype
    if (parent_type == 0) {
        face14.marked_edge.nodes.insert(enode);
        if (t1.refinement_edge.nodes.count(cnode)) {
            face14.marked_edge.nodes.insert(cnode);
        } else {
            face14.marked_edge.nodes.insert(dnode);
        }
    } else {
        face14.marked_edge.nodes = {cnode, dnode};
    }
    t1.faces.push_back(face14);

    face24.nodes = {enode, cnode, dnode};
    // Pftype
    if (parent_type == 0) {
        face24.marked_edge.nodes.insert(enode);
        if (t2.refinement_edge.nodes.count(cnode)) {
            face24.marked_edge.nodes.insert(cnode);
        } else {
            face24.marked_edge.nodes.insert(dnode);
        }
    } else {
        face24.marked_edge.nodes = {cnode, dnode};
    }
    t2.faces.push_back(face24);

    // tet.flag
    // Putype
    if (parent_type == 1) {
        t1.flag = true;
        t2.flag = true;
    } else {
        t1.flag = false;
        t2.flag = false;
    }

    // tet.material
    t1.material = tet.material;
    t2.material = tet.material;

    // tet.original
    t1.original = false;
    t2.original = false;

    t.insert(t1);
    t.insert(t2);

    std::vector<Edge> child_edges;
    collect_tet_edges(t1, child_edges);
    for (const auto &edge : child_edges) {
        edge_to_tets[edge].insert(t1);
    }
    collect_tet_edges(t2, child_edges);
    for (const auto &edge : child_edges) {
        edge_to_tets[edge].insert(t2);
    }
    // std::cout << "bisected edge (" << anode << ", " << bnode << ") to create enode " << enode << " in process " << m_myid << std::endl;
}

void Refinement_scheme::refine_to_conformity(std::set<Tetrahedron> &t,
                                            std::set<Tetrahedron> &s,
                                            std::map<Edge, int> &edges_cut,
                                            std::map<Edge, std::set<Tetrahedron>> &edge_to_tets,
                                            std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges) {
    // check hanging nodes (local)
    s.clear();
    for (const auto &p : edges_cut) {
        const Edge &edge = p.first;
        if (!edge_to_tets.count(edge)) {
            std::cerr << "ERROR: edge_to_tets does not contain edge" << std::endl;
            std::cerr << "edge nodes: ";
            for (int node_id : edge.nodes) {
                std::cerr << node_id << " ";
            }
            // ////////////////// 
            // std::vector<Edge> edges;
            // for (const auto &tet : t) {
            //     collect_tet_edges(tet, edges);
            //     for (const auto &edge_tmp : edges) {
            //         if (edge_tmp.nodes == edge.nodes) {
            //             std::cerr << "edge found in tets" << std::endl;
            //         }
            //     }
            // }
            // //////////////////
            std::exit(1);
        }
        // //////////////
        // std::cout << "edge: ";
        // for (int node_id : edge.nodes) {
        //     std::cout << node_id << " ";
        // }
        // std::cout << std::endl;
        // std::vector<double> edge_center(3, 0.);
        // for (int node_id : edge.nodes) {
        //     for (int idim = 0; idim < 3; idim++) {
        //         edge_center[idim] += m_coor_map[node_id][idim];
        //     }
        // }
        // for (int idim = 0; idim < 3; idim++) {
        //     edge_center[idim] /= 2.0;
        // }
        // std::cout << "edge center: ";
        // for (int idim = 0; idim < 3; idim++) {
        //     std::cout << edge_center[idim] << " ";
        // }
        // std::cout << std::endl;
        // //////////////
        for (const auto &tet : edge_to_tets.at(edge)) {
            s.insert(tet);

            // std::vector<double> center(3, 0.);
            // for (int node_id : tet.nodes) {
            //     for (int idim = 0; idim < 3; idim++) {
            //         center[idim] += m_coor_map[node_id][idim];
            //     }
            // }
            // for (int idim = 0; idim < 3; idim++) {
            //     center[idim] /= 4.0;
            // }
            // std::cout << "center: ";
            // for (int idim = 0; idim < 3; idim++) {
            //     std::cout << center[idim] << " ";
            // }
            // std::cout << std::endl;
        }
    }
    int n_refine_local, n_refine_global;
    n_refine_local = s.size();
    MPI_Allreduce(&n_refine_local, &n_refine_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (m_myid == 0) std::cout << "number of elements with hanging nodes: " << n_refine_global << std::endl;

    // tmp_output(s);
    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Finalize();
    // std::exit(0);


    if (n_refine_global == 0) {
        return;
    } else {
        bisect_tets(t, s, edges_cut, edge_to_tets, neighbor_edges);
        refine_to_conformity(t, s, edges_cut, edge_to_tets, neighbor_edges);
    }
}

 void Refinement_scheme::collect_tet_edges(const Tetrahedron &tet, std::vector<Edge> &edges) {
    edges.resize(6);
    std::vector<int> nodes(tet.nodes.begin(), tet.nodes.end());
    edges[0].nodes = {nodes[0], nodes[1]};
    edges[1].nodes = {nodes[1], nodes[2]};
    edges[2].nodes = {nodes[2], nodes[0]};
    edges[3].nodes = {nodes[0], nodes[3]};
    edges[4].nodes = {nodes[1], nodes[3]};
    edges[5].nodes = {nodes[2], nodes[3]};
}

void Refinement_scheme::output_new_mesh(const std::set<Tetrahedron> &t,
                                        std::vector<std::vector<int>> &new_conn,
                                        std::vector<std::vector<double>> &new_coor,
                                        std::vector<int> &new_matid_arr,
                                        std::vector<int> &original) {
    // output new mesh
    new_coor.resize(m_coor_map.size());
    for (int node_id = 0; node_id < static_cast<int>(m_coor_map.size()); node_id++) {
        new_coor[node_id] = m_coor_map[node_id];
    }

    new_conn.resize(t.size());
    new_matid_arr.resize(t.size());
    original.resize(t.size());
    int ielem = 0;
    for (const auto &tet : t) {
        for (int node_id : tet.nodes) {
            new_conn[ielem].push_back(node_id);
        }
        std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
        for (int i = 0; i < 4; i++) {
            tetra[i] = new_coor[new_conn[ielem][i]];
        }
        double volume = findTetraVolume(tetra);
        while (volume < 1.e-10) {
            std::swap(new_conn[ielem][2], new_conn[ielem][3]);
            for (int i = 0; i < 4; i++) {
                tetra[i] = new_coor[new_conn[ielem][i]];
            }
            volume = findTetraVolume(tetra);
        }
        new_matid_arr[ielem] = tet.material;
        if (tet.original) {
            original[ielem] = 1;
        } else {
            original[ielem] = 0;
        }

        ielem += 1;
    }

    std::map<std::set<int>, std::vector<int>> face_to_elems;
    for (int elem_id = 0; elem_id < static_cast<int>(new_conn.size()); elem_id++) {
        std::vector<int> &elem = new_conn[elem_id];
        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(elem[(iface + inode + 1) % 4]);
            }
            face_to_elems[face].push_back(elem_id);
        }
    }
    for (const auto &p : face_to_elems) {
        const auto &face = p.first;
        if (p.second.size() > 2) {
            std::cerr << "ERROR: face_to_elems more than 2" << std::endl;
            std::cerr << "face: ";
            for (int inode : face) {
                std::cerr << inode << " ";
            }
            std::cerr << std::endl;
            std::cerr << "elems: ";
            for (int elem_id : p.second) {
                std::cerr << elem_id << " ";
            }
            std::cerr << std::endl;
            std::exit(1);
        }
        if (p.second.size() == 2) {
            continue;
        }
    }
}

void Refinement_scheme::build_neighbor_edges(const std::set<Tetrahedron> &t,
                                             std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges) {
    // neighbor_edges: neighbor_id -> vector of (edge, cut flag)
    std::set<Edge> edges_all;
    for (const auto &tet : t) {
        std::vector<Edge> edges;
        collect_tet_edges(tet, edges);
        for (const auto &edge : edges) {
            edges_all.insert(edge);
        }
    }

    std::map<int, std::set<int>> neighbor_nodes_set;
    for (const auto &p : m_neighbor_nodes) {
        int node_id = p.first;
        neighbor_nodes_set[node_id] = std::set<int>(p.second.begin(), p.second.end());
    }

    std::map<int, std::set<std::vector<int>>> neighbor_edges_tmp;
    // for (const auto &p : neighbor_nodes_set) {
    for (const auto &p : m_neighbor_nodes) {
        int neighbor_id = p.first;
        neighbor_edges_tmp[neighbor_id] = std::set<std::vector<int>>();
        const auto &nodes = p.second;
        for (const auto &edge : edges_all) {
            std::vector<int> com_ids;
            for (int node_id : edge.nodes) {
                auto it = std::find(nodes.begin(), nodes.end(), node_id);
                if (it != nodes.end()) {
                    com_ids.push_back(std::distance(nodes.begin(), it));
                }
            }
            if (com_ids.size() < 2) continue;
            std::sort(com_ids.begin(), com_ids.end());
            neighbor_edges_tmp.at(neighbor_id).insert(com_ids);
        }
    }

    // communication to find shared edges
    int nneighbor = m_neighbor_nodes.size();
    MPI_Request reqs[nneighbor][2];
    MPI_Status stats[nneighbor][2];
    std::map<int, int> send_counts, recv_counts;
    std::map<int, std::vector<int>> send_buf, recv_buf;
    int ineighbor = 0;

    for (const auto &p : neighbor_edges_tmp) {
        int neighbor_id = p.first;
        send_counts[neighbor_id] = p.second.size() * 2;
    }
    ineighbor = 0;
    for (const auto &p : neighbor_edges_tmp) {
        int rank = p.first;
        int tag = 0;
        MPI_Isend(&send_counts[rank], 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][0]);
        MPI_Irecv(&recv_counts[rank], 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);

    for (const auto &p : neighbor_edges_tmp) {
        int neighbor_id = p.first;
        send_buf[neighbor_id].resize(p.second.size() * 2);
        int cnt = 0;
        for (const auto &edge : p.second) {
            for (int node_id : edge) {
                send_buf[neighbor_id][cnt] = node_id;
                cnt++;
            }
        }
        recv_buf[neighbor_id].resize(recv_counts[neighbor_id]); 
    }

    ineighbor = 0;
    for (const auto &p : neighbor_edges_tmp) {
        int rank = p.first;
        int tag = 0;
        MPI_Isend(&send_buf[rank][0], send_buf[rank].size(), MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][0]);
        MPI_Irecv(&recv_buf[rank][0], recv_buf[rank].size(), MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);

    std::map<int, std::set<std::vector<int>>> neighbor_edges_other;
    for (const auto &p : neighbor_edges_tmp) {
        int neighbor_id = p.first;
        neighbor_edges_other[neighbor_id] = std::set<std::vector<int>>();
        for (int i = 0; i < static_cast<int>(recv_buf[neighbor_id].size()); i += 2) {
            std::vector<int> edge = {recv_buf[neighbor_id][i], recv_buf[neighbor_id][i + 1]};
            neighbor_edges_other[neighbor_id].insert(edge);
        }
    }

    for (const auto &p : neighbor_edges_tmp) {
        int neighbor_id = p.first;
        neighbor_edges[neighbor_id] = std::vector<std::tuple<std::vector<int>, int>>();
        const auto &edges = p.second;
        const auto &edges_other = neighbor_edges_other.at(neighbor_id);
        for (const auto &edge : edges) {
            if (!edges_other.count(edge)) continue;
            std::vector<int> edge_obj;
            for (int node_id : edge) {
                edge_obj.push_back(m_neighbor_nodes.at(neighbor_id).at(node_id));
            }
            neighbor_edges.at(neighbor_id).push_back(std::make_tuple(edge_obj, 0));
        }
    }

}

void Refinement_scheme::check_neighbor_edges(const std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges,
                                             const std::map<Edge, std::set<Tetrahedron>> &edge_to_tets) {
    // std::cout << "checking neighbor edges... " << std::endl;
    int nneighbor = neighbor_edges.size();
    MPI_Request reqs[nneighbor][2];
    MPI_Status stats[nneighbor][2];
    std::map<int, int> send_counts, recv_counts;
    std::map<int, std::vector<double>> send_buf, recv_buf;
    int ineighbor;

    for (const auto &p : neighbor_edges) {
        int neighbor_id = p.first;
        send_counts[neighbor_id] = p.second.size();
    }

    ineighbor = 0;
    for (const auto &p : neighbor_edges) {
        int rank = p.first;
        int tag = 0;
        MPI_Isend(&send_counts[rank], 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][0]);
        MPI_Irecv(&recv_counts[rank], 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);

    for (const auto &p : neighbor_edges) {
        int neighbor_id = p.first;
        if (recv_counts[neighbor_id] != send_counts[neighbor_id]) {
            std::cerr << "(626) ERROR: neighbor_edges size mismatch with neighbor " << neighbor_id << std::endl;
            std::exit(1);
        }
    }

    for (const auto &p : neighbor_edges) {
        int neighbor_id = p.first;
        send_buf[neighbor_id] = std::vector<double>(0);
        for (const auto &edge_tuple : p.second) {
            const auto &edge = std::get<0>(edge_tuple);
            for (int node_id : edge) {
                for (int idim = 0; idim < 3; idim++) {
                    send_buf[neighbor_id].push_back(m_coor_map[node_id][idim]);
                }
            }
        }
        recv_buf[neighbor_id].resize(recv_counts[neighbor_id] * 6);
    }

    ineighbor = 0;
    for (const auto &p : neighbor_edges) {
        int rank = p.first;
        int tag = 0;
        MPI_Isend(&send_buf[rank][0], send_buf[rank].size(), MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][0]);
        MPI_Irecv(&recv_buf[rank][0], recv_buf[rank].size(), MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);

    for (const auto &p : neighbor_edges) {
        int neighbor_id = p.first;
        const auto &edges = p.second;
        for (int i = 0; i < static_cast<int>(recv_buf[neighbor_id].size()); i++) {
            double v1 = recv_buf[neighbor_id][i];
            double v2 = send_buf[neighbor_id][i];
            if (std::abs(v1 - v2) > 100.0) {
                std::cerr << "ERROR: neighbor_edges mismatch with neighbor " << neighbor_id << std::endl;
                std::cerr << "recv_buf: " << v1 << ", send_buf: " << v2 << std::endl;
                std::exit(1);
            }
        }
    }

//     // edge_to_tets check
//     for (const auto &p : neighbor_edges) {
//         for (const auto &edge_tuple : p.second) {
//             const auto &edge_nodes = std::get<0>(edge_tuple);
//             Edge edge_obj;
//             edge_obj.nodes.insert(edge_nodes[0]);
//             edge_obj.nodes.insert(edge_nodes[1]);
//             if (!edge_to_tets.count(edge_obj)) {
//                 std::cerr << "ERROR: edge_to_tets does not contain edge" << std::endl;
//                 std::cerr << "edge nodes: ";
//                 for (int node_id : edge_obj.nodes) {
//                     std::cerr << node_id << " ";
//                 }
//                 std::cerr << std::endl;
//                 std::exit(1);
//             }
//         }
//     }
    // std::cout << "neighbor edges check passed" << std::endl;
}

void Refinement_scheme::sync_edges_cut(std::map<Edge, int> &edges_cut,
                                       std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges) {
    // std::cout << "synchronizing cut edges... " << std::endl;
    int nneighbor = neighbor_edges.size();
    MPI_Request reqs[nneighbor][2];
    MPI_Status stats[nneighbor][2];
    std::map<int, int> send_counts, recv_counts;
    std::map<int, std::vector<int>> send_buf, recv_buf;
    int ineighbor;

    for (const auto &p : edges_cut) {
        const Edge &edge = p.first;
        for (auto &q : neighbor_edges) {
            int neighbor_id = q.first;
            auto &edges = q.second;
            for (auto &edge_tuple : edges) {
                if (edge.nodes.count(std::get<0>(edge_tuple).at(0)) && edge.nodes.count(std::get<0>(edge_tuple).at(1))) {
                    std::get<1>(edge_tuple) = 1; // mark cut edge
                }
            }
        }
    }

    // communication for edges_cut
    for (const auto &p : neighbor_edges) {
        int neighbor_id = p.first;
        send_counts[neighbor_id] = p.second.size();
    }
    ineighbor = 0;
    for (const auto &p : neighbor_edges) {
        int rank = p.first;
        int tag = 0;
        MPI_Isend(&send_counts[rank], 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][0]);
        MPI_Irecv(&recv_counts[rank], 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);
    
    for (const auto &p : neighbor_edges) {
        int neighbor_id = p.first;
        if (recv_counts[neighbor_id] != send_counts[neighbor_id]) {
            std::cerr << "(731) ERROR: neighbor_edges size mismatch with process " << neighbor_id << " and " << m_myid << std::endl;
            std::cout << "recv_counts: " << recv_counts[neighbor_id] << ", send_counts: " << send_counts[neighbor_id] << std::endl;
            std::exit(1);
        }
        send_buf[neighbor_id].resize(send_counts[neighbor_id]);
        for (int i = 0; i < static_cast<int>(send_counts[neighbor_id]); i++) {
            const auto &edge_tuple = p.second.at(i);
            int cut_flag = std::get<1>(edge_tuple);
            send_buf[neighbor_id][i] = cut_flag;
        }
        recv_buf[neighbor_id].resize(recv_counts[neighbor_id]);
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);
    // MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Finalize();
    // std::exit(0);

    ineighbor = 0;
    for (const auto &p : neighbor_edges) {
        int rank = p.first;
        int tag = 0;
        MPI_Isend(&send_buf[rank][0], send_buf[rank].size(), MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][0]);
        MPI_Irecv(&recv_buf[rank][0], recv_buf[rank].size(), MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);

    for (auto &p : neighbor_edges) {
        int neighbor_id = p.first;
        auto &edges = p.second;
        // std::cout << "neighbor_edge(" << m_myid << ", " << neighbor_id << ") orgsize: " << edges.size() << std::endl;
        for (int i = 0; i < static_cast<int>(recv_buf[neighbor_id].size()); i++) {
            // std::cout << "send_buf[" << neighbor_id << "][" << i << "]: " << send_buf[neighbor_id][i] << ", recv_buf[" << neighbor_id << "][" << i << "]: " << recv_buf[neighbor_id][i]
            //           << " (edge: " << std::get<0>(edges.at(i)).at(0) << ", " << std::get<0>(edges.at(i)).at(1) << ")" << std::endl;

            if (send_buf[neighbor_id][i] == 0 && recv_buf[neighbor_id][i] == 1) {
                // new node
                const auto edge_to_cut = std::get<0>(edges.at(i));
                Edge edge_obj;
                edge_obj.nodes.insert(edge_to_cut[0]);
                edge_obj.nodes.insert(edge_to_cut[1]);
                if (!edges_cut.count(edge_obj)) {
                    std::vector<double> mid_point(3);
                    for (int node_id : edge_to_cut) {
                        for (int idim = 0; idim < 3; idim++) {
                            mid_point.at(idim) += m_coor_map.at(node_id).at(idim);
                        }
                    }
                    for (int idim = 0; idim < 3; idim++) {
                        mid_point.at(idim) /= 2.0;
                    }
                    int enode = m_coor_map.size();
                    m_coor_map[enode] = mid_point;
                    edges_cut[edge_obj] = enode;
                    // std::cout << "(iedge: " << i << ") new node " << enode << " created at edge (" << edge_to_cut[0] << ", " << edge_to_cut[1] << ") between process " << m_myid << " and " << neighbor_id << std::endl;
                }
            }
            if (send_buf[neighbor_id][i] == 1 || recv_buf[neighbor_id][i] == 1) {
                // update neighbor_edges
                const auto &edge = std::get<0>(edges.at(i));
                Edge edge_obj;
                edge_obj.nodes.insert(edge[0]);
                edge_obj.nodes.insert(edge[1]);

                int enode = edges_cut.at(edge_obj);
                std::vector<std::vector<int>> child_edges;
                child_edges.push_back({edge[0], enode});
                child_edges.push_back({edge[1], enode});
                for (const auto &child_edge : child_edges) {
                    bool found0 = (std::find(edges.begin(), edges.end(), std::make_tuple(child_edge, 0)) != edges.end());
                    bool found1 = (std::find(edges.begin(), edges.end(), std::make_tuple(child_edge, 1)) != edges.end());
                    if (!found0 && !found1) {
                        edges.push_back(std::make_tuple(child_edge, 0));
                        // std::cout << "(iedge: " << i << ") add edge (" << child_edge[0] << ", " << child_edge[1] << ") to neighbor_edges of process " << m_myid << " for neighbor " << neighbor_id << std::endl;
                    }
                }
            }
        }
        // std::cout << "neighbor_edge(" << m_myid << ", " << neighbor_id << ") update size: " << edges.size() << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // std::cout << "synchronization of cut edges done" << std::endl;
}

void Refinement_scheme::tmp_output(std::set<Tetrahedron> &s) {  
    std::string rank_id = std::to_string(m_myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "rdata/" + rank_id + ".tmp.txt";
    std::ofstream ofs(filename);
    for (const auto &tet : s) {
        std::vector<double> center(3, 0.);
        for (int node_id : tet.nodes) {
            for (int idim = 0; idim < 3; idim++) {
                center[idim] += m_coor_map[node_id][idim];
            }
        }
        for (int idim = 0; idim < 3; idim++) {
            center[idim] /= 4.0;
        }
        ofs << center[0] << " " << center[1] << " " << center[2] << std::endl;
    }
}

void Refinement_scheme::build_neighbor_nodes(const std::map<int, std::vector<std::tuple<std::vector<int>, int>>> &neighbor_edges,
                                             std::map<int, std::vector<int>> &new_neighbor_nodes) {
    for (const auto &p : neighbor_edges) {
        std::map<int, int> node_map;
        std::set<int> nodes_set;
        int neighbor_id = p.first;
        int cnt = 0;
        const auto &nodes = m_neighbor_nodes.at(neighbor_id);
        for (int node_id : nodes) {
            node_map[node_id] = cnt;
            nodes_set.insert(node_id);
            cnt++;
        }

        const auto &edges = p.second;
        for (const auto &edge_tuple : edges) {
            const auto &edge_nodes = std::get<0>(edge_tuple);
            for (int node_id : edge_nodes) {
                if (!node_map.count(node_id)) {
                    node_map[node_id] = cnt;
                    nodes_set.insert(node_id);
                    cnt++;
                }
            }
        }
        std::vector<int> nodes_vec(nodes_set.begin(), nodes_set.end());
        std::sort(nodes_vec.begin(), nodes_vec.end(), [&node_map](int a, int b) { return node_map[a] < node_map[b]; });
        new_neighbor_nodes[neighbor_id] = nodes_vec;
    }
}

void Refinement_scheme::check_neighbor_nodes(const std::map<int, std::vector<int>> &neighbor_nodes) {
    int nneighbor = neighbor_nodes.size();
    std::map<int, int> send_count, recv_count;
    std::map<int, std::vector<double>> send_buf, recv_buf;
    MPI_Request reqs[nneighbor][2];
    MPI_Status stats[nneighbor][2];
    int ineighbor;
    
    for (const auto &p : neighbor_nodes) {
        int neighbor_id = p.first;
        send_count[neighbor_id] = p.second.size();
    }

    ineighbor = 0;
    for (const auto &p : neighbor_nodes) {
        int rank = p.first;
        int tag = 0;
        MPI_Isend(&send_count[rank], 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][0]);
        MPI_Irecv(&recv_count[rank], 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);

    for (const auto &p : neighbor_nodes) {
        int neighbor_id = p.first;
        if (recv_count[neighbor_id] != send_count[neighbor_id]) {
            std::cerr << "(834) ERROR: neighbor_nodes size mismatch with neighbor " << neighbor_id << std::endl;
            std::exit(1);
        }
    }

    for (const auto &p : neighbor_nodes) {
        int neighbor_id = p.first;
        for (const auto &node_id : p.second) {
            for (int idim = 0; idim < 3; idim++) {
                send_buf[neighbor_id].push_back(m_coor_map[node_id][idim]);
            }
        }
        recv_buf[neighbor_id].resize(recv_count[neighbor_id] * 3);
    }

    ineighbor = 0;
    for (const auto &p : neighbor_nodes) {
        int rank = p.first;
        int tag = 0;
        MPI_Isend(&send_buf[rank][0], send_buf[rank].size(), MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][0]);
        MPI_Irecv(&recv_buf[rank][0], recv_buf[rank].size(), MPI_DOUBLE, rank, tag, MPI_COMM_WORLD, &reqs[ineighbor][1]);
        ineighbor++;
    }
    MPI_Waitall(2 * nneighbor, &reqs[0][0], &stats[0][0]);

    for (const auto &p : neighbor_nodes) {
        int neighbor_id = p.first;
        for (int i = 0; i < static_cast<int>(recv_buf[neighbor_id].size()); i++) {
            double v1 = recv_buf[neighbor_id][i];
            double v2 = send_buf[neighbor_id][i];
            if (std::abs(v1 - v2) > 100.0) {
                std::cerr << "ERROR: neighbor_nodes mismatch with neighbor " << neighbor_id << std::endl;
                std::cerr << "recv_buf: " << v1 << ", send_buf: " << v2 << std::endl;
                std::exit(1);
            }
        }
    }
}
