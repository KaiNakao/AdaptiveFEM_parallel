#include "refinement_scheme.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <tuple>

#include "geo_util.hpp"


Refinement_scheme::Refinement_scheme(
    const std::set<int> &elem_refine,
    const std::vector<std::vector<int>> &connectivity,
    const std::vector<std::vector<double>> &coordinates,
    const std::vector<int> &matid_arr)
    : m_elem_refine(elem_refine), m_connectivity(connectivity) {
    for (int ielem = 0; ielem < static_cast<int>(connectivity.size()); ielem++) {
        m_matid_map[ielem] = matid_arr[ielem];
    }
    for (int inode = 0; inode < static_cast<int>(coordinates.size()); inode++) {
        m_coor_map[inode] = coordinates[inode];
    }
}

Refinement_scheme::~Refinement_scheme() {}

void Refinement_scheme::executeRefinement_bisect(
    std::vector<std::vector<int>> &new_conn,
    std::vector<std::vector<double>> &new_coor, std::vector<int> &new_matid_arr,
    std::vector<int> &original) {

    std::set<Tetrahedron> t, s;
    std::map<Edge, std::set<Tetrahedron>> edge_to_tets;
    std::cout << "initialization... " << std::endl;
    initialize_tets(t, s, edge_to_tets);
    std::cout << "local refinement... " << std::endl;
    local_refine(t, s, edge_to_tets);
    std::cout << "output new mesh... " << std::endl;
    output_new_mesh(t, new_conn, new_coor, new_matid_arr, original);
    std::cout << "done" << std::endl;
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
                                     std::map<Edge, std::set<Tetrahedron>> &edge_to_tets) {
    std::map<Edge, int> edges_cut;
    bisect_tets(t, s, edges_cut, edge_to_tets);
    refine_to_conformity(t, s, edges_cut, edge_to_tets);
    return;
}

void Refinement_scheme::bisect_tets(std::set<Tetrahedron> &t,
                                     std::set<Tetrahedron> &s,
                                     std::map<Edge, int> &edges_cut,
                                     std::map<Edge, std::set<Tetrahedron>> &edge_to_tets) {
    for (const Tetrahedron &tet : s) {
        bisect_tet(t, tet, edges_cut, edge_to_tets);
        // remove
        t.erase(tet);
    }
    s.clear();
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
        m_newnode_edge[enode] = tet.refinement_edge;
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
        // auto it = edge_to_tets.find(edge);
        // if (it == edge_to_tets.end()) {
        //     continue;
        // }
        // it->second.erase(tet.nodes);
        // if (it->second.empty()) {
        //     edge_to_tets.erase(it);
        // }
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

}

void Refinement_scheme::refine_to_conformity(std::set<Tetrahedron> &t,
                                            std::set<Tetrahedron> &s,
                                            std::map<Edge, int> &edges_cut,
                                            std::map<Edge, std::set<Tetrahedron>> &edge_to_tets) {
    // check hanging nodes (local)
    s.clear();
    for (const auto &p : edges_cut) {
        const Edge &edge = p.first;
        // /////////////
        // std::vector<double> edge_center(3, 0.);
        // for (int node_id : edge.nodes) {
        //     for (int idim = 0; idim < 3; idim++) {
        //         edge_center[idim] += m_coor_map[node_id][idim];
        //     }
        // }
        // for (int idim = 0; idim < 3; idim++) {
        //     edge_center[idim] /= 2.0;
        // }
        // std::cout << "edge center: " << edge_center[0] << " " << edge_center[1] << " " << edge_center[2] << std::endl;
        // /////////////
        for (const auto &tet : edge_to_tets[edge]) {
            s.insert(tet);
        }
    }
    std::cout << "number of elements with hanging nodes: " << s.size() << std::endl;
    // tmp_output(s);
    // std::exit(0);

    if (s.size() == 0) {
        return;
    } else {
        bisect_tets(t, s, edges_cut, edge_to_tets);
        refine_to_conformity(t, s, edges_cut, edge_to_tets);
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

void Refinement_scheme::tmp_output(const std::set<Tetrahedron> &s) {
    std::ofstream ofs("tmp.txt");
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
