#include "../inc/refinement_scheme.hpp"
#include "../inc/geo_util.hpp"
#include <iostream>

Refinement_scheme::Refinement_scheme(const std::set<int> &elem_refine,
                                     std::vector<std::vector<int>> &connectivity,
                                     std::vector<std::vector<double>> &coordinates,
                                     std::vector<int> &matid_arr,
                                     std::map<std::set<int>, std::vector<int>> &face_to_elems) {
    for (int ielem = 0; ielem < connectivity.size(); ielem++) {
        m_elems_in_mesh[ielem] = connectivity[ielem];
        m_matid_map[ielem] = matid_arr[ielem];
    }
    m_max_elem_id = m_elems_in_mesh.size() - 1;
    for (int inode = 0; inode < coordinates.size(); inode++) {
        m_coor_map[inode] = coordinates[inode];
    }

    std::map<std::set<int>, std::vector<int>> edge_to_elems;
    std::set<std::set<int>> edge_refine;
    for (const auto &elem : m_elems_in_mesh) {
        int elem_id = elem.first;
        std::vector<int> nodes = elem.second;
        std::vector<std::vector<int>> edge_id;
        edge_id = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
        for (int iedge = 0; iedge < 6; iedge++) {
            std::set<int> edge;
            for (int inode = 0; inode < 2; inode++) {
                edge.insert(nodes[edge_id[iedge][inode]]);
            }
            edge_to_elems[edge].push_back(elem_id);
            if (elem_refine.count(elem_id)) {
                edge_refine.insert(edge);
            }
        }
    }

    for (const auto &edge : edge_refine) {
        std::vector<double> mid_point(3);
        for (int inode : edge) {
            for (int idim = 0; idim < 3; idim++) {
                mid_point[idim] += m_coor_map[inode][idim];
            }
        }
        for (int idim = 0; idim < 3; idim++) {
            mid_point[idim] /= 2.0;
        }
        m_points.push_back(mid_point);
        for (int elem_id : edge_to_elems[edge]) {
            m_candidate_elems.insert(elem_id);
        }
    }
    m_face_to_elems = face_to_elems;
}

Refinement_scheme::~Refinement_scheme() {
}

void Refinement_scheme::executeRefinement(std::vector<std::vector<int>> &new_conn, 
                                          std::vector<std::vector<double>> &new_coor,
                                          std::vector<int> &new_matid_arr) {
    std::cout << "executeRefinement" << std::endl;
    std::cout << "m_points.size(): " << m_points.size() << std::endl;
    for (int ipoint = 0; ipoint < m_points.size(); ipoint++) {
        std::cout << std::endl;
        m_coor_map[m_coor_map.size()] = m_points[ipoint];
        std::cout << "add point: " << m_coor_map.size() - 1 << std::endl;
        // find elements that contains the point
        std::vector<int> elems_containing_point;
        for (int ielem : m_candidate_elems) {
            std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
            for (int inode = 0; inode < 4; inode++) {
                tetra[inode] = m_coor_map[m_elems_in_mesh[ielem][inode]];
            }
            std::vector<double> r_loc = normalizeLocTetra(tetra, m_points[ipoint]);

            double r1 = r_loc[0], r2 = r_loc[1], r3 = r_loc[2];
            if (r1 >= -1e-8 && r2 >= -1e-8 && r3 >= -1e-8 && r1 + r2 + r3 <= 1.0 + 1e-8) {
                elems_containing_point.push_back(ielem);
            }
        }

        if (elems_containing_point.size() == 0) {
            std::cout << "ERROR: No element contains the point" << std::endl;
            std::exit(1);
        }

        std::set<int> elems_for_flip;
        for (int ielem = 0; ielem < elems_containing_point.size(); ielem++) {
            int elem_id = elems_containing_point[ielem];
            split14(elems_containing_point[ielem], m_points[ipoint], elems_for_flip);
        }

        for (int elem_id : elems_for_flip) {
            performFlip(elem_id);
        }
        std::cout << "number of elements: " << m_elems_in_mesh.size() << std::endl;
    }

    // for (const auto &p : m_face_to_elems) {
    //     const auto &face = p.first;
    //     const auto &elems = p.second;
    //     if (elems.size() != 1) {
    //         continue;
    //     }
    //     std::vector<std::vector<double>> tri;
    //     for (int inode : face) {
    //         tri.push_back(m_coor_map[inode]);
    //     }
    //     std::vector<double> center(3);
    //     for (int idim = 0; idim < 3; idim++) {
    //         center[idim] = (tri[0][idim] + tri[1][idim] + tri[2][idim]) / 3.0;
    //     }
    //     if (fabs(center[0] < 1e-8) || fabs(center[0] - 340000) < 1e-8) {
    //         continue;
    //     }
    //     if (fabs(center[1] < 1e-8) || fabs(center[1] - 440000) < 1e-8) {
    //         continue;
    //     }
    //     if (fabs(center[2] < 1e-8)) {
    //         continue;
    //     }

    //     std::cout << "face: ";
    //     for (int inode : face) {
    //         std::cout << inode << " ";
    //     }
    //     std::cout << "center: ";
    //     for (int idim = 0; idim < 3; idim++) {
    //         std::cout << center[idim] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // // check
    // for (const auto &p : m_face_to_elems) {
    //     if (p.second.size() > 2) {
    //         const auto &face = p.first;
    //         std::cout << "ERROR: face_to_elems more than 2" << std::endl;
    //         std::cout << "face: ";
    //         for (int inode : face) {
    //             std::cout << inode << " ";
    //         }
    //         std::cout << std::endl; 
    //         std::cout << "elems: ";
    //         for (int elem_id : p.second) {
    //             std::cout << elem_id << " ";
    //         }
    //         std::cout << std::endl;
    //         std::exit(1);
    //     }
    // }

    // output result
    new_coor.resize(m_coor_map.size());
    for (int inode = 0; inode < m_coor_map.size(); inode++) {
        new_coor[inode] = m_coor_map[inode];
    }

    new_conn.resize(m_elems_in_mesh.size());
    new_matid_arr.resize(m_matid_map.size());
    std::vector<int> elem_tmp;
    for (const auto &p : m_elems_in_mesh) {
        elem_tmp.push_back(p.first);
    }
    std::sort(elem_tmp.begin(), elem_tmp.end());
    for (int ielem = 0; ielem < elem_tmp.size(); ielem++) {
        new_conn[ielem] = m_elems_in_mesh[elem_tmp[ielem]];
        new_matid_arr[ielem] = m_matid_map[elem_tmp[ielem]];
    }
}

// split root element into 4 tetras
void Refinement_scheme::split14(int _tetra_id, std::vector<double> &_pt, std::set<int> &elems_for_flip)
{
    // Store root tetra
    std::vector<int> simplex = m_elems_in_mesh[_tetra_id];

    // Add new node
    int new_node_id = m_coor_map.size() - 1;

    // create new tetra
    std::vector<int> t1, t2, t3, t4;
    t1 = {simplex[0], simplex[1], simplex[2], new_node_id};
    t2 = {simplex[0], simplex[2], simplex[3], new_node_id};
    t3 = {simplex[0], simplex[3], simplex[1], new_node_id};
    t4 = {simplex[3], simplex[2], simplex[1], new_node_id};

    std::vector<std::vector<int>> t_arr;
    std::vector<std::set<int>> faces_to_remove;
    std::vector<int> t_id_arr;
    std::vector<std::vector<int>> t_arr_all = {t1, t2, t3, t4};
    std::vector<double> vol_arr_all(4);
    for (int ielem = 0; ielem < 4; ielem++) {
        std::vector<int> &t = t_arr_all[ielem];
        std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
        for (int inode = 0; inode < 4; inode++) {
            tetra[inode] = m_coor_map[t[inode]];
        }
        vol_arr_all[ielem] = findTetraVolume(tetra);
    }
    double mean_vol = 0.0;
    for (double vol : vol_arr_all) {
        mean_vol += vol;
    }
    mean_vol /= 4.0;

    for (int ielem = 0; ielem < 4; ielem++) {
        std::vector<int> &t = t_arr_all[ielem];
        if (vol_arr_all[ielem]/mean_vol < 1e-8) {
            faces_to_remove.push_back({t[0], t[1], t[2]});
            continue;
        }

        int t_id = m_max_elem_id + 1;
        m_elems_in_mesh[t_id] = t;
        m_candidate_elems.insert(t_id);
        m_matid_map[t_id] = m_matid_map[_tetra_id];
        t_arr.push_back(t);
        t_id_arr.push_back(t_id);
        elems_for_flip.insert(t_id);

        m_max_elem_id++;
    }

    std::cout << "split: " << _tetra_id << " -> ";
    for (int i = 0; i < t_arr.size(); i++) {
        std::cout << t_id_arr[i] << " ";
    }
    std::cout << std::endl;

    if (t_arr.size() == 0) {
        std::cout << "ERROR: No tetra generated" << std::endl;
        std::cout << "original tetra: " << _tetra_id << std::endl;
        for (int inode = 0; inode < 4; inode++) {
            std::cout << simplex[inode] << " ";
            for (int idim = 0; idim < 3; idim++) {
                std::cout << m_coor_map[simplex[inode]][idim] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << "new point: ";
        for (int idim = 0; idim < 3; idim++) {
            std::cout << _pt[idim] << " ";
        }
        std::exit(1);
    }

    for (const auto &t : t_arr) {
        std::set<int> face = {t[0], t[1], t[2]};
        std::vector<int> &elems = m_face_to_elems[face];
        elems.erase(std::remove(elems.begin(), elems.end(), _tetra_id), elems.end());
    }
    for (const auto &face : faces_to_remove) {
        m_face_to_elems.erase(face);
    }

    for (int ielem = 0; ielem < t_id_arr.size(); ielem++) {
        int elem_id = t_id_arr[ielem];
        std::vector<int> &elem = m_elems_in_mesh[elem_id];
        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(elem[(iface + inode + 1) % 4]);
            }
            m_face_to_elems[face].push_back(elem_id);
        }
    }

    // delete root element
    m_elems_in_mesh.erase(_tetra_id);
    m_candidate_elems.erase(_tetra_id);
    m_matid_map.erase(_tetra_id);
}

void Refinement_scheme::performFlip(int _elem_id)
{
    bool flippable = checkFlippable(_elem_id);
    if (!flippable)
    {
        return;
    }
    std::set<int> elems_for_flip;
    flip23(_elem_id, elems_for_flip);
    for (int elem_id : elems_for_flip) {
        performFlip(elem_id);
    }
}

// check if flippable
bool Refinement_scheme::checkFlippable(int _o_tetra_id)
{
    if (_o_tetra_id == -1)
    {
        return false;
    }
    if (m_elems_in_mesh.find(_o_tetra_id) == m_elems_in_mesh.end())
    {
        return false;
    }
    // Find paired tetra: newly added node is the last node in connectivity
    std::vector<int> simplex_o = m_elems_in_mesh.at(_o_tetra_id);
    std::set<int> face = {simplex_o.at(0), simplex_o.at(1), simplex_o.at(2)};
    if (m_face_to_elems[face].size() != 2) {
        return false;
    }
    int p_tetra_id = m_face_to_elems[face][0] == _o_tetra_id ? m_face_to_elems[face][1] : m_face_to_elems[face][0];

    // material boundary 
    if (m_matid_map[_o_tetra_id] != m_matid_map[p_tetra_id]) {
        return false;
    }

    // Fetch tetras to be flipped
    std::vector<int> simplex_p = m_elems_in_mesh[p_tetra_id];
    int aid = findVertexPos(p_tetra_id, _o_tetra_id, 0);
    int bid = findVertexPos(p_tetra_id, _o_tetra_id, 1);
    int cid = findVertexPos(p_tetra_id, _o_tetra_id, 2);
    int did = 6 - (aid+bid+cid);
    std::vector<std::vector<double>> o_tetra(4, std::vector<double>(3));
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            o_tetra[i][j] = m_coor_map[simplex_o[i]][j];
        }
    }
    std::vector<std::vector<double>> p_tetra(4, std::vector<double>(3));
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            p_tetra[i][j] = m_coor_map[simplex_p[i]][j];
        }
    }

    std::vector<std::vector<double>> tmp_tri(3, std::vector<double>(3));
    std::vector<double> tmp_pt1(3), tmp_pt2(3);
    tmp_tri[0] = o_tetra[0];
    tmp_tri[1] = o_tetra[1];
    tmp_tri[2] = o_tetra[2];
    tmp_pt1 = o_tetra[3];
    tmp_pt2 = p_tetra[did];
    if (!checkConcave(tmp_tri, tmp_pt1, tmp_pt2))
    {
        return false;
    }

    // insphere check
    std::vector<double> d = p_tetra[did];
    std::vector<double> circum_center = findCircumCenter(o_tetra);
    double circum_radius = findCircumRadius(o_tetra);
    if (findLength(circum_center, d) > circum_radius)  // temprally modified condition
    {
        return false;
    }
    return true;
}

// find paired tetra and flip
void Refinement_scheme::flip23(int _o_tetra_id, std::set<int> &elems_for_flip)
{
    // Find paired tetra: newly added node is the last node in connectivity
    std::vector<int> simplex_o = m_elems_in_mesh[_o_tetra_id];
    std::set<int> face = {simplex_o[0], simplex_o[1], simplex_o[2]};
    int p_tetra_id = m_face_to_elems[face][0] == _o_tetra_id ? m_face_to_elems[face][1] : m_face_to_elems[face][0];

    // Fetch tetras to be flipped
    std::vector<int> simplex_p = m_elems_in_mesh[p_tetra_id];
    int aid = findVertexPos(p_tetra_id, _o_tetra_id, 0);
    int bid = findVertexPos(p_tetra_id, _o_tetra_id, 1);
    int cid = findVertexPos(p_tetra_id, _o_tetra_id, 2);
    int did = 6 - (aid+bid+cid);

    // create new tetra
    std::vector<int> t1, t2, t3;
    t1 = {simplex_o[0], simplex_o[1], simplex_p[did], simplex_o[3]};
    t2 = {simplex_o[1], simplex_o[2], simplex_p[did], simplex_o[3]};
    t3 = {simplex_o[2], simplex_o[0], simplex_p[did], simplex_o[3]};
    std::vector<std::vector<int>> t_arr_all = {t1, t2, t3};
    std::vector<double> vol_arr_all(3);
    for (int ielem = 0; ielem < 3; ielem++) {
        std::vector<int> &t = t_arr_all[ielem];
        std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
        for (int inode = 0; inode < 4; inode++) {
            tetra[inode] = m_coor_map[t[inode]];
        }
        vol_arr_all[ielem] = findTetraVolume(tetra);
    }
    double mean_vol = 0.0;
    for (double vol : vol_arr_all) {
        mean_vol += vol;
    }
    mean_vol /= 3.0;

    std::vector<std::vector<int>> t_arr;
    std::vector<std::set<int>> faces_to_remove;
    std::vector<int> t_id_arr;
    faces_to_remove.push_back({simplex_o[0], simplex_o[1], simplex_o[2]});

    // 4-4 flip
    for (int ielem = 0; ielem < 3; ielem++) {
        std::vector<int> &t = t_arr_all[ielem];
        if (vol_arr_all[ielem]/mean_vol > 1e-8) {
            continue;
        }

        std::set<int> p_face = {t[0], t[1], t[2]};
        std::set<int> o_face = {t[0], t[1], t[3]};
        if (m_face_to_elems[p_face].size() == 0) {
            std::cout << "ERROR: No element contains the p_face" << std::endl;
        }
        if (m_face_to_elems[o_face].size() == 0) {
            std::cout << "ERROR: No element contains the o_face" << std::endl;
        }
        const auto &p_elems = m_face_to_elems[p_face];
        const auto &o_elems = m_face_to_elems[o_face];

        if ((p_elems.size() == 1 && o_elems.size() == 2) ||
            (p_elems.size() == 2 && o_elems.size() == 1)) {
            // std::cout << "unflippable" << std::endl;
            return;
        }
        
        faces_to_remove.push_back(p_face);
        faces_to_remove.push_back(o_face);

        if (p_elems.size() == 1 && o_elems.size() == 1) {
            continue;
        } 

        int q_tetra_id = p_elems.at(0) == p_tetra_id ? p_elems.at(1) : p_elems.at(0);
        int r_tetra_id = o_elems.at(0) == _o_tetra_id ? o_elems.at(1) : o_elems.at(0);
        if (m_matid_map[q_tetra_id] != m_matid_map[r_tetra_id]) {
            // std::cout << "unflippable" << std::endl;
            return;
        }
        std::vector<int> simplex_q = m_elems_in_mesh[q_tetra_id];
        std::vector<int> simplex_r = m_elems_in_mesh[r_tetra_id];
        for (int inode = 0; inode < 3; inode++) {
            if (std::find(simplex_q.begin(), simplex_q.end(), simplex_r[inode]) == simplex_q.end()) {
                // std::cout << "unflippable" << std::endl;
                return;
            }
        }

        int aid = findVertexPos(q_tetra_id, r_tetra_id, 0);
        int bid = findVertexPos(q_tetra_id, r_tetra_id, 1);
        int cid = findVertexPos(q_tetra_id, r_tetra_id, 2);
        int did = 6 - (aid+bid+cid);

        int eid;
        for (int inode = 0; inode < 3; inode++) {
            if (std::find(simplex_o.begin(), simplex_o.end(), simplex_r[inode]) == simplex_o.end()) {
                eid = inode;
                break;
            }
        }

        std::vector<int> t4, t5;
        std::vector<int> t_id_arr_add;
        t4 = {t[0], simplex_r[eid], simplex_q[did], t[3]};
        t5 = {simplex_r[eid], t[1], simplex_q[did], t[3]};
        for (const auto &t : {t4, t5}) {
            int t_id = m_max_elem_id + 1;
            m_elems_in_mesh[t_id] = t;
            m_candidate_elems.insert(t_id);
            m_matid_map[t_id] = m_matid_map[r_tetra_id];
            elems_for_flip.insert(t_id);
            m_max_elem_id++;
            t_id_arr_add.push_back(t_id);

            // update face_to_elems
            std::set<int> face;
            face = {t[0], t[1], t[2]};
            std::vector<int> &elems_q = m_face_to_elems[face];
            elems_q.erase(std::remove(elems_q.begin(), elems_q.end(), q_tetra_id), elems_q.end());
            face = {t[0], t[1], t[3]};
            std::vector<int> &elems_r = m_face_to_elems[face];
            elems_r.erase(std::remove(elems_r.begin(), elems_r.end(), r_tetra_id), elems_r.end());
        }
        std::cout << "flip(4) : " << q_tetra_id << " " << r_tetra_id << " -> ";
        for (int i = 0; i < t_id_arr_add.size(); i++) {
            std::cout << t_id_arr_add[i] << " ";
        }
        std::cout << std::endl;

        for (int elem_id : t_id_arr_add) {
            std::vector<int> &elem = m_elems_in_mesh[elem_id];
            for (int iface = 0; iface < 4; iface++) {
                std::set<int> face;
                for (int inode = 0; inode < 3; inode++) {
                    face.insert(elem[(iface + inode + 1) % 4]);
                }
                m_face_to_elems[face].push_back(elem_id);
            }
        }

        // delete root element
        m_elems_in_mesh.erase(r_tetra_id);
        m_elems_in_mesh.erase(q_tetra_id);
        m_candidate_elems.erase(r_tetra_id);
        m_candidate_elems.erase(q_tetra_id);
        m_matid_map.erase(r_tetra_id);
        m_matid_map.erase(q_tetra_id);
    }

    for (int ielem = 0; ielem < 3; ielem++) {
        std::vector<int> &t = t_arr_all[ielem];
        if (vol_arr_all[ielem]/mean_vol < 1e-8) {
            continue;
        }
        int t_id = m_max_elem_id + 1;
        m_elems_in_mesh[t_id] = t;
        m_candidate_elems.insert(t_id);
        m_matid_map[t_id] = m_matid_map[_o_tetra_id];
        t_arr.push_back(t);
        t_id_arr.push_back(t_id);
        elems_for_flip.insert(t_id);
        m_max_elem_id++;
    }

    std::cout << "flip : " << _o_tetra_id << " " << p_tetra_id << " -> ";
    for (int i = 0; i < t_arr.size(); i++) {
        std::cout << t_id_arr[i] << " ";
    }
    std::cout << std::endl;

    // update face_to_elems
    for (const auto &t : t_arr) {
        std::set<int> face;
        face = {t[0], t[1], t[2]};
        std::vector<int> &elems_p = m_face_to_elems[face];
        elems_p.erase(std::remove(elems_p.begin(), elems_p.end(), p_tetra_id), elems_p.end());

        face = {t[0], t[1], t[3]};
        std::vector<int> &elems_o = m_face_to_elems[face];
        elems_o.erase(std::remove(elems_o.begin(), elems_o.end(), _o_tetra_id), elems_o.end());
    }
    for (const auto &face : faces_to_remove) {
        m_face_to_elems.erase(face);
    }

    for (int ielem = 0; ielem < t_id_arr.size(); ielem++) {
        int elem_id = t_id_arr[ielem];
        std::vector<int> &elem = m_elems_in_mesh[elem_id];
        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(elem[(iface + inode + 1) % 4]);
            }
            m_face_to_elems[face].push_back(elem_id);
        }
    }

    // delete root elements
    m_elems_in_mesh.erase(_o_tetra_id);
    m_elems_in_mesh.erase(p_tetra_id);
    m_candidate_elems.erase(_o_tetra_id);
    m_candidate_elems.erase(p_tetra_id);
    m_matid_map.erase(_o_tetra_id);
    m_matid_map.erase(p_tetra_id);
    
}

// find position of vertex in _vert_pos in _known_tetra w.r.t _target_tetra
int Refinement_scheme::findVertexPos(int _target_tetra, int _known_tetra, int _vert_pos)
{
    for (int i=0; i<4; ++i)
    {
        if (m_elems_in_mesh[_target_tetra][i] == m_elems_in_mesh[_known_tetra][_vert_pos])
        {
            return i;
        }
    }
    std::cout << "ERROR: Position not found" << std::endl;
    std::exit(1);
}

bool Refinement_scheme::checkConcave(const std::vector<std::vector<double>> &_tri,
                                     const std::vector<double> &_pt1,
                                     const std::vector<double> &_pt2) {
    // find projection of _pt1-_pt2 on _tri
    std::vector<double> v1 = {_tri[1][0] - _tri[0][0], _tri[1][1] - _tri[0][1], _tri[1][2] - _tri[0][2]};
    std::vector<double> v2 = {_tri[2][0] - _tri[0][0], _tri[2][1] - _tri[0][1], _tri[2][2] - _tri[0][2]};
    std::vector<double> nvec = {v1[1]*v2[2] - v1[2]*v2[1],
                                v1[2]*v2[0] - v1[0]*v2[2],
                                v1[0]*v2[1] - v1[1]*v2[0]};
    double nnorm = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);
    nvec[0] /= nnorm; nvec[1] /= nnorm; nvec[2] /= nnorm;

    double alpha = (nvec[0]*(_tri[0][0] - _pt1[0]) + nvec[1]*(_tri[0][1] - _pt1[1]) + nvec[2]*(_tri[0][2] - _pt1[2])) /
                   (nvec[0]*(_pt2[0] - _pt1[0]) + nvec[1]*(_pt2[1] - _pt1[1]) + nvec[2]*(_pt2[2] - _pt1[2]));
    std::vector<double> proj = {_pt1[0] + alpha*(_pt2[0] - _pt1[0]),
                                _pt1[1] + alpha*(_pt2[1] - _pt1[1]),
                                _pt1[2] + alpha*(_pt2[2] - _pt1[2])};

    // check if proj is inside _tri
    std::vector<std::vector<double>> dxdr(2, std::vector<double>(2));
    dxdr[0][0] = _tri[1][0] - _tri[0][0];
    dxdr[0][1] = _tri[2][0] - _tri[0][0];
    dxdr[1][0] = _tri[1][1] - _tri[0][1];
    dxdr[1][1] = _tri[2][1] - _tri[0][1];
    double detj = dxdr[0][0]*dxdr[1][1] - dxdr[0][1]*dxdr[1][0];
    std::vector<std::vector<double>> drdx(2, std::vector<double>(2));
    drdx[0][0] = dxdr[1][1]/detj;
    drdx[0][1] = -dxdr[0][1]/detj;
    drdx[1][0] = -dxdr[1][0]/detj;
    drdx[1][1] = dxdr[0][0]/detj;

    std::vector<double> dx(2);
    dx[0] = proj[0] - _tri[0][0];
    dx[1] = proj[1] - _tri[0][1];

    double r1 = drdx[0][0]*dx[0] + drdx[0][1]*dx[1];
    double r2 = drdx[1][0]*dx[0] + drdx[1][1]*dx[1];

    if (r1 >= -1e-8 && r2 >= -1e-8 && r1 + r2 <= 1.0 + 1e-8) {
        return true;
    } else {
        return false;
    }
}