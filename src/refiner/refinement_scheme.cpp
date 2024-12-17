#include "../inc/refinement_scheme.hpp"
#include "../inc/geo_util.hpp"
#include <iostream>

Refinement_scheme::Refinement_scheme(const std::set<int> &elem_refine,
                                     std::vector<std::vector<int>> &connectivity,
                                     std::vector<std::vector<double>> &coordinates,
                                     std::vector<int> &matid_arr,
                                     std::map<std::set<int>, std::vector<int>> &face_to_elems) {
    // std::cout << "constructor" << std::endl;
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

    // for (int node_id : m_elems_in_mesh[68498]) {
    //     std::cout << node_id << std::endl;
    // }
    // std::exit(0);

    // for (int elem_id : edge_to_elems[{13500, 13595}]) {
    //     std::cout << elem_id << std::endl;
    // }
    // std::exit(0);

    // for (int elem_id : elem_refine) {
    //     std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
    //     for (int inode = 0; inode < 4; inode++) {
    //         tetra[inode] = m_coor_map[m_elems_in_mesh[elem_id][inode]];
    //     }
    //     std::vector<double> point = findInCenter(tetra);
    //     m_points.push_back(point);
    //     m_candidate_elems.insert(elem_id);
    // }
    // std::cout << "constructor end" << std::endl;

    m_face_to_elems = face_to_elems;
}

Refinement_scheme::~Refinement_scheme() {
    // std::cout << "destructor" << std::endl;
}

void Refinement_scheme::executeRefinement(std::vector<std::vector<int>> &new_conn, 
                                          std::vector<std::vector<double>> &new_coor) {
    std::cout << "executeRefinement" << std::endl;
    std::cout << "m_points.size(): " << m_points.size() << std::endl;
    for (int ipoint = 0; ipoint < m_points.size(); ipoint++) {

        std::cout << "add point: ";
        for (int idim = 0; idim < 3; idim++) {
            std::cout << m_points[ipoint][idim] << " ";
        }
        std::cout << std::endl;

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
                std::cout << "(r1, r2, r3) : " << r1 << " " << r2 << " " << r3 << std::endl;
            }
            // 面上と辺上はあとで
        }

        if (elems_containing_point.size() == 0) {
            std::cout << "ERROR: No element contains the point" << std::endl;
            std::exit(1);
        }
        // if (elems_containing_point.size() > 1) {
        //     std::cout << "ERROR: Multiple elements contain the point" << std::endl;
        //     for (int ielem : elems_containing_point) {
        //         std::cout << ielem << std::endl;
        //         for (int inode = 0; inode < 4; inode++) {
        //             std::cout << m_elems_in_mesh[ielem][inode] << " ";
        //         }
        //         std::cout << std::endl;
        //     }
        //     std::exit(1);
        // }

        std::set<int> elems_for_flip;
        for (int ielem = 0; ielem < elems_containing_point.size(); ielem++) {
            int elem_id = elems_containing_point[ielem];
            split14(elems_containing_point[ielem], m_points[ipoint], elems_for_flip);
        }

        // for (int ielem = 0; ielem < 4; ielem++) {
        //     performFlip(m_max_elem_id - 3 + ielem);
        // }
        for (int elem_id : elems_for_flip) {
            performFlip(elem_id);
        }
        std::cout << "number of elements: " << m_elems_in_mesh.size() << std::endl;
    }

    std::cout << "elem 177084" << std::endl;
    for (int inode = 0; inode < 4; inode++) {
        int node_id = m_elems_in_mesh[177084][inode];
        for (int idim = 0; idim < 3; idim++) {
            std::cout << m_coor_map[node_id][idim] << " ";
        }
        std::cout << std::endl;
    }

    // output result
    new_coor.resize(m_coor_map.size());
    for (int inode = 0; inode < m_coor_map.size(); inode++) {
        new_coor[inode] = m_coor_map[inode];
    }

    new_conn.resize(m_elems_in_mesh.size());
    std::vector<int> elem_tmp;
    for (const auto &p : m_elems_in_mesh) {
        elem_tmp.push_back(p.first);
    }
    std::sort(elem_tmp.begin(), elem_tmp.end());
    for (int ielem = 0; ielem < elem_tmp.size(); ielem++) {
        new_conn[ielem] = m_elems_in_mesh[elem_tmp[ielem]];
        std::cout << "elem_id: " << elem_tmp[ielem] << " -> " << ielem << std::endl; 
    }
}

// split root element into 4 tetras
void Refinement_scheme::split14(int _tetra_id, std::vector<double> &_pt, std::set<int> &elems_for_flip)
{
    // Store root tetra
    std::vector<int> simplex = m_elems_in_mesh[_tetra_id];

    // Add new node
    int new_node_id = m_coor_map.size();
    m_coor_map[new_node_id] = _pt;

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

    // int t1_id, t2_id, t3_id, t4_id;
    // t1_id = m_max_elem_id + 1;
    // t2_id = m_max_elem_id + 2;
    // t3_id = m_max_elem_id + 3;
    // t4_id = m_max_elem_id + 4;
    // m_max_elem_id += 4;

    // m_elems_in_mesh[t1_id] = t1;
    // m_elems_in_mesh[t2_id] = t2;
    // m_elems_in_mesh[t3_id] = t3;
    // m_elems_in_mesh[t4_id] = t4;
    // m_candidate_elems.insert(t1_id);
    // m_candidate_elems.insert(t2_id);
    // m_candidate_elems.insert(t3_id);
    // m_candidate_elems.insert(t4_id);
    // m_matid_map[t1_id] = m_matid_map[_tetra_id];
    // m_matid_map[t2_id] = m_matid_map[_tetra_id];
    // m_matid_map[t3_id] = m_matid_map[_tetra_id];
    // m_matid_map[t4_id] = m_matid_map[_tetra_id];

    // std::cout << t1_id << " " << t2_id << " " << t3_id << " " << t4_id << " are generated" << std::endl;
    std::cout << "split: " << _tetra_id << " -> ";
    for (int i = 0; i < t_arr.size(); i++) {
        std::cout << t_id_arr[i] << " ";
    }
    std::cout << std::endl;
    //  << t1_id << " " << t2_id << " " << t3_id << " " << t4_id << std::endl;

    // update face_to_elems
    // std::vector<std::set<int>> face_arr;
    // face_arr.push_back({simplex[0], simplex[1], simplex[2]});
    // face_arr.push_back({simplex[0], simplex[2], simplex[3]});
    // face_arr.push_back({simplex[0], simplex[3], simplex[1]});
    // face_arr.push_back({simplex[1], simplex[2], simplex[3]});
    // for (int iface = 0; iface < 4; iface++) {
    //     std::set<int> &face = face_arr[iface];
    //     std::vector<int> &elems = m_face_to_elems[face];
    //     elems.erase(std::remove(elems.begin(), elems.end(), _tetra_id), elems.end());
    // }
    for (const auto &t : t_arr) {
        std::set<int> face = {t[0], t[1], t[2]};
        std::vector<int> &elems = m_face_to_elems[face];
        elems.erase(std::remove(elems.begin(), elems.end(), _tetra_id), elems.end());
    }
    for (const auto &face : faces_to_remove) {
        m_face_to_elems.erase(face);
    }

    // std::vector<int> id_arr = {t1_id, t2_id, t3_id, t4_id};
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
    if (checkFlippable(_elem_id) == false)
    {
        // std::cout << "not flippable" <<std::endl;
        return;
    }
    // std::cout << "judged flippable" <<std::endl;
    flip23(_elem_id);
    // int new_id = m_connectivity.size()-1;
    // performFlip(new_id-2);
    // performFlip(new_id-1);
    // performFlip(new_id);
    performFlip(m_max_elem_id - 2);
    performFlip(m_max_elem_id - 1);
    performFlip(m_max_elem_id);
}

// check if flippable
bool Refinement_scheme::checkFlippable(int _o_tetra_id)
{
    if (_o_tetra_id == -1)
    {
        return false;
    }
    // Find paired tetra: newly added node is the last node in connectivity
    std::vector<int> simplex_o = m_elems_in_mesh[_o_tetra_id];
    std::set<int> face = {simplex_o[0], simplex_o[1], simplex_o[2]};
    if (m_face_to_elems[face].size() != 2) {
        return false;
    }
    int p_tetra_id = m_face_to_elems[face][0] == _o_tetra_id ? m_face_to_elems[face][1] : m_face_to_elems[face][0];

    // material boundary 
    if (m_matid_map[_o_tetra_id] != m_matid_map[p_tetra_id]) {
        return false;
    }

    // std::cout << "checking flippable " << _o_tetra_id << " " << p_tetra_id <<std::endl;
    // if (p_tetra_id == -1)
    // {
    //     return false;
    // }
    // std::cout << "connectivity size " << m_connectivity.size() << std::endl;
    // Fetch tetras to be flipped
    std::vector<int> simplex_p = m_elems_in_mesh[p_tetra_id];

    // std::cout << simplex_o[0] << " " << simplex_o[1] << " " << simplex_o[2] << " " << simplex_o[3] << std::endl;
    // std::cout << simplex_p[0] << " " << simplex_p[1] << " " << simplex_p[2] << " " << simplex_p[3] << std::endl;

    int aid = findVertexPos(p_tetra_id, _o_tetra_id, 0);
    int bid = findVertexPos(p_tetra_id, _o_tetra_id, 1);
    int cid = findVertexPos(p_tetra_id, _o_tetra_id, 2);
    int did = 6 - (aid+bid+cid);
    //std::cout << aid << " " << bid << " " << cid << " " << did << std::endl;
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
    std::vector<double> tmp_pt(3);
    tmp_tri[0] = o_tetra[0];
    tmp_tri[1] = o_tetra[1];
    tmp_tri[2] = o_tetra[2];
    tmp_pt = o_tetra[3];
    if (checkConcave(tmp_tri, tmp_pt) == false) {
        return false;
    }
    tmp_tri[0] = o_tetra[0];
    tmp_tri[1] = o_tetra[2];
    tmp_tri[2] = o_tetra[1];
    tmp_pt = p_tetra[did];
    if (checkConcave(tmp_tri, tmp_pt) == false) {
        return false;
    }

    // insphere check
    std::vector<double> d = p_tetra[did];
    std::vector<double> circum_center = findCircumCenter(o_tetra);
    double circum_radius = findCircumRadius(o_tetra);
    // std::cout << "DISTANCES: " << circum_radius << " " << findLength(circum_center, d) << std::endl;
    if (findLength(circum_center, d) > circum_radius)  // temprally modified condition
    {
        return false;
    }
    return true;
}

// find paired tetra and flip
void Refinement_scheme::flip23(int _o_tetra_id)
{
    // Find paired tetra: newly added node is the last node in connectivity
    std::vector<int> simplex_o = m_elems_in_mesh[_o_tetra_id];
    std::set<int> face = {simplex_o[0], simplex_o[1], simplex_o[2]};
    int p_tetra_id = m_face_to_elems[face][0] == _o_tetra_id ? m_face_to_elems[face][1] : m_face_to_elems[face][0];
    // int p_tetra_id = m_adj_elements[_o_tetra_id][3];

    // Fetch tetras to be flipped
    // std::vector<int> simplex_o = m_connectivity[_o_tetra_id];
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
    int t1_id, t2_id, t3_id;
    t1_id = m_max_elem_id + 1;
    t2_id = m_max_elem_id + 2;
    t3_id = m_max_elem_id + 3;
    m_max_elem_id += 3;

    m_elems_in_mesh[t1_id] = t1;
    m_elems_in_mesh[t2_id] = t2;
    m_elems_in_mesh[t3_id] = t3;
    m_candidate_elems.insert(t1_id);
    m_candidate_elems.insert(t2_id);
    m_candidate_elems.insert(t3_id);
    m_matid_map[t1_id] = m_matid_map[_o_tetra_id];
    m_matid_map[t2_id] = m_matid_map[_o_tetra_id];
    m_matid_map[t3_id] = m_matid_map[_o_tetra_id];

    std::cout << "before flip23" << std::endl;
    std::cout << simplex_o[0] << " " << simplex_o[1] << " " << simplex_o[2] << " " << simplex_o[3] << std::endl;
    std::cout << simplex_p[0] << " " << simplex_p[1] << " " << simplex_p[2] << " " << simplex_p[3] << std::endl;
    std::cout << "after flip23" << std::endl;
    std::cout << t1_id << ":" << t1[0] << " " << t1[1] << " " << t1[2] << " " << t1[3] << std::endl;
    for (int i = 0; i < 4; i++) {
        std::cout << m_coor_map[t1[i]][0] << " " << m_coor_map[t1[i]][1] << " " << m_coor_map[t1[i]][2] << std::endl;
    }
    std::cout << t2_id << ":" << t2[0] << " " << t2[1] << " " << t2[2] << " " << t2[3] << std::endl;
    for (int i = 0; i < 4; i++) {
        std::cout << m_coor_map[t2[i]][0] << " " << m_coor_map[t2[i]][1] << " " << m_coor_map[t2[i]][2] << std::endl;
    }
    std::cout << t3_id << ":" << t3[0] << " " << t3[1] << " " << t3[2] << " " << t3[3] << std::endl;
    for (int i = 0; i < 4; i++) {
        std::cout << m_coor_map[t3[i]][0] << " " << m_coor_map[t3[i]][1] << " " << m_coor_map[t3[i]][2] << std::endl;
    }
    std::cout << t1_id << " " << t2_id << " " << t3_id << " are generated" << std::endl;
    std::cout << "flip : " << _o_tetra_id << " " << p_tetra_id << " -> " << t1_id << " " << t2_id << " " << t3_id << std::endl;

    // //update_adj_elements_list
    // int n1, n2, n3, n4, n5, n6;
    // n1 = m_adj_elements[_o_tetra_id][2];
    // n2 = m_adj_elements[_o_tetra_id][0];
    // n3 = m_adj_elements[_o_tetra_id][1];
    // n4 = m_adj_elements[p_tetra_id][cid];
    // n5 = m_adj_elements[p_tetra_id][aid];    
    // n6 = m_adj_elements[p_tetra_id][bid];
    // replaceAdjElem(n1, _o_tetra_id, t1_id);
    // replaceAdjElem(n2, _o_tetra_id, t2_id);
    // replaceAdjElem(n3, _o_tetra_id, t3_id);
    // replaceAdjElem(n4, _o_tetra_id, t1_id);
    // replaceAdjElem(n5, _o_tetra_id, t2_id);
    // replaceAdjElem(n6, _o_tetra_id, t3_id);
    // m_adj_elements.push_back({t2_id, t3_id, n1, n4});
    // m_adj_elements.push_back({t3_id, t1_id, n2, n5});
    // m_adj_elements.push_back({t1_id, t2_id, n3, n6});
    // std::cout << "flip23 neighbours" << std::endl;
    //std::cout << n1 << " " << m_connectivity[n1][0] << " " << m_connectivity[n1][1] << " " << m_connectivity[n1][2] << " " << m_connectivity[n1][3] << std::endl;
    //std::cout << n2 << " " << m_connectivity[n2][0] << " " << m_connectivity[n2][1] << " " << m_connectivity[n2][2] << " " << m_connectivity[n2][3] << std::endl;
    //std::cout << n3 << " " << m_connectivity[n3][0] << " " << m_connectivity[n3][1] << " " << m_connectivity[n3][2] << " " << m_connectivity[n3][3] << std::endl;
    //std::cout << n4 << " " << m_connectivity[n4][0] << " " << m_connectivity[n4][1] << " " << m_connectivity[n4][2] << " " << m_connectivity[n4][3] << std::endl;
    //std::cout << n5 << " " << m_connectivity[n5][0] << " " << m_connectivity[n5][1] << " " << m_connectivity[n5][2] << " " << m_connectivity[n5][3] << std::endl;
    //std::cout << n6 << " " << m_connectivity[n6][0] << " " << m_connectivity[n6][1] << " " << m_connectivity[n6][2] << " " << m_connectivity[n6][3] << std::endl;

    // update face_to_elems
    std::vector<std::set<int>> face_arr;
    face_arr.push_back({simplex_o[0], simplex_o[1], simplex_p[did]});
    face_arr.push_back({simplex_o[1], simplex_o[2], simplex_p[did]});
    face_arr.push_back({simplex_o[2], simplex_o[0], simplex_p[did]});
    for (int iface = 0; iface < 3; iface++) {
        std::set<int> &face = face_arr[iface];
        std::vector<int> &elems = m_face_to_elems[face];
        elems.erase(std::remove(elems.begin(), elems.end(), p_tetra_id), elems.end());
    }
    face_arr.clear();
    face_arr.push_back({simplex_o[0], simplex_o[1], simplex_o[2]});
    face_arr.push_back({simplex_o[0], simplex_o[2], simplex_o[3]});
    face_arr.push_back({simplex_o[0], simplex_o[3], simplex_o[1]});
    for (int iface = 0; iface < 3; iface++) {
        std::set<int> &face = face_arr[iface];
        std::vector<int> &elems = m_face_to_elems[face];
        elems.erase(std::remove(elems.begin(), elems.end(), _o_tetra_id), elems.end());
    }
    m_face_to_elems.erase({simplex_o[0], simplex_o[1], simplex_o[2]});

    std::vector<int> id_arr = {t1_id, t2_id, t3_id};
    for (int ielem = 0; ielem < id_arr.size(); ielem++) {
        int elem_id = id_arr[ielem];
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

    // std::cout << "flipped!" << std::endl;
}

// find position of vertex in _vert_pos in _known_tetra w.r.t _target_tetra
int Refinement_scheme::findVertexPos(int _target_tetra, int _known_tetra, int _vert_pos)
{
    //std::cout << "searching for " << m_connectivity[_known_tetra][_vert_pos] << std::endl;
    for (int i=0; i<4; ++i)
    {
        if (m_elems_in_mesh[_target_tetra][i] == m_elems_in_mesh[_known_tetra][_vert_pos])
        {
            //std::cout << "position found " << i << std::endl;
            return i;
        }
    }
    std::cout << "ERROR: Position not found" << std::endl;
}

bool Refinement_scheme::checkConcave(const std::vector<std::vector<double>> &_tri,
                                const std::vector<double> &_pt) {
    // find projection of _pt on _tri
    std::vector<double> v1 = {_tri[1][0] - _tri[0][0], _tri[1][1] - _tri[0][1], _tri[1][2] - _tri[0][2]};
    std::vector<double> v2 = {_tri[2][0] - _tri[0][0], _tri[2][1] - _tri[0][1], _tri[2][2] - _tri[0][2]};
    std::vector<double> wvec = {_pt[0] - _tri[0][0], _pt[1] - _tri[0][1], _pt[2] - _tri[0][2]};

    std::vector<double> nvec = {v1[1]*v2[2] - v1[2]*v2[1],
                                v1[2]*v2[0] - v1[0]*v2[2],
                                v1[0]*v2[1] - v1[1]*v2[0]};
    double nnorm = sqrt(nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2]);
    nvec[0] /= nnorm; nvec[1] /= nnorm; nvec[2] /= nnorm;

    double dist = wvec[0]*nvec[0] + wvec[1]*nvec[1] + wvec[2]*nvec[2];
    if (dist < 0) {
        std::cout << "ERROR: dist is negative" << std::endl;
    }

    std::vector<double> proj = {_pt[0] - dist*nvec[0], _pt[1] - dist*nvec[1], _pt[2] - dist*nvec[2]};

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