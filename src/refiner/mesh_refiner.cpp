#include "../inc/mesh_refiner.hpp"
#include "../inc/adj_elem.hpp"

Refiner::Refiner()
{
    m_connectivity.clear();
    m_coordinates.clear();
    m_error.clear();
}

Refiner::Refiner(const std::string &data_dir)
{
    std::cout << " Reading shape..." << std::endl;
    int num_elem, num_node_linear, num_node_quad, num_material, num_elem_marked;
    read_shape(data_dir, num_elem, num_node_linear, num_node_quad, num_material, num_elem_marked);

    std::cout << " Reading mesh..." << std::endl;
    std::vector<std::vector<int>> cny(num_elem, std::vector<int>(10));
    std::vector<std::vector<double>> coor(num_node_quad, std::vector<double>(3));
    std::vector<std::vector<int>> adj_elems(num_elem, std::vector<int>(4));
    std::vector<int> matid_arr(num_elem);
    std::map<std::set<int>, std::vector<int>> face_to_elems;
    read_mesh(data_dir, num_elem, num_node_quad, cny, coor, matid_arr);
    search_adj_element(cny, coor, adj_elems, face_to_elems);

    m_coordinates = coor;
    m_adj_elements = adj_elems;
    std::vector<std::vector<int>> connectivity_tmp(num_elem, std::vector<int>(4));
    for (int ielem=0; ielem<num_elem; ++ielem)
    {
        for (int inode=0; inode<4; ++inode)
        {
            connectivity_tmp[ielem][inode] = cny[ielem][inode];
        }
    }
    m_connectivity = connectivity_tmp;


    // std::cout << " Reading eta..." << std::endl;
    // std::vector<double> eta(num_elem, sizeof(double));
    // read_eta(data_dir, num_elem, eta);
    // m_error = eta;

    std::cout << " Reading marked elements..." << std::endl;
    std::vector<int> marked_elems(num_elem_marked, sizeof(int));
    read_marked_elem(data_dir, num_elem_marked, marked_elems);
    m_marked_elems_id = marked_elems;
}

// MAIN ROUTINE for Refiner
void Refiner::executeRefinement()
{
    // int num_points_added = m_points.size();
    // for (int ipoint=0; ipoint<num_points_added; ++ipoint)
    // {
        
    // }

    int num_elem_marked = m_marked_elems_id.size();
    for (int i=0; i<num_elem_marked; ++i)
    {
        int elem_id = m_marked_elems_id[i];
        std::cout << elem_id << " refinement is executed" << std::endl;
        int scheme = switchScheme(elem_id);
        std::cout << scheme << " is selected" << std::endl;
        if (scheme == 1)
        {
            m_elem_refine.insert(elem_id);
            addPoints(elem_id);
        } 
        else if (scheme == 2)
        {
            m_elem_smooth.insert(elem_id);
        }
    }

    // 1. DO node movement
    // 2. DO mesh refinement
}

int Refiner::switchScheme(int _elem_id)
{
    std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
    for (int i=0; i<4; ++i)
    {
        tetra[i] = m_coordinates[m_connectivity[_elem_id][i]];
    }
    double aspectratio = findTetraAspectRatio(tetra);
    std::cout << "aspect_ratio: " << aspectratio << std::endl;
    if (aspectratio < 7.)
    {
        return 1;
    }
    return 2;
}

void Refiner::addPoints(int _elem_id)
{
    std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
    for (int i=0; i<4; ++i)
    {
        tetra[i] = m_coordinates[m_connectivity[_elem_id][i]];
    }
    std::vector<double> point;
    point = findInCenter(tetra);
    m_points.push_back(point);
    // some another point adding scheme for interface
}

//--------------------------------------------//
//            REFINEMENT SCHEME               //
//--------------------------------------------//

// MAIN ROUTINE for refinement
void Refiner::elementRefine(int _elem_id)
{
    std::vector<std::vector<double>> verts(4, std::vector<double>(3));
    for (int i=0; i<4; ++i)
    {
        verts[i] = m_coordinates[m_connectivity[_elem_id][i]];
    }
    std::vector<double> in_center = findInCenter(verts);

    std::cout << "root element is splitted" << std::endl;
    split14(_elem_id, in_center);
    int new_id = m_connectivity.size();
    for (int i=0; i<4; ++i)
    {
        std::cout << "recursive flipping********************************** " << new_id-4+i << std::endl;
        performFlip(new_id-4+i);
        updateMeshData();
    }

    // deletion of root element
    std::cout << "deleting root element" << std::endl;
    for (int j=0; j<m_adj_elements.size(); ++j)
    {
        for (int k=0; k<4; ++k)
        {
            if (m_adj_elements[j][k]>_elem_id)
            {
                m_adj_elements[j][k] =  m_adj_elements[j][k]-1;
            }
        }
    }
    m_connectivity.erase(m_connectivity.begin()+_elem_id);
    m_adj_elements.erase(m_adj_elements.begin()+_elem_id);

    for (int j=0; j<m_marked_elems_id.size(); ++j)
    {
        if (m_marked_elems_id[j]>_elem_id)
        {
            m_marked_elems_id[j] = m_marked_elems_id[j]-1;
        }
    }
}

// check if flippable
bool Refiner::checkFlippable(int _o_tetra_id)
{
    if (_o_tetra_id == -1)
    {
        return false;
    }
    // Find paired tetra: newly added node is the last node in connectivity
    int p_tetra_id = m_adj_elements[_o_tetra_id][3];
    std::cout << "checking flippable " << _o_tetra_id << " " << p_tetra_id <<std::endl;
    if (p_tetra_id == -1)
    {
        return false;
    }
    std::cout << "connectivity size " << m_connectivity.size() << std::endl;
    // Fetch tetras to be flipped
    std::vector<int> simplex_o = m_connectivity[_o_tetra_id];
    std::vector<int> simplex_p = m_connectivity[p_tetra_id];
    std::cout << simplex_o[0] << " " << simplex_o[1] << " " << simplex_o[2] << " " << simplex_o[3] << std::endl;
    std::cout << simplex_p[0] << " " << simplex_p[1] << " " << simplex_p[2] << " " << simplex_p[3] << std::endl;
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
            o_tetra[i][j] = m_coordinates[simplex_o[i]][j];
        }
    }
    std::vector<std::vector<double>> p_tetra(4, std::vector<double>(3));
    for (int i=0; i<4; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            p_tetra[i][j] = m_coordinates[simplex_p[i]][j];
        }
    }

    // insphere check
    std::vector<double> d = p_tetra[did];
    std::vector<double> circum_center = findCircumCenter(o_tetra);
    double circum_radius = findCircumRadius(o_tetra);
    std::cout << "DISTANCES: " << circum_radius << " " << findLength(circum_center, d) << std::endl;
    if (findLength(circum_center, d) > circum_radius)  // temprally modified condition
    {
        return false;
    }
    return true;
}

// recursively perform flip23, if the pair is fluppable
void Refiner::performFlip(int _elem_id)
{
    if (checkFlippable(_elem_id) == false)
    {
        std::cout << "not flippable" <<std::endl;
        return;
    }
    std::cout << "judged flippable" <<std::endl;
    flip23(_elem_id);
    int new_id = m_connectivity.size()-1;
    performFlip(new_id-2);
    performFlip(new_id-1);
    performFlip(new_id);
}

// split root element into 4 tetras
void Refiner::split14(int _tetra_id, std::vector<double> &_pt)
{
    // Store root tetra
    std::vector<int> simplex = m_connectivity[_tetra_id];
    //m_deleted_elems.push_back(_tetra_id);

    // Add new node
    m_coordinates.push_back(_pt);
    int new_node_id = m_coordinates.size()-1;

    // create new tetra
    std::vector<int> t1, t2, t3, t4;
    t1 = {simplex[0], simplex[1], simplex[2], new_node_id};
    t2 = {simplex[0], simplex[2], simplex[3], new_node_id};
    t3 = {simplex[0], simplex[3], simplex[1], new_node_id};
    t4 = {simplex[3], simplex[2], simplex[1], new_node_id};
    int t1_id, t2_id, t3_id, t4_id;
    m_connectivity.push_back(t1);
    t1_id = m_connectivity.size()-1;
    m_connectivity.push_back(t2);
    t2_id = m_connectivity.size()-1;
    m_connectivity.push_back(t3);
    t3_id = m_connectivity.size()-1;
    m_connectivity.push_back(t4);
    t4_id = m_connectivity.size()-1;
    //std::cout << t1_id << " " << t1[0] << " " << t1[1] << " " << t1[2] << " " << t1[3] << std::endl;
    //std::cout << t2_id << " " << t2[0] << " " << t2[1] << " " << t2[2] << " " << t2[3] << std::endl;
    //std::cout << t3_id << " " << t3[0] << " " << t3[1] << " " << t3[2] << " " << t3[3] << std::endl;
    //std::cout << t4_id << " " << t4[0] << " " << t4[1] << " " << t4[2] << " " << t4[3] << std::endl;
    std::cout << t1_id << " " << t2_id << " " << t3_id << " " << t4_id << " are generated" << std::endl;

    // update adj_elements list
    int n1, n2, n3, n4;
    n1 = m_adj_elements[_tetra_id][0];
    n2 = m_adj_elements[_tetra_id][1];
    n3 = m_adj_elements[_tetra_id][2];
    n4 = m_adj_elements[_tetra_id][3];
    replaceAdjElem(n1, _tetra_id, t4_id);
    replaceAdjElem(n2, _tetra_id, t2_id);
    replaceAdjElem(n3, _tetra_id, t3_id);
    replaceAdjElem(n4, _tetra_id, t1_id);
    m_adj_elements.push_back({t4_id, t2_id, t3_id, n4});
    m_adj_elements.push_back({t4_id, t3_id, t1_id, n2});
    m_adj_elements.push_back({t4_id, t1_id, t2_id, n3});
    m_adj_elements.push_back({t1_id, t3_id, t2_id, n1});
    //std::cout << "flip14 neighbours" << std::endl;
    //std::cout << n1 << " " << m_connectivity[n1][0] << " " << m_connectivity[n1][1] << " " << m_connectivity[n1][2] << " " << m_connectivity[n1][3] << std::endl;
    //std::cout << n2 << " " << m_connectivity[n2][0] << " " << m_connectivity[n2][1] << " " << m_connectivity[n2][2] << " " << m_connectivity[n2][3] << std::endl;
    //std::cout << n3 << " " << m_connectivity[n3][0] << " " << m_connectivity[n3][1] << " " << m_connectivity[n3][2] << " " << m_connectivity[n3][3] << std::endl;
    //std::cout << n4 << " " << m_connectivity[n4][0] << " " << m_connectivity[n4][1] << " " << m_connectivity[n4][2] << " " << m_connectivity[n4][3] << std::endl;
}

// find paired tetra and flip
void Refiner::flip23(int _o_tetra_id)
{
    // Find paired tetra: newly added node is the last node in connectivity
    int p_tetra_id = m_adj_elements[_o_tetra_id][3];

    // Fetch tetras to be flipped
    std::vector<int> simplex_o = m_connectivity[_o_tetra_id];
    std::vector<int> simplex_p = m_connectivity[p_tetra_id];
    int aid = findVertexPos(p_tetra_id, _o_tetra_id, 0);
    int bid = findVertexPos(p_tetra_id, _o_tetra_id, 1);
    int cid = findVertexPos(p_tetra_id, _o_tetra_id, 2);
    int did = 6 - (aid+bid+cid);
    m_deleted_elems.push_back(_o_tetra_id);
    m_deleted_elems.push_back(p_tetra_id);

    // create new tetra
    std::vector<int> t1, t2, t3;
    t1 = {simplex_o[0], simplex_o[1], simplex_p[did], simplex_o[3]};
    t2 = {simplex_o[1], simplex_o[2], simplex_p[did], simplex_o[3]};
    t3 = {simplex_o[2], simplex_o[0], simplex_p[did], simplex_o[3]};
    int t1_id, t2_id, t3_id;
    m_connectivity.push_back(t1);
    t1_id = m_connectivity.size()-1;
    m_connectivity.push_back(t2);
    t2_id = m_connectivity.size()-1;
    m_connectivity.push_back(t3);
    t3_id = m_connectivity.size()-1;
    //std::cout << t1_id << " " << t1[0] << " " << t1[1] << " " << t1[2] << " " << t1[3] << std::endl;
    //std::cout << t2_id << " " << t2[0] << " " << t2[1] << " " << t2[2] << " " << t2[3] << std::endl;
    //std::cout << t3_id << " " << t3[0] << " " << t3[1] << " " << t3[2] << " " << t3[3] << std::endl;
    std::cout << t1_id << " " << t2_id << " " << t3_id << " are generated" << std::endl;

    //update_adj_elements_list
    int n1, n2, n3, n4, n5, n6;
    n1 = m_adj_elements[_o_tetra_id][2];
    n2 = m_adj_elements[_o_tetra_id][0];
    n3 = m_adj_elements[_o_tetra_id][1];
    n4 = m_adj_elements[p_tetra_id][cid];
    n5 = m_adj_elements[p_tetra_id][aid];    
    n6 = m_adj_elements[p_tetra_id][bid];
    replaceAdjElem(n1, _o_tetra_id, t1_id);
    replaceAdjElem(n2, _o_tetra_id, t2_id);
    replaceAdjElem(n3, _o_tetra_id, t3_id);
    replaceAdjElem(n4, _o_tetra_id, t1_id);
    replaceAdjElem(n5, _o_tetra_id, t2_id);
    replaceAdjElem(n6, _o_tetra_id, t3_id);
    m_adj_elements.push_back({t2_id, t3_id, n1, n4});
    m_adj_elements.push_back({t3_id, t1_id, n2, n5});
    m_adj_elements.push_back({t1_id, t2_id, n3, n6});
    // std::cout << "flip23 neighbours" << std::endl;
    //std::cout << n1 << " " << m_connectivity[n1][0] << " " << m_connectivity[n1][1] << " " << m_connectivity[n1][2] << " " << m_connectivity[n1][3] << std::endl;
    //std::cout << n2 << " " << m_connectivity[n2][0] << " " << m_connectivity[n2][1] << " " << m_connectivity[n2][2] << " " << m_connectivity[n2][3] << std::endl;
    //std::cout << n3 << " " << m_connectivity[n3][0] << " " << m_connectivity[n3][1] << " " << m_connectivity[n3][2] << " " << m_connectivity[n3][3] << std::endl;
    //std::cout << n4 << " " << m_connectivity[n4][0] << " " << m_connectivity[n4][1] << " " << m_connectivity[n4][2] << " " << m_connectivity[n4][3] << std::endl;
    //std::cout << n5 << " " << m_connectivity[n5][0] << " " << m_connectivity[n5][1] << " " << m_connectivity[n5][2] << " " << m_connectivity[n5][3] << std::endl;
    //std::cout << n6 << " " << m_connectivity[n6][0] << " " << m_connectivity[n6][1] << " " << m_connectivity[n6][2] << " " << m_connectivity[n6][3] << std::endl;

    std::cout << "flipped!" << std::endl;
}

// replace _old_id tetra in adj_elements[_elem_id] with _new_id tetra 
void Refiner::replaceAdjElem(int _elem_id, int _old_id, int _new_id)
{
    for (int i=0; i<4; ++i)
    {
        if (m_adj_elements[_elem_id][i] == _old_id)
        {
            m_adj_elements[_elem_id][i] = _new_id;
            break;
        }
    }
}

// find position of vertex in _vert_pos in _known_tetra w.r.t _target_tetra
int Refiner::findVertexPos(int _target_tetra, int _known_tetra, int _vert_pos)
{
    //std::cout << "searching for " << m_connectivity[_known_tetra][_vert_pos] << std::endl;
    for (int i=0; i<4; ++i)
    {
        if (m_connectivity[_target_tetra][i] == m_connectivity[_known_tetra][_vert_pos])
        {
            //std::cout << "position found " << i << std::endl;
            return i;
        }
    }
    std::cout << "ERROR: Position not found" << std::endl;
}

// updating mesh
void Refiner::updateMeshData()
{
    std::cout << "updating mesh" << std::endl;
    int num_deleted_elem = m_deleted_elems.size();
    for (int i=0; i<num_deleted_elem; ++i)
    {
        int deleted_id = m_deleted_elems[i];
        std::cout << deleted_id << std::endl;
        for (int j=0; j<m_adj_elements.size(); ++j)
        {
            for (int k=0; k<4; ++k)
            {
                if (m_adj_elements[j][k]>deleted_id)
                {
                    m_adj_elements[j][k] =  m_adj_elements[j][k]-1;
                }
            }
        }
        m_connectivity.erase(m_connectivity.begin()+deleted_id);
        m_adj_elements.erase(m_adj_elements.begin()+deleted_id);

        for (int j=0; j<m_marked_elems_id.size(); ++j)
        {
            if (m_marked_elems_id[j]>deleted_id)
            {
                m_marked_elems_id[j] = m_marked_elems_id[j]-1;
            }
        }
    }
    m_deleted_elems.clear();
}


//--------------------------------------------//
//             SMOOTHING SCHEME               //
//--------------------------------------------//

// MAIN ROUTINE for node smoothing


// calculating moving direction



// Destructor
Refiner::~Refiner()
{
    ;
}