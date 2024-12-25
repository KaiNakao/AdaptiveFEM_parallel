#include "../inc/smoothing_scheme.hpp"

SmoothingScheme::SmoothingScheme(const std::set<int> &elem_smooth,
                                 std::vector<std::vector<int>> &connectivity,
                                 std::vector<std::vector<double>> &coordinates,
                                 std::vector<int> &matid_arr,
                                 std::map<std::set<int>, std::vector<int>> &face_to_elems,
                                 std::map<int, std::vector<std::set<int>>> &node_to_faces)
{
    std::cout << "smoothing constructor" << std::endl;
    std::vector<std::vector<int>> connectivity_tmp(connectivity.size(), std::vector<int>(4));
    std::vector<std::vector<double>> coordinates_tmp(coordinates.size(), std::vector<double>(3));
    for (int ielem = 0; ielem < connectivity.size(); ielem++) 
    {
        connectivity_tmp[ielem] = connectivity[ielem];
    }
    for (int inode = 0; inode < coordinates.size(); inode++) 
    {
        coordinates_tmp[inode] = coordinates[inode];
    }
    m_connectivity = connectivity_tmp;
    m_coordinates = coordinates_tmp;
    for (int elem_id : elem_smooth) {
        m_elem_smooth.push_back(elem_id);
    }
    for (int ielem = 0; ielem < connectivity.size(); ielem++) {
        m_matid_map[ielem] = matid_arr[ielem];
    }
    m_node_to_faces = node_to_faces;
    m_face_to_elems = face_to_elems;
}

SmoothingScheme::~SmoothingScheme()
{

}

// MAIN ROUTINE for node smoothing
void SmoothingScheme::executeSmoothing(std::vector<std::vector<int>> &new_conn, 
                                       std::vector<std::vector<double>> &new_coor)
{
    std::cout << "executing smoothing" << std::endl;
    int max_iter = 200;
    search_adj_nodes();
    fetchNodes();
    createMovableTabel();
    double prev_aspect_ratio = 100;
    double prev_norm = 0;
    for (int iter=0; iter<max_iter; ++iter)
    {
        moveNodes();
        double max_aspect_ratio = findMaxAspectRatio();
        double euclidian_norm = findEuclidianNorm();
        if (iter%1 == 0)
        {
            std::cout << iter << " " <<"max_aspect_ratio: " << max_aspect_ratio;
            std::cout << " total distance: " << euclidian_norm << std::endl;
        }
        if (max_aspect_ratio > prev_aspect_ratio)
        {
            break;
        }
        prev_aspect_ratio = max_aspect_ratio;
        prev_norm = euclidian_norm;
    }

    // output result
    new_coor.resize(m_coordinates.size());
    for (int inode = 0; inode < m_coordinates.size(); inode++) 
    {
        new_coor[inode] = m_coordinates[inode];
    }
    new_conn.resize(m_connectivity.size());
    for (int ielem = 0; ielem < m_connectivity.size(); ielem++) {
        new_conn[ielem] = m_connectivity[ielem];
    }
}

// fetching nodes to be laplacian smoothed
void SmoothingScheme::fetchNodes()
{
    for (int ielem=0; ielem<m_elem_smooth.size(); ++ielem)
    {
        int elem = m_elem_smooth[ielem];
        for (int inode=0; inode<4; ++inode)
        {
            m_node_smooth.insert(m_connectivity[elem][inode]);
        }
    }
}

// calculating moving direction
std::vector<double> SmoothingScheme::calculateMovingDirection(int _node_id)
{
    std::set<int> adj_nodes = m_adj_nodes[_node_id];
    std::vector<double> moving_direction(3, sizeof(double));
    for (int idim=0; idim<3; ++idim)
    {
        moving_direction[idim] = 0.0;
    }
    int num_adj_nodes = 0;
    for (int inode : adj_nodes)
    {
        //std::cout << "node_id: " << inode << std::endl;
        for (int idim=0; idim<3; ++idim)
        {
            moving_direction[idim] += (m_coordinates[inode][idim] - m_coordinates[_node_id][idim]);
        }
        num_adj_nodes += 1;
    }
    //std::cout << "adj_nodes: " << num_adj_nodes << std::endl;
    for (int idim=0; idim<3; ++idim)
    {
        moving_direction[idim] = moving_direction[idim]/num_adj_nodes;
        //std::cout << moving_direction[idim] << num_adj_nodes << std::endl;
    }
    return moving_direction;
}

// move nodes
void SmoothingScheme::moveNodes()
{
    std::vector<std::vector<double>> new_coordinates(m_coordinates.size(), std::vector<double>(3));
    std::copy(m_coordinates.begin(), m_coordinates.end(), new_coordinates.begin());
    std::vector<double> moving_direction;
    double lmd = 1e-2;
    for (int inode : m_node_smooth)
    {
        if (m_node_movable[inode])
        {
            moving_direction = calculateMovingDirection(inode);
            for (int idim=0; idim<3; ++idim)
            {
              //std::cout << "mv " << moving_direction[idim] << std::endl;
              new_coordinates[inode][idim] += lmd * moving_direction[idim];
            }
            //new_coordinates[inode][2] += lmd * moving_direction[2];
        }
    }
    m_coordinates_prev = m_coordinates;
    m_coordinates = new_coordinates;
}

// search adjacent nodes
void SmoothingScheme::search_adj_nodes()
{
    int num_elem = m_connectivity.size();
    for (int ielem=0; ielem<num_elem; ++ielem)
    {
        const std::vector<int> elem = m_connectivity[ielem];
        for (int inode=0; inode<4; ++inode)
        {
            for (int jnode=0; jnode<3; ++jnode)
            {
                m_adj_nodes[elem[inode]].insert(elem[(inode+jnode+1)%4]);
            }
        }
    }
}

// judge if the node is movable
void SmoothingScheme::createMovableTabel()
{
    int count = 0;
    for (int node : m_node_smooth)
    {
        std::vector<std::set<int>> faces = m_node_to_faces[node];
        std::set<int> materials;
        int flag_boundary = 0;
        for (int iface=0; iface<faces.size(); ++iface)
        {
            std::set<int> face = faces[iface];
            std::vector<int> elems = m_face_to_elems[face];
            for (int ielem=0; ielem<elems.size(); ++ielem)
            {
                materials.insert(m_matid_map[elems[ielem]]);
            }
            if (elems.size() == 1)
            {
                flag_boundary = 1;
            }
        }
        
        if (materials.size()>1) // material boundary
        {
            m_node_movable[node] = false;
        }
        else if (flag_boundary == 1) // regional boundary
        {
            m_node_movable[node] = false;
        }
        else
        {
            m_node_movable[node] = true;
            count++;
        }
    }
    std::cout << "movable point count: " << count << " / " << m_node_smooth.size() <<std::endl;
}

// calculate max aspect_ratio
double SmoothingScheme::findMaxAspectRatio()
{
    double max_aspect_ratio = 0.0;
    double aspect_ratio, prev_aspect_ratio;
    int num_elem_smooth = m_elem_smooth.size();
    for (int ielem=0; ielem<num_elem_smooth; ++ielem)
    {
        int elem_id = m_elem_smooth[ielem];
        std::vector<std::vector<double>> verts(4, std::vector<double>(3));
        for (int inode=0; inode<4; ++inode)
        {
            verts[inode] = m_coordinates[m_connectivity[elem_id][inode]];
        }

        std::vector<std::vector<double>> prev_verts(4, std::vector<double>(3));
        for (int inode=0; inode<4; ++inode)
        {
            prev_verts[inode] = m_coordinates_prev[m_connectivity[elem_id][inode]];
        }

        aspect_ratio = findTetraAspectRatio(verts);
        prev_aspect_ratio = findTetraAspectRatio(prev_verts);
        
        if (aspect_ratio > max_aspect_ratio && prev_aspect_ratio != aspect_ratio)
        {
            max_aspect_ratio = aspect_ratio;
        } 
    }
    return max_aspect_ratio;
}

double SmoothingScheme::findEuclidianNorm()
{
    double euclidian_norm = 0.0;
    double dx, dy, dz;
    for (int inode=0; inode<m_coordinates.size(); ++inode)
    {
        dx = m_coordinates[inode][0] - m_coordinates_prev[inode][0];
        dy = m_coordinates[inode][1] - m_coordinates_prev[inode][1];
        dz = m_coordinates[inode][2] - m_coordinates_prev[inode][2];
        euclidian_norm += sqrt(dx*dx + dy*dy + dz*dz);
    }
    return euclidian_norm;
}