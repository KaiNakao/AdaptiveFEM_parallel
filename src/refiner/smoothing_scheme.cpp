#include "../inc/smoothing_scheme.hpp"

SmoothingScheme::SmoothingScheme(const std::set<int> &elem_smooth,
                                 std::vector<std::vector<int>> &connectivity,
                                 std::vector<std::vector<double>> &coordinates)
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
    double prev_aspect_ratio = 100;
    for (int iter=0; iter<max_iter; ++iter)
    {
        moveNodes();
        double max_aspect_ratio = findMaxAspectRatio();
        std::cout <<"max_aspect_ratio" << max_aspect_ratio << std::endl;
        if (prev_aspect_ratio<max_aspect_ratio)
        {
            break;
        }
        prev_aspect_ratio = max_aspect_ratio;
        //bool comp = compareMaxAspectRatio();
        //if (!comp)
        //{
        //    break;
        //}
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
        //std::cout << inode << std::endl;
        moving_direction = calculateMovingDirection(inode);
        //for (int idim=0; idim<3; ++idim)
        //{
            //std::cout << "mv " << moving_direction[idim] << std::endl;
            //new_coordinates[inode][idim] += lmd * moving_direction[idim];
        //}
        new_coordinates[inode][2] += lmd * moving_direction[2];
    }
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

// calculate max aspect_ratio
double SmoothingScheme::findMaxAspectRatio()
{
    double max_aspect_ratio = 0.0;
    double aspect_ratio;
    int num_elem_smooth = m_elem_smooth.size();
    for (int ielem=0; ielem<num_elem_smooth; ++ielem)
    {
        int elem_id = m_elem_smooth[ielem];
        std::vector<std::vector<double>> verts(4, std::vector<double>(3));
        for (int inode=0; inode<4; ++inode)
        {
            verts[inode] = m_coordinates[m_connectivity[elem_id][inode]];
        }

        aspect_ratio = findTetraAspectRatio(verts);
        
        if (aspect_ratio > max_aspect_ratio)
        {
            max_aspect_ratio = aspect_ratio;
        } 
    }
    return max_aspect_ratio;
}

bool SmoothingScheme::compareMaxAspectRatio()
{
    int num_elem = m_connectivity.size();
    double max_aspect_ratio = 0.0;
    double external_max_aspect_ratio = 0.0;
    double aspect_ratio;
    for (int ielem=0; ielem<num_elem; ++ielem)
    {
        int flag = 0;
        for (int j=0; j<m_elem_smooth.size(); ++j)
        {
            if (ielem == m_elem_smooth[j])
            {
                flag = 1;
                break;
            }
        }

        std::vector<std::vector<double>> verts(4, std::vector<double>(3));
        for (int inode=0; inode<4; ++inode)
        {
            verts[inode] = m_coordinates[m_connectivity[ielem][inode]];
        }
        aspect_ratio = findTetraAspectRatio(verts);
        if (flag == 0)
        {
            if (aspect_ratio > external_max_aspect_ratio)
            {
                external_max_aspect_ratio = aspect_ratio;
            } 
        }
        else
        {
            if (aspect_ratio > max_aspect_ratio)
            {
                max_aspect_ratio = aspect_ratio;
            } 
        }
    }

    std::cout << max_aspect_ratio << " " << external_max_aspect_ratio << std::endl;
    if (max_aspect_ratio > external_max_aspect_ratio)
    {
        return false;
    }
    return true;
}