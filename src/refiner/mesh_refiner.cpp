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
    m_data_dir = data_dir;
    std::cout << "data_dir: " << data_dir << std::endl;
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
    m_face_to_elems = face_to_elems;
    m_coordinates = coor;
    m_adj_elements = adj_elems;
    m_matid_arr = matid_arr;
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
        // std::cout << elem_id << " refinement is executed" << std::endl;
        int scheme = switchScheme(elem_id);
        // std::cout << scheme << " is selected" << std::endl;
        if (scheme == 1)
        {
            m_elem_refine.insert(elem_id);
            // addPoints(elem_id);
        } 
        else if (scheme == 2)
        {
            m_elem_smooth.insert(elem_id);
        }
    }

    Refinement_scheme refinement_scheme = Refinement_scheme(m_elem_refine, m_connectivity, 
                                                            m_coordinates, m_matid_arr, 
                                                            m_face_to_elems);
    std::vector<std::vector<int>> new_connectivity;
    std::vector<std::vector<double>> new_coordinates;
    std::vector<int> new_matid_arr;
    // 1. DO node movement
    // 2. DO mesh refinement
    refinement_scheme.executeRefinement(new_connectivity, new_coordinates, new_matid_arr);

    // output new mesh
    std::ofstream ofs;
    ofs.open(m_data_dir + "new_coordinates.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_coordinates" << std::endl;
        return;
    }
    for (int inode = 0; inode < new_coordinates.size(); inode++) {
        ofs.write(reinterpret_cast<const char*>(new_coordinates[inode].data()), 3 * sizeof(double));
    }
    ofs.close();

    ofs.open(m_data_dir + "new_connectivity.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_connectivity" << std::endl;
        return;
    }
    for (int ielem = 0; ielem < new_connectivity.size(); ielem++) {
        for (int inode = 0; inode < 4; inode++) {
            new_connectivity[ielem][inode] += 1;
        }
        new_matid_arr[ielem] += 1;
        ofs.write(reinterpret_cast<const char*>(new_connectivity[ielem].data()), 4 * sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&new_matid_arr[ielem]), sizeof(int));
    }
    ofs.close();
}

int Refiner::switchScheme(int _elem_id)
{
    std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
    for (int i=0; i<4; ++i)
    {
        tetra[i] = m_coordinates[m_connectivity[_elem_id][i]];
    }
    double aspectratio = findTetraAspectRatio(tetra);
    // std::cout << "aspect_ratio: " << aspectratio << std::endl;
    if (aspectratio < 7.)
    {
        return 1;
    }
    return 2;
}

// Destructor
Refiner::~Refiner()
{
    ;
}