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
    std::cout << " Reading shape..." << std::endl;
    int num_elem, num_node_linear, num_node_quad, num_material, num_elem_marked;
    read_shape(data_dir, num_elem, num_node_linear, num_node_quad, num_material, num_elem_marked);

    std::cout << " Reading mesh..." << std::endl;
    std::vector<std::vector<int>> cny(num_elem, std::vector<int>(10));
    std::vector<std::vector<double>> coor(num_node_linear, std::vector<double>(3));
    std::vector<std::vector<int>> adj_elems(num_elem, std::vector<int>(4));
    std::vector<int> matid_arr(num_elem);
    std::map<std::set<int>, std::vector<int>> face_to_elems;
    std::map<int, std::vector<std::set<int>>> node_to_faces;
    read_mesh(data_dir, num_elem, num_node_linear, cny, coor, matid_arr);
    search_adj_element(cny, coor, adj_elems, face_to_elems, node_to_faces);
    m_face_to_elems = face_to_elems;
    m_node_to_faces = node_to_faces;
    m_coordinates = coor;
    m_adj_elements = adj_elems;
    m_matid_arr = matid_arr;
    std::vector<std::vector<int>> connectivity_tmp(num_elem, std::vector<int>(4));
    std::vector<int> elem_to_scheme(num_elem);
    for (int ielem=0; ielem<num_elem; ++ielem)
    {
        for (int inode=0; inode<4; ++inode)
        {
            connectivity_tmp[ielem][inode] = cny[ielem][inode];
        }
        elem_to_scheme[ielem] = 0;
    }
    m_connectivity = connectivity_tmp;
    m_elem_to_scheme = elem_to_scheme;

    std::cout << " Reading marked elements..." << std::endl;
    std::vector<int> marked_elems(num_elem_marked, sizeof(int));
    read_marked_elem(data_dir, num_elem_marked, marked_elems);
    m_marked_elems_id = marked_elems;

}

// MAIN ROUTINE for Refiner
void Refiner::executeRefinement()
{
    int num_points_added = m_points.size();
    for (int ipoint=0; ipoint<num_points_added; ++ipoint)
    {
        
    }
    std::set<int> elem_smooth_tmp;
    int num_elem_marked = m_marked_elems_id.size();

    for (int i=0; i<num_elem_marked; ++i)
    {
        int elem_id = m_marked_elems_id[i];
        int scheme = switchScheme(elem_id);
        if (scheme == 1) // refinement
        {
            m_elem_refine.insert(elem_id);
            elem_smooth_tmp.insert(elem_id);
            m_elem_to_scheme[elem_id] = 1;
        } 
        else if (scheme == 2) // smoothing
        {
            m_elem_refine.insert(elem_id);
            elem_smooth_tmp.insert(elem_id);
            m_elem_to_scheme[elem_id] = 2;
        }
    }

    for (int elem : elem_smooth_tmp)
    {
        m_elem_smooth.insert(elem);
        for (int i=0; i<4; ++i)
        {
            int adj_elem = m_adj_elements[elem][i];
            if (adj_elem != -1)
            {
                m_elem_smooth.insert(adj_elem);
            }
        }
    }

    // 1. Do smoothing
    std::vector<std::vector<int>> tmp_connectivity;
    std::vector<std::vector<double>> tmp_coordinates;
    SmoothingScheme smoothing(m_elem_smooth, m_connectivity, m_coordinates, m_matid_arr, m_face_to_elems, m_node_to_faces);
    smoothing.executeSmoothing(tmp_connectivity, tmp_coordinates);
    // 2. DO mesh refinement
    Refinement_scheme refinement_scheme = Refinement_scheme(m_elem_refine, tmp_connectivity, 
                                                            tmp_coordinates, m_matid_arr, 
                                                            m_face_to_elems);
    std::vector<int> new_matid_arr;
    std::vector<std::vector<int>> new_connectivity;
    std::vector<std::vector<double>> new_coordinates;
    std::vector<int> original;  // flag  0 を拾う
    refinement_scheme.executeRefinement_bisect(new_connectivity, new_coordinates, new_matid_arr, original);
    // 3. Do smoothing w.r.t. new elements
    std::set<int> elem_resmooth;
    for (int ielem=0; ielem<original.size(); ++ielem)
    {
        if (original[ielem] == 0)
        {
            elem_resmooth.insert(ielem);
        }
    }
    std::vector<std::vector<int>> new_adj_elements(new_connectivity.size(), std::vector<int>(4));
    std::map<std::set<int>, std::vector<int>> new_face_to_elems;
    std::map<int, std::vector<std::set<int>>> new_node_to_faces;
    search_adj_element( new_connectivity, new_coordinates, new_adj_elements, new_face_to_elems, new_node_to_faces);
    std::vector<std::vector<int>> renew_connectivity;
    std::vector<std::vector<double>> renew_coordinates;
    SmoothingScheme re_smoothing(elem_resmooth, new_connectivity, new_coordinates, new_matid_arr, new_face_to_elems, new_node_to_faces);
    re_smoothing.executeSmoothing(renew_connectivity, renew_coordinates);


    std::vector<double> volume_arr(renew_connectivity.size());
    std::vector<double> aspectratio_arr(renew_connectivity.size());
    for (int ielem = 0; ielem < renew_connectivity.size(); ielem++)
    {
        std::vector<std::vector<double>> tetra(4, std::vector<double>(3));
        for (int i=0; i<4; ++i) {
            tetra[i] = renew_coordinates[renew_connectivity[ielem][i]];
        }
        volume_arr[ielem] = findTetraVolume(tetra);
        aspectratio_arr[ielem] = findTetraAspectRatio(tetra);
    }
    std::cout << "minimum volume: " << *std::min_element(volume_arr.begin(), volume_arr.end()) << std::endl;
    std::cout << "maximum aspect ratio: " << *std::max_element(aspectratio_arr.begin(), aspectratio_arr.end()) << std::endl;

    // output new mesh
    std::ofstream ofs;
    ofs.open(m_data_dir + "new_coordinates.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_coordinates" << std::endl;
        return;
    }
    for (int inode = 0; inode < renew_coordinates.size(); inode++) {
        ofs.write(reinterpret_cast<const char*>(renew_coordinates[inode].data()), 3 * sizeof(double));
    }
    ofs.close();

    ofs.open(m_data_dir + "new_connectivity.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_connectivity" << std::endl;
        return;
    }
    for (int ielem = 0; ielem < renew_connectivity.size(); ielem++) {
        for (int inode = 0; inode < 4; inode++) {
            renew_connectivity[ielem][inode] += 1;
        }
        new_matid_arr[ielem] += 1;
        ofs.write(reinterpret_cast<const char*>(renew_connectivity[ielem].data()), 4 * sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&new_matid_arr[ielem]), sizeof(int));
    }
    ofs.close();

    ofs.open(m_data_dir + "new_volume.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_volume" << std::endl;
        return;
    }
    for (int ielem = 0; ielem < volume_arr.size(); ielem++) {
        ofs.write(reinterpret_cast<const char*>(&volume_arr[ielem]), sizeof(double));
    }
    ofs.close();

    ofs.open(m_data_dir + "new_aspectratio.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_aspectratio" << std::endl;
        return;
    }
    for (int ielem = 0; ielem < aspectratio_arr.size(); ielem++) {
        ofs.write(reinterpret_cast<const char*>(&aspectratio_arr[ielem]), sizeof(double));
    }
    ofs.close();

    ofs.open(m_data_dir + "new_original.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_original" << std::endl;
        return;
    }
    for (int ielem = 0; ielem < original.size(); ielem++) {
        ofs.write(reinterpret_cast<const char*>(&original[ielem]), sizeof(int));
    }
    ofs.close();

    //output aspect ratio
    ofs.open(m_data_dir + "aspect_ratio.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for aspect_ratio" << std::endl;
        return;
    }
    for (int ielem = 0; ielem < m_connectivity.size(); ielem++) {
        double aspect_ratio;
        std::vector<std::vector<double>> verts(4, std::vector<double>(3));
        for (int inode=0; inode<4; ++inode){
            verts[inode] = m_coordinates[m_connectivity[ielem][inode]];
        }
        aspect_ratio = findTetraAspectRatio(verts);
        ofs.write(reinterpret_cast<const char*>(&aspect_ratio), sizeof(double));
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
    //std::cout << "aspect_ratio: " << aspectratio << std::endl;
    if (aspectratio < 4.)
    {
        return 1;
    }
    return 2;
}

Refiner::~Refiner()
{
    
}