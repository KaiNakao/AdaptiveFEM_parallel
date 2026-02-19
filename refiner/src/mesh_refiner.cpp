#include "mesh_refiner.hpp"
#include "adj_elem.hpp"
#include <fstream>
#include <iostream>

Refiner::Refiner() {
    m_connectivity.clear();
    m_coordinates.clear();
}

Refiner::Refiner(const std::string &data_dir) {
    m_data_dir = data_dir;
}

void Refiner::readData() {
    std::cout << "Refiner::readData() start" << std::endl;
    int num_elem, num_node_linear, num_node_quad, num_material, num_elem_marked;
    read_shape(m_data_dir, num_elem, num_node_linear, num_node_quad, num_material, num_elem_marked);

    read_mesh(m_data_dir, num_elem, num_node_linear, m_connectivity, m_coordinates, m_matid_arr);
    read_marked_elem(m_data_dir, num_elem_marked, m_marked_elems_id);
    std::cout << "Refiner::readData() end" << std::endl;
}

// MAIN ROUTINE for Refiner
void Refiner::executeRefinement() {
    // Use bisection refinement only (no smoothing)
    Refinement_scheme refinement_scheme = Refinement_scheme(m_marked_elems_id, m_connectivity,
                                                            m_coordinates, m_matid_arr);
    refinement_scheme.executeRefinement_bisect(m_new_connectivity, m_new_coordinates, m_new_matid_arr, m_original);
}

void Refiner::writeNewMesh() {
    // output new mesh
    std::ofstream ofs(m_data_dir + "new_coordinates.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_coordinates" << std::endl;
        return;
    }
    for (int inode = 0; inode < static_cast<int>(m_new_coordinates.size()); inode++) {
        ofs.write(reinterpret_cast<const char*>(m_new_coordinates[inode].data()), 3 * sizeof(double));
    }
    ofs.close();

    ofs.open(m_data_dir + "new_connectivity.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_connectivity" << std::endl;
        return;
    }
    for (int ielem = 0; ielem < static_cast<int>(m_new_connectivity.size()); ielem++) {
        int nodes[4] = {
            m_new_connectivity[ielem][0] + 1,
            m_new_connectivity[ielem][1] + 1,
            m_new_connectivity[ielem][2] + 1,
            m_new_connectivity[ielem][3] + 1
        };
        int matid = m_new_matid_arr[ielem] + 1;
        ofs.write(reinterpret_cast<const char*>(nodes), 4 * sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&matid), sizeof(int));
    }
    ofs.close();

    ofs.open(m_data_dir + "new_original.bin", std::ios::binary);
    if (!ofs) {
        std::cerr << "Error opening file for new_original" << std::endl;
        return;
    }
    for (int ielem = 0; ielem < static_cast<int>(m_original.size()); ielem++) {
        ofs.write(reinterpret_cast<const char*>(&m_original[ielem]), sizeof(int));
    }
    ofs.close();
}

Refiner::~Refiner() {}
