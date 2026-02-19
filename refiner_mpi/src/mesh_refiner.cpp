#include "mesh_refiner.hpp"
#include "adj_elem.hpp"

Refiner::Refiner() {
    m_connectivity.clear();
    m_coordinates.clear();
}

// Refiner::Refiner(const std::string &data_dir) {
//     m_data_dir = data_dir;
// }
Refiner::Refiner(const int myid, const int numprocs) {
    m_myid = myid;
    m_numprocs = numprocs;
}

void Refiner::readData() {
    int num_elem, num_node_linear, num_material, num_elem_marked;
    read_shape(m_myid, num_elem, num_node_linear, num_material, num_elem_marked);
    m_nmaterial = num_material;

    int num_elem_total, num_elem_marked_total;
    MPI_Allreduce(&num_elem, &num_elem_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&num_elem_marked, &num_elem_marked_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    read_mesh(m_myid, num_elem, num_node_linear, m_connectivity, m_coordinates, m_matid_arr);
    read_marked_elem(m_myid, num_elem_marked, num_elem, m_marked_elems_id);
    read_neighbor(m_myid, m_neighbor_nodes);
}

// MAIN ROUTINE for Refiner
void Refiner::executeRefinement() {
    // Use bisection refinement only (no smoothing)
    Refinement_scheme refinement_scheme = Refinement_scheme(
        m_myid, m_numprocs,
        m_marked_elems_id, m_connectivity,
        m_coordinates, m_matid_arr, m_neighbor_nodes);
    refinement_scheme.executeRefinement_bisect(m_new_connectivity, m_new_coordinates, m_new_matid_arr, m_original, m_new_neighbor_nodes);
}

void Refiner::writeData() {
    write_shape(m_myid, m_new_coordinates.size(), m_new_connectivity.size(), m_nmaterial);
    write_mesh(m_myid, m_new_connectivity, m_new_coordinates, m_new_matid_arr);
    write_neighbor(m_myid, m_new_neighbor_nodes);
    // // output new mesh
    // std::ofstream ofs(m_data_dir + "new_coordinates.bin", std::ios::binary);
    // if (!ofs) {
    //     std::cerr << "Error opening file for new_coordinates" << std::endl;
    //     return;
    // }
    // for (int inode = 0; inode < static_cast<int>(m_new_coordinates.size()); inode++) {
    //     ofs.write(reinterpret_cast<const char*>(m_new_coordinates[inode].data()), 3 * sizeof(double));
    // }
    // ofs.close();

    // ofs.open(m_data_dir + "new_connectivity.bin", std::ios::binary);
    // if (!ofs) {
    //     std::cerr << "Error opening file for new_connectivity" << std::endl;
    //     return;
    // }
    // for (int ielem = 0; ielem < static_cast<int>(m_new_connectivity.size()); ielem++) {
    //     int nodes[4] = {
    //         m_new_connectivity[ielem][0] + 1,
    //         m_new_connectivity[ielem][1] + 1,
    //         m_new_connectivity[ielem][2] + 1,
    //         m_new_connectivity[ielem][3] + 1
    //     };
    //     int matid = m_new_matid_arr[ielem] + 1;
    //     ofs.write(reinterpret_cast<const char*>(nodes), 4 * sizeof(int));
    //     ofs.write(reinterpret_cast<const char*>(&matid), sizeof(int));
    // }
    // ofs.close();

    // ofs.open(m_data_dir + "new_original.bin", std::ios::binary);
    // if (!ofs) {
    //     std::cerr << "Error opening file for new_original" << std::endl;
    //     return;
    // }
    // for (int ielem = 0; ielem < static_cast<int>(m_original.size()); ielem++) {
    //     ofs.write(reinterpret_cast<const char*>(&m_original[ielem]), sizeof(int));
    // }
    // ofs.close();
}

Refiner::~Refiner() {}
