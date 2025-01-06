#include <iostream>
#include <map>
#include <numeric>
#include <set>

#include "extend_mesh.hpp"
#include "mark_elem.hpp"
#include "mpi.h"
#include "posterior_error.hpp"
#include "reader.hpp"
#include "writer.hpp"

int main(int argc, char **argv) {
    int myid, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int nelem, nnode, nmaterial, nload;
    read_shape(myid, nelem, nnode, nmaterial, nload);

    std::vector<std::vector<int>> cny;
    std::vector<int> matid_arr;
    std::vector<std::vector<double>> coor;
    read_mesh(myid, nelem, nnode, cny, coor, matid_arr);

    std::vector<int> load_elem;
    read_load_elem(myid, load_elem, nelem);

    std::vector<std::vector<double>> material;
    read_material(myid, nmaterial, material);

    std::map<int, std::vector<int>> neighbor_map;
    read_neighbor(myid, neighbor_map);

    std::map<std::set<int>, std::set<int>> face_to_elems;
    std::map<int, std::vector<int>> import_index, export_index;
    int nelem_org = cny.size();
    int nnode_org = coor.size();
    extend_mesh(myid, neighbor_map, cny, coor, matid_arr, face_to_elems,
                import_index, export_index);

    std::vector<double> eta_arr;
    posterior_error_estimation(myid, nnode_org, nelem_org, import_index,
                               export_index, cny, coor, material, matid_arr,
                               face_to_elems, nload, eta_arr);

    std::vector<int> marked_elem;
    mark_elem_to_refine(eta_arr, load_elem, marked_elem, 0.01);

    write_eta(myid, eta_arr);
    write_marked_elem(myid, marked_elem);

    MPI_Finalize();
}
