#include "error_estimator.hpp"

Error_Estimator::Error_Estimator(int myid, int numprocs) {
    // constructor
    this->myid = myid;
    this->numprocs = numprocs;
}

void Error_Estimator::prep() {
    nsp_error_estimator::read_shape(myid, nelem, nnode, nmaterial, nload);

    nsp_error_estimator::read_mesh(myid, nelem, nnode, cny, coor, matid_arr);

    // nsp_error_estimator::read_load_elem(myid, load_elem, nelem);
    // load_elem.resize(nelem);

    nsp_error_estimator::read_material(myid, nmaterial, material);

    std::map<int, std::vector<int>> neighbor_map;
    nsp_error_estimator::read_neighbor(myid, neighbor_map);

    nelem_org = cny.size();
    nnode_org = coor.size();
    nsp_error_estimator::extend_mesh(myid, neighbor_map, cny, coor, matid_arr, face_to_elems,
                import_index, export_index);

    for (const auto &p : face_to_elems) {
        face_arr.push_back(p.first);
    }

    eta_arr.resize(nelem_org);
}

void Error_Estimator::post() {
    std::vector<int> marked_elem;
    nsp_error_estimator::mark_elem_to_refine(eta_arr, load_elem, marked_elem, 0.01);

    nsp_error_estimator::write_eta(myid, eta_arr);
    nsp_error_estimator::write_marked_elem(myid, marked_elem);
}

void Error_Estimator::exec(int iload, const std::vector<double> &uvec) {
    std::vector<std::vector<double>> displacement(coor.size(),
                                                  std::vector<double>(3));
    for (int inode = 0; inode < nnode; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            displacement[inode][idim] = uvec[inode * 3 + idim];
        }
    }
    nsp_error_estimator::comm_displacement(myid, import_index, export_index, displacement);

    std::vector<double> rk_arr;
    std::map<std::set<int>, double> re_map;
    nsp_error_estimator::calc_rk(nelem_org, cny, coor, displacement, material, matid_arr,
            rk_arr);
    nsp_error_estimator::calc_re(cny, coor, displacement, material, matid_arr, face_to_elems,
            face_arr, re_map);
    for (int ielem = 0; ielem < nelem_org; ielem++) {
        double eta = 0.0;

        std::vector<int> elem = cny[ielem];
        std::vector<std::vector<double>> xnode(10, std::vector<double>(3));
        for (int inode = 0; inode < 10; inode++) {
            for (int idim = 0; idim < 3; idim++) {
                xnode[inode][idim] = coor[elem[inode]][idim];
            }
        }

        double ave_disp = 0.;
        for (int inode = 0; inode < 10; inode++) {
            double disp_norm = 0;
            for (int idim = 0; idim < 3; idim++) {
                disp_norm += displacement[elem[inode]][idim] * displacement[elem[inode]][idim];   ////節点変位を取ってくる
            }
            disp_norm = sqrt(disp_norm);
            ave_disp += disp_norm;
        }
        ave_disp = ave_disp / 10;


        double rk = rk_arr[ielem];
        double hk = nsp_error_estimator::calc_hk(xnode);

        eta += pow(hk, 2) * rk;

        for (int iface = 0; iface < 4; iface++) {
            std::set<int> face;
            for (int inode = 0; inode < 3; inode++) {
                face.insert(elem[(iface + inode) % 4]);
            }
            double beta = 1.0 / face_to_elems.at(face).size();
            double re = re_map.at(face);
            double he = nsp_error_estimator::calc_he(coor, face);

            eta += beta * he * re;
        }

        eta = sqrt(eta) / ave_disp;
        eta_arr[ielem] = std::max(eta_arr[ielem], eta); ////////////
    }
}

// int main(int argc, char **argv) {
//     int myid, numprocs;
//     MPI_Init(&argc, &argv);
//     MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//     MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

//     int nelem, nnode, nmaterial, nload;
//     read_shape(myid, nelem, nnode, nmaterial, nload);

//     std::vector<std::vector<int>> cny;
//     std::vector<int> matid_arr;
//     std::vector<std::vector<double>> coor;
//     read_mesh(myid, nelem, nnode, cny, coor, matid_arr);

//     std::vector<int> load_elem;
//     read_load_elem(myid, load_elem, nelem);

//     std::vector<std::vector<double>> material;
//     read_material(myid, nmaterial, material);

//     std::map<int, std::vector<int>> neighbor_map;
//     read_neighbor(myid, neighbor_map);

//     std::map<std::set<int>, std::set<int>> face_to_elems;
//     std::map<int, std::vector<int>> import_index, export_index;
//     int nelem_org = cny.size();
//     int nnode_org = coor.size();
//     extend_mesh(myid, neighbor_map, cny, coor, matid_arr, face_to_elems,
//                 import_index, export_index);

//     std::vector<double> eta_arr;
//     posterior_error_estimation(myid, nnode_org, nelem_org, import_index,
//                                export_index, cny, coor, material, matid_arr,
//                                face_to_elems, nload, eta_arr);

//     std::vector<int> marked_elem;
//     mark_elem_to_refine(eta_arr, load_elem, marked_elem, 0.01);

//     write_eta(myid, eta_arr);
//     write_marked_elem(myid, marked_elem);

//     MPI_Finalize();
// }
