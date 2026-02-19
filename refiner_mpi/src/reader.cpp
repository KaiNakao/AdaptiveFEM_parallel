#include "reader.hpp"


// void read_shape(const std::string &data_dir, int &nelem, int &nnode_linear,
//                 int &nnode_quad, int &nmaterial, int &nelem_marked) {
//     std::ifstream ifs(data_dir + "shape.dat");
//     std::cout << "path: " << data_dir + "shape.dat" << std::endl;
//     if (!ifs) {
//         std::cerr << "Error: cannot open file shape.dat" << std::endl;
//         return;
//     }

//     std::string buf;
//     ifs >> buf;  // read "nelem"
//     ifs >> nelem;
//     ifs >> buf;  // read "nnode_linear"
//     ifs >> nnode_linear;
//     ifs >> buf;  // read "nnode_quad"
//     ifs >> nnode_quad;
//     ifs >> buf;  // read "nmaterial"
//     ifs >> nmaterial;
//     ifs >> buf;  // read "nelem_marked"
//     ifs >> nelem_marked;

//     std::cout << "nelem: " << nelem << std::endl;
//     std::cout << "nnode_linear: " << nnode_linear << std::endl;
//     std::cout << "nnode_quad: " << nnode_quad << std::endl;
//     std::cout << "nmaterial: " << nmaterial << std::endl;
//     std::cout << "nelem_marked: " << nelem_marked << std::endl;
// }

void read_shape(const int &myid,
                int &nelem, int &nnode_linear,
                int &nmaterial, int &nelem_marked) {
    long fid = 0, count, size;
    // read cdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "cdata/" + rank_id + ".data.h5";

    hdf_open_file_(filename.c_str(), &fid);
    size = 3;
    int settings[size];
    hdf_read_int_array_(&fid, "/setting", settings, &size, &count);
    nnode_linear = settings[0];
    nelem = settings[1];
    hdf_close_file_(&fid);

    std::ifstream ifs;
    std::string buf;
    ifs.open("data/setting.dat");
    for (int iter = 0; iter < 5; iter++) {
        std::getline(ifs, buf);
    }
    std::getline(ifs, buf);
    nmaterial = std::stoi(buf);
    ifs.close();

    std::vector<int> marked_elems(nelem);
    ifs.open("mdata/" + rank_id + ".marked_elem.bin");
    if (!ifs) {
        std::cerr << "Error: cannot open file marked_elem.bin" << std::endl;
        return;
    }
    ifs.read(reinterpret_cast<char*>(marked_elems.data()), nelem * sizeof(int));
    ifs.close();
    nelem_marked = 0;
    for (int ielem = 0; ielem < nelem; ielem++) {
        if (marked_elems[ielem] == 1) {
            nelem_marked++;
        }
    }
}

// void read_mesh(const std::string &data_dir, const int &nelem, const int &nnode,
//                std::vector<std::vector<int>> &cny,
//                std::vector<std::vector<double>> &coor,
//                std::vector<int> &matid_arr) {

//     std::vector<int> buf_cny(5 * nelem);
//     std::vector<double> buf_coor(3 * nnode);
//     cny.resize(nelem, std::vector<int>(4));
//     coor.resize(nnode, std::vector<double>(3));
//     matid_arr.resize(nelem);

//     FILE *fp;
//     // Read cny
//     if ((fp = fopen((data_dir + "cny_linear.bin").c_str(), "r")) == NULL) {
//         std::cerr << "Error: cannot open file cny_linear.bin" << std::endl;
//         return;
//     }
//     fread(&buf_cny[0], sizeof(int), 5 * nelem, fp);
//     fclose(fp);
//     for (int ielem = 0; ielem < nelem; ielem++) {
//         for (int inode = 0; inode < 4; inode++) {
//             cny.at(ielem).at(inode) = buf_cny.at(ielem * 5 + inode) - 1;
//         }
//         matid_arr.at(ielem) = buf_cny.at(ielem * 5 + 4) - 1;
//     }

//     // Read coor
//     if ((fp = fopen((data_dir + "coor_linear.bin").c_str(), "r")) == NULL) {
//         std::cerr << "Error: cannot open file coor_linear.bin" << std::endl;
//         return;
//     }
//     fread(&buf_coor[0], sizeof(double), 3 * nnode, fp);
//     fclose(fp);
//     for (int inode = 0; inode < nnode; inode++) {
//         for (int idim = 0; idim < 3; idim++) {
//             coor.at(inode).at(idim) = buf_coor.at(inode * 3 + idim);
//         }
//     }
// }

void read_mesh(const int &myid, const int &nelem, const int &nnode,
               std::vector<std::vector<int>> &cny,
               std::vector<std::vector<double>> &coor,
               std::vector<int> &matid_arr) {
    long fid = 0, count, size;
    // read cdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "cdata/" + rank_id + ".data.h5";

    hdf_open_file_(filename.c_str(), &fid);
    cny.resize(nelem, std::vector<int>(4));
    matid_arr.resize(nelem);
    size = 5 * nelem;
    int conn_buf[size];
    hdf_read_int_array_(&fid, "/conn", conn_buf, &size, &count);
    for (int ielem = 0; ielem < nelem; ielem++) {
        for (int inode = 0; inode < 4; inode++) {
            cny.at(ielem).at(inode) = conn_buf[ielem * 5 + inode] - 1;
        }
        matid_arr.at(ielem) = conn_buf[ielem * 5 + 4] - 1;
    }

    size = 3 * nnode;
    coor.resize(nnode, std::vector<double>(3));
    double coor_buf[size];
    hdf_read_double_array_(&fid, "/coor", coor_buf, &size, &count);
    for (int inode = 0; inode < nnode; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            coor[inode][idim] = coor_buf[inode * 3 + idim];
        }
    }
    hdf_close_file_(&fid);
}

void read_displacement(const std::string &data_dir, const int &nnode,
                       std::vector<std::vector<double>> &displacement,
                       int iload) {
    std::vector<double> buf_displacement(3 * nnode);

    FILE *fp;
    std::string load_id = std::to_string(iload + 1);
    load_id.insert(load_id.begin(), 4 - load_id.size(), '0');
    // Read displacement
    // if ((fp = fopen((data_dir + "displacement_quad.bin").c_str(), "r")) ==
    if ((fp =
             fopen((data_dir + "displacement_quad_" + load_id + ".bin").c_str(),
                   "r")) == NULL) {
        std::cerr << "Error: cannot open file displacement.bin" << std::endl;
        return;
    }
    fread(&buf_displacement[0], sizeof(double), 3 * nnode, fp);
    fclose(fp);
    for (int inode = 0; inode < nnode; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            displacement[inode][idim] = buf_displacement[inode * 3 + idim];
        }
    }
}

void read_material(const std::string &data_dir,
                   std::vector<std::vector<double>> &material) {
    std::vector<double> buf_material(3 * material.size());
    FILE *fp;
    if ((fp = fopen((data_dir + "material.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file material.bin" << std::endl;
        return;
    }
    fread(&buf_material[0], sizeof(double), 3 * material.size(), fp);
    fclose(fp);

    for (int imaterial = 0; imaterial < static_cast<int>(material.size()); imaterial++) {
        double vp = buf_material[imaterial * 3];
        double vs = buf_material[imaterial * 3 + 1];
        double rho = buf_material[imaterial * 3 + 2];
        double lam = rho * (vp * vp - 2.0 * vs * vs);
        double mu = rho * vs * vs;
        material[imaterial][0] = lam;
        material[imaterial][1] = mu;
    }

    for (int imaterial = 0; imaterial < static_cast<int>(material.size()); imaterial++) {
        std::cout << "material: " << imaterial << " ";
        for (int idim = 0; idim < 2; idim++) {
            std::cout << material[imaterial][idim] << " ";
        }
        std::cout << std::endl;
    }
}

void read_load_elem(const std::string &data_dir, const int &nelem,
                    std::vector<int> &load_elem) {
    FILE *fp;
    if ((fp = fopen((data_dir + "load_elem.bin").c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file load_elem.bin" << std::endl;
        return;
    }
    fread(&load_elem[0], sizeof(int), nelem, fp);
    fclose(fp);
}

// void read_marked_elem(const std::string &data_dir, const int &nelem_marked,
//                       std::set<int> &marked_elem) {
//     FILE *fp;
//     if ((fp = fopen((data_dir + "marked_elem_quad.bin").c_str(), "r")) ==
//         NULL) {
//         std::cerr << "Error: cannot open file marked_elem_quad.bin"
//                   << std::endl;
//         return;
//     }
//     std::vector<int> buf_marked_elem(nelem_marked);
//     fread(&buf_marked_elem[0], sizeof(int), nelem_marked, fp);
//     fclose(fp);
//     for (int ielem = 0; ielem < nelem_marked; ielem++) {
//         marked_elem.insert(buf_marked_elem[ielem] - 1);
//     }
// }

void read_marked_elem(const int &myid, const int &nelem_marked,
                      const int &nelem, std::set<int> &marked_elem) {
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::ifstream ifs("mdata/" + rank_id + ".marked_elem.bin");
    if (!ifs) {
        std::cerr << "Error: cannot open file marked_elem.bin" << std::endl;
        return;
    }
    std::vector<int> buf_marked_elem(nelem);
    ifs.read(reinterpret_cast<char*>(buf_marked_elem.data()), nelem * sizeof(int));
    ifs.close();
    marked_elem.clear();
    for (int ielem = 0; ielem < nelem; ielem++) {
        if (buf_marked_elem.at(ielem) == 1) {
            marked_elem.insert(ielem);
        }
    }
    // verify
    if (marked_elem.size() != static_cast<size_t>(nelem_marked)) {
        std::cerr << "Error: number of marked elements does not match. "
                  << "Expected: " << nelem_marked << ", Actual: " << marked_elem.size() << std::endl;
    }
}

void read_refinement_edge(const std::string &data_dir,
                          const int &nelem,
                          std::vector<std::vector<int>> &refinement_edge) {
    std::vector<int> buf_refienment_edge(2 * nelem, -1);
    FILE *fp;
    if ((fp = fopen((data_dir + "refienment_edge.bin").c_str(), "r")) != NULL) {
        fread(&buf_refienment_edge[0], sizeof(int), 2 * nelem, fp);
        fclose(fp);
    } else {
        std::cout << "this is 1st iteration of refinement" << std::endl;
    }

    for (int ielem = 0; ielem < nelem; ielem++) {
        for (int inode = 0; inode < 2; inode++) {
            refinement_edge[ielem][inode] = buf_refienment_edge[ielem * 2 + inode];
        }
    }
}

void read_marked_edge(const std::string &data_dir,
                      const int &nelem,
                      std::vector<std::vector<std::vector<int>>> &marked_edge) {
    std::vector<int> buf_marked_edge(2 * 4 * nelem, -1);
    FILE *fp;
    if ((fp = fopen((data_dir + "marked_edge.bin").c_str(), "r")) != NULL) {
        fread(&buf_marked_edge[0], sizeof(int), 2 * 4 * nelem, fp);
        fclose(fp);
    } else {
        std::cout << "this is 1st iteration of refinement" << std::endl;
    }

    for (int ielem = 0; ielem < nelem; ielem++) {
        for (int iface = 0; iface < 4; iface++) {
            for (int inode = 0; inode < 2; inode++) {
                marked_edge[ielem][iface][inode] = buf_marked_edge[ielem * 8 + iface * 2 + inode];
            }
        }
    }
}

void read_tet_flag(const std::string &data_dir,
                   const int &nelem,
                   std::vector<bool> &tet_flag) {
    std::vector<bool> buf_tet_flag(nelem, false);
    FILE *fp;
    if ((fp = fopen((data_dir + "refienment_edge.bin").c_str(), "r")) != NULL) {
        std::vector<char> temp_buf(nelem, 0);
        size_t read_count = fread(temp_buf.data(), sizeof(char), nelem, fp);
        for (int i = 0; i < static_cast<int>(read_count); ++i) {
            buf_tet_flag[i] = static_cast<bool>(temp_buf[i]);
        }
        fclose(fp);
    } else {
        std::cout << "this is 1st iteration of refinement" << std::endl;
    }
    tet_flag = buf_tet_flag;
}

void read_nload(int &nload) {
    std::ifstream ifs("data/pointload.dat");
    if (!ifs) {
        std::cerr << "Error: cannot open file nload.dat" << std::endl;
        return;
    }

    std::string buf;
    getline(ifs, buf);  // read "nload"
    getline(ifs, buf);  // read nload
    nload = std::stoi(buf);
}

void read_neighbor(const int &myid,
                   std::map<int, std::vector<int>> &neighbor_nodes) {
    neighbor_nodes.clear();
    long fid = 0, count, size, start;
    // read cdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "cdata/" + rank_id + ".data.h5";

    hdf_open_file_(filename.c_str(), &fid);
    start = 0;
    size = 1;
    std::vector<int> mpinode_buf;
    mpinode_buf.resize(size);
    hdf_read_int_array_part_(&fid, "/MPInode", &start, &size,
                             &(mpinode_buf[0]));
    int nneighbor = mpinode_buf[0];

    std::vector<int> neighbor_size(nneighbor);
    std::vector<int> neighbor_rank(nneighbor);
    start += size;
    size = 2 * nneighbor;
    mpinode_buf.resize(size);
    hdf_read_int_array_part_(&fid, "/MPInode", &start, &size,
                             &(mpinode_buf[0]));
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        neighbor_rank.at(ineighbor) = mpinode_buf[2 * ineighbor];
        neighbor_size.at(ineighbor) = mpinode_buf[2 * ineighbor + 1];
    }

    start += size;
    size = std::accumulate(neighbor_size.begin(), neighbor_size.end(), 0);
    mpinode_buf.resize(size);
    hdf_read_int_array_part_(&fid, "/MPInode", &start, &size,
                             &(mpinode_buf[0]));
    int pt = 0;
    for (int ineighbor = 0; ineighbor < nneighbor; ineighbor++) {
        neighbor_nodes[neighbor_rank.at(ineighbor)].resize(
            neighbor_size.at(ineighbor));
        for (int i = 0; i < neighbor_size.at(ineighbor); i++) {
            neighbor_nodes.at(neighbor_rank.at(ineighbor)).at(i) =
                mpinode_buf[pt + i] - 1;
        }
        pt += neighbor_size.at(ineighbor);
    }
    hdf_close_file_(&fid);
}