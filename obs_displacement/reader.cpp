#include "reader.hpp"
extern "C" {
#include "hdf_lib.c"
}

void read_shape(const int &myid, int &nelem, int &nnode, int &nmaterial,
                int &nload) {
    long fid = 0, count, size;
    // read mdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "mdata/" + rank_id + ".data.h5";

    hdf_open_file_(filename.c_str(), &fid);
    size = 3;
    int settings[size];
    hdf_read_int_array_(&fid, "/setting", settings, &size, &count);
    nnode = settings[0];
    nelem = settings[1];

    std::ifstream ifs;
    std::string buf;
    ifs.open("data/setting.dat");
    for (int iter = 0; iter < 5; iter++) {
        std::getline(ifs, buf);
    }
    std::getline(ifs, buf);
    nmaterial = std::stoi(buf);
    ifs.close();

    ifs.open("data/pointload.dat");
    std::getline(ifs, buf);
    std::getline(ifs, buf);
    nload = std::stoi(buf);

    hdf_close_file_(&fid);
}

void read_mesh(const int &myid, const int &nelem, const int &nnode,
               std::vector<std::vector<int>> &cny,
               std::vector<std::vector<double>> &coor,
               std::vector<int> &matid_arr) {
    long fid = 0, count, size;
    // read mdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "mdata/" + rank_id + ".data.h5";

    hdf_open_file_(filename.c_str(), &fid);
    cny.resize(nelem, std::vector<int>(10));
    matid_arr.resize(nelem);
    size = 11 * nelem;
    int conn_buf[size];
    hdf_read_int_array_(&fid, "/conn", conn_buf, &size, &count);
    for (int ielem = 0; ielem < nelem; ielem++) {
        for (int inode = 0; inode < 10; inode++) {
            cny.at(ielem).at(inode) = conn_buf[ielem * 11 + inode] - 1;
        }
        matid_arr.at(ielem) = conn_buf[ielem * 11 + 10] - 1;
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

void read_load_elem(const int &myid, std::vector<int> &load_elem, int nelem) {
    FILE *fp;
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "mdata/" + rank_id + ".load_elem.bin";
    if ((fp = fopen(filename.c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file load_elem.bin" << std::endl;
        return;
    }
    load_elem.resize(nelem);
    fread(&load_elem[0], sizeof(int), nelem, fp);
    fclose(fp);
}

void read_material(const int &myid, const int &nmaterial,
                   std::vector<std::vector<double>> &material) {
    material.resize(nmaterial, std::vector<double>(2));
    double material_buf[2 * nmaterial];
    if (myid == 0) {
        std::string buf;
        std::ifstream ifs;
        ifs.open("data/material.dat");
        for (int imat = 0; imat < nmaterial; imat++) {
            std::getline(ifs, buf);
            std::getline(ifs, buf);
            double vp = std::stod(buf);
            std::getline(ifs, buf);
            double vs = std::stod(buf);
            std::getline(ifs, buf);
            double rho = std::stod(buf);
            for (int iter = 0; iter < 7; iter++) {
                std::getline(ifs, buf);
            }

            double lam = rho * (vp * vp - 2.0 * vs * vs);
            double mu = rho * vs * vs;
            material_buf[imat * 2] = lam;
            material_buf[imat * 2 + 1] = mu;
        }
    }
    MPI_Bcast(material_buf, 2 * nmaterial, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int imat = 0; imat < nmaterial; imat++) {
        material.at(imat).at(0) = material_buf[imat * 2];
        material.at(imat).at(1) = material_buf[imat * 2 + 1];
    }
}

void read_neighbor(const int &myid,
                   std::map<int, std::vector<int>> &neighbor_map) {
    long fid = 0, count, size, start;
    // read mdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "mdata/" + rank_id + ".data.h5";

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
        neighbor_map[neighbor_rank.at(ineighbor)].resize(
            neighbor_size.at(ineighbor));
        for (int i = 0; i < neighbor_size.at(ineighbor); i++) {
            neighbor_map.at(neighbor_rank.at(ineighbor)).at(i) =
                mpinode_buf[pt + i] - 1;
        }
        pt += neighbor_size.at(ineighbor);
    }
    hdf_close_file_(&fid);
}

void read_displacement(const int &myid, const int &iload, const int &nnode,
                       std::vector<std::vector<double>> &displacement) {
    std::string rank_id, load_id;
    rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    load_id = std::to_string(iload + 1);
    load_id.insert(load_id.begin(), 4 - load_id.size(), '0');
    std::string filename = "displacement/" + load_id + "_" + rank_id + ".bin";
    FILE *fp;
    double displacement_buf[3 * nnode];
    if ((fp = fopen(filename.c_str(), "r")) == NULL) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }
    fread(&displacement_buf[0], sizeof(double), 3 * nnode, fp);
    fclose(fp);
    for (int inode = 0; inode < nnode; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            displacement[inode][idim] = displacement_buf[inode * 3 + idim];
        }
    }
}

void read_obs_point(const int &myid,
                    std::vector<std::vector<double>> &obs_point) {
    int nobs;
    std::vector<double> obs_point_buf;
    if (myid == 0) {
        std::string filename, buf;
        filename = "data/obs_points.dat";
        std::ifstream ifs;
        ifs.open(filename);
        if (!ifs) {
            std::cerr << "Error: cannot open file " << filename << std::endl;
            return;
        }
        getline(ifs, buf);  // number of observation points
        getline(ifs, buf);
        nobs = std::stoi(buf);
        obs_point_buf.resize(3 * nobs);
        getline(ifs, buf);  // x, y, z
        for (int iobs = 0; iobs < nobs; iobs++) {
            getline(ifs, buf);
            std::stringstream ss(buf);
            for (int idim = 0; idim < 3; idim++) {
                getline(ss, buf, ' ');
                obs_point_buf[iobs * 3 + idim] = std::stod(buf);
            }
        }
    }

    MPI_Bcast(&nobs, 1, MPI_INT, 0, MPI_COMM_WORLD);
    obs_point_buf.resize(3 * nobs);
    obs_point.resize(nobs, std::vector<double>(3));

    MPI_Bcast(&obs_point_buf[0], obs_point_buf.size(), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
    for (int iobs = 0; iobs < nobs; iobs++) {
        for (int idim = 0; idim < 3; idim++) {
            obs_point[iobs][idim] = obs_point_buf[iobs * 3 + idim];
        }
    }
}
