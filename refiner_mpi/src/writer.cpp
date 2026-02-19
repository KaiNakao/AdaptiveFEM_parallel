#include "writer.hpp"

void write_shape(const int &myid, const int &nnode, const int &nelem, const int &nmaterial) {
    long fid = 0, count, size;
    // write rdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "rdata/" + rank_id + ".data.h5";

    // check if file exists to open / create
    hdf_create_file_(filename.c_str(), &fid);
    // if (FILE *fp = fopen(filename.c_str(), "r")) {
    //     fclose(fp);
    //     hdf_open_file_(filename.c_str(), &fid);
    // } else {
    //     hdf_create_file_(filename.c_str(), &fid);
    // }

    int settings[3] = {nnode, nelem, nmaterial};
    size = 3;
    hdf_write_int_array_(&fid, "/setting", settings, &size);

    hdf_close_file_(&fid);
}

void write_mesh(const int &myid, 
                const std::vector<std::vector<int>> &cny,
                const std::vector<std::vector<double>> &coor,
                const std::vector<int> &matid_arr) {
    long fid = 0, count, size;
    // write rdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "rdata/" + rank_id + ".data.h5";

    hdf_open_file_(filename.c_str(), &fid);                
    int nelem = cny.size();
    size = 5 * nelem;
    int conn_buf[size];
    for (int ielem = 0; ielem < nelem; ielem++) {
        for (int inode = 0; inode < 4; inode++) {
            conn_buf[ielem * 5 + inode] = cny.at(ielem).at(inode) + 1;
        }
        conn_buf[ielem * 5 + 4] = matid_arr.at(ielem) + 1;
    }
    hdf_write_int_array_(&fid, "/conn", conn_buf, &size);

    int nnode = coor.size();
    size = 3 * nnode;
    double coor_buf[size];
    for (int inode = 0; inode < nnode; inode++) {
        for (int idim = 0; idim < 3; idim++) {
            coor_buf[inode * 3 + idim] = coor[inode][idim];
        }
    }
    hdf_write_double_array_(&fid, "/coor", coor_buf, &size);

    hdf_close_file_(&fid);
}

void write_neighbor(const int &myid, const std::map<int, std::vector<int>> &neighbor_nodes) {
    long fid = 0, count, size, start;
    // write rdata
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "rdata/" + rank_id + ".data.h5";

    hdf_open_file_(filename.c_str(), &fid);
    std::vector<int> mpinode_buf;
    mpinode_buf.push_back(neighbor_nodes.size());
    for (const auto &p : neighbor_nodes) {
        int neighbor_id = p.first;
        mpinode_buf.push_back(neighbor_id);
        mpinode_buf.push_back(p.second.size());
    }
    for (const auto &p : neighbor_nodes) {
        for (int node_id : p.second) {
            mpinode_buf.push_back(node_id + 1); 
        }
    }
    size = mpinode_buf.size();
    hdf_write_int_array_(&fid, "/MPInode", mpinode_buf.data(), &size);
    hdf_close_file_(&fid);
}
