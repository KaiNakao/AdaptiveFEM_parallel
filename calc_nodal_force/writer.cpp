#include "writer.hpp"
void write_nodal_force(const int &myid, const std::vector<int> &node_id_obs,
                       const std::vector<double> &fvec) {
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "mdata/" + rank_id + ".nodal_force.dat";
    std::ofstream ofs;
    ofs.open(filename);
    if (!ofs) {
        std::cerr << "Error: file open failed: " << filename << std::endl;
        MPI_Finalize();
        std::exit(1);
    }
    int cnt = 0;
    for (int ipoint = 0; ipoint < 8; ipoint++) {
        if (node_id_obs[ipoint] != -1) {
            cnt++;
        }
    }
    ofs << "number of nodes\n";
    ofs << cnt << "\n";
    ofs << "node_id fx fy fz\n";
    for (int ipoint = 0; ipoint < 8; ipoint++) {
        if (node_id_obs[ipoint] == -1) {
            continue;
        }
        ofs << node_id_obs[ipoint] + 1 << " ";
        ofs << std::scientific << std::setprecision(20);
        for (int idim = 0; idim < 3; idim++) {
            ofs << fvec[ipoint * 3 + idim] << " ";
        }
        ofs << "\n";
    }
    if (cnt > 0) {
        std::cout << "nodal force is written to " << filename << std::endl;
    }
}