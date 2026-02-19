#include "surf.hpp"

namespace nsp_obs_displacement {
void calc_projection(const int &myid,
                     const std::vector<std::vector<double>> &coor,
                     std::vector<std::vector<double>> &obs_point) {
    std::vector<std::vector<double>> load_arr;
    int nload;
    if (myid == 0) {
        read_loads(nload, load_arr);
    }
    MPI_Bcast(&nload, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (myid != 0) {
        load_arr.resize(nload, std::vector<double>(5));
    }
    for (int iload = 0; iload < nload; iload++) {
        MPI_Bcast(&load_arr[iload][0], 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    obs_point.resize(nload, std::vector<double>(3));

    std::vector<std::vector<int>> surf_cny;
    read_surf_cny(myid, surf_cny);
    if (myid == 3) {
        std::cout << "surf_cny.size(): " << surf_cny.size() << std::endl;
        for (int ielem = 0; ielem < surf_cny.size(); ielem++) {
            for (int inode = 0; inode < 4; inode++) {
                std::cout << surf_cny[ielem][inode] << " ";
            }
            std::cout << std::endl;
        }
    }

    std::vector<double> load_point(2), load_point_proj(3), load_point_proj_g(3),
        dx(2);
    std::vector<std::vector<double>> xnode(3, std::vector<double>(3));
    std::vector<std::vector<double>> dxdr(2, std::vector<double>(2));
    std::vector<std::vector<double>> drdx(2, std::vector<double>(2));
    int found = 0;
    for (int iload = 0; iload < nload; iload++) {
        found = 0;
        load_point[0] = load_arr[iload][0];
        load_point[1] = load_arr[iload][1];

        for (int ielem = 0; ielem < surf_cny.size(); ielem++) {
            const auto &nodes = surf_cny[ielem];
            for (int inode = 0; inode < 3; inode++) {
                for (int idim = 0; idim < 3; idim++) {
                    xnode[inode][idim] = coor[nodes[inode]][idim];
                }
            }
            dxdr[0][0] = xnode[1][0] - xnode[0][0];
            dxdr[0][1] = xnode[2][0] - xnode[0][0];
            dxdr[1][0] = xnode[1][1] - xnode[0][1];
            dxdr[1][1] = xnode[2][1] - xnode[0][1];
            double detj = dxdr[0][0] * dxdr[1][1] - dxdr[0][1] * dxdr[1][0];
            drdx[0][0] = dxdr[1][1] / detj;
            drdx[0][1] = -dxdr[0][1] / detj;
            drdx[1][0] = -dxdr[1][0] / detj;
            drdx[1][1] = dxdr[0][0] / detj;

            dx[0] = load_point[0] - xnode[0][0];
            dx[1] = load_point[1] - xnode[0][1];

            double r1 = drdx[0][0] * dx[0] + drdx[0][1] * dx[1];
            double r2 = drdx[1][0] * dx[0] + drdx[1][1] * dx[1];

            if (r1 >= -1e-8 && r2 >= -1e-8 && r1 + r2 <= 1.0 + 1e-8) {
                load_point_proj[0] = (1.0 - r1 - r2) * xnode[0][0] +
                                     r1 * xnode[1][0] + r2 * xnode[2][0];
                load_point_proj[1] = (1.0 - r1 - r2) * xnode[0][1] +
                                     r1 * xnode[1][1] + r2 * xnode[2][1];
                load_point_proj[2] = (1.0 - r1 - r2) * xnode[0][2] +
                                     r1 * xnode[1][2] + r2 * xnode[2][2];
                found = 1;
                break;
            }
        }
        for (int idim = 0; idim < 3; idim++) {
            load_point_proj[idim] *= found;
        }
        MPI_Allreduce(&load_point_proj[0], &load_point_proj_g[0], 3, MPI_DOUBLE,
                      MPI_SUM, MPI_COMM_WORLD);
        if (myid == 0) {
            std::cout << std::scientific << std::setprecision(20);
            std::cout << "projection: " << load_point_proj_g[0] << " "
                      << load_point_proj_g[1] << " " << load_point_proj_g[2]
                      << std::endl;
        }
        obs_point[iload] = load_point_proj_g;
    }
}

void read_loads(int &nload, std::vector<std::vector<double>> &load_arr) {
    std::ifstream ifs;
    std::string buf;
    ifs.open("data/pointload.dat");
    if (!ifs) {
        std::cerr << "pointload.dat is not found." << std::endl;
        exit(1);
    }
    std::getline(ifs, buf);  // number of loads
    std::getline(ifs, buf);  // nload
    nload = std::stoi(buf);
    std::getline(ifs, buf);  // x y ex ey ez
    load_arr.resize(nload, std::vector<double>(5));
    for (int iload = 0; iload < nload; iload++) {
        std::getline(ifs, buf);
        std::istringstream iss(buf);
        iss >> load_arr[iload][0] >> load_arr[iload][1] >> load_arr[iload][2] >>
            load_arr[iload][3] >> load_arr[iload][4];
    }
}

void read_surf_cny(const int &myid, std::vector<std::vector<int>> &surf_cny) {
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "surf_mesh/" + rank_id + "_nelem.dat";
    std::ifstream ifs;
    std::string buf;
    ifs.open(filename);
    if (!ifs) {
        std::cerr << filename << " is not found." << std::endl;
        exit(1);
    }
    std::getline(ifs, buf);
    int surf_nelem = std::stoi(buf);

    if (surf_nelem > 0) {
        filename = "surf_mesh/" + rank_id + "_cny.bin";
        FILE *fp;
        if ((fp = fopen(filename.c_str(), "rb")) == NULL) {
            std::cerr << filename << " is not found." << std::endl;
            exit(1);
        }
        std::vector<int> surf_cny_buf(4 * surf_nelem);
        fread(&surf_cny_buf[0], sizeof(int), 4 * surf_nelem, fp);
        fclose(fp);

        surf_cny.resize(surf_nelem, std::vector<int>(4));
        for (int ielem = 0; ielem < surf_nelem; ielem++) {
            for (int inode = 0; inode < 4; inode++) {
                surf_cny[ielem][inode] = surf_cny_buf[4 * ielem + inode] - 1;
            }
        }
    } else {
        surf_cny.resize(0);
    }
}
}
