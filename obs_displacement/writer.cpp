#include "writer.hpp"

void write_eta(const int &myid, const std::vector<double> &eta_arr) {
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "mdata/" + rank_id + ".eta.bin";
    FILE *fp;
    if ((fp = fopen(filename.c_str(), "wb")) == NULL) {
        std::cerr << "cannot open file: " << filename << std::endl;
        exit(1);
    }

    fwrite(&eta_arr[0], sizeof(double), eta_arr.size(), fp);
}

void write_marked_elem(const int &myid, const std::vector<int> &marked_elem) {
    std::string rank_id = std::to_string(myid);
    rank_id.insert(rank_id.begin(), 6 - rank_id.size(), '0');
    std::string filename = "mdata/" + rank_id + ".marked_elem.bin";
    FILE *fp;
    if ((fp = fopen(filename.c_str(), "wb")) == NULL) {
        std::cerr << "cannot open file: " << filename << std::endl;
        exit(1);
    }

    fwrite(&marked_elem[0], sizeof(int), marked_elem.size(), fp);
}

void write_displacement_obs(
    const std::vector<std::vector<double>> &displacement_obs,
    const int &iload) {
    std::vector<double> displacement_obs_buf(displacement_obs.size() * 3);
    for (int iobs = 0; iobs < displacement_obs.size(); iobs++) {
        for (int idim = 0; idim < 3; idim++) {
            displacement_obs_buf[iobs * 3 + idim] =
                displacement_obs[iobs][idim];
        }
    }

    std::string load_id = std::to_string(iload + 1);
    load_id.insert(load_id.begin(), 4 - load_id.size(), '0');
    std::string filename = "displacement/" + load_id + "_obs.bin";
    FILE *fp;
    if ((fp = fopen(filename.c_str(), "wb")) == NULL) {
        std::cerr << "cannot open file: " << filename << std::endl;
        exit(1);
    }
    fwrite(&displacement_obs_buf[0], sizeof(double),
           displacement_obs_buf.size(), fp);
}