#include "writer.hpp"

namespace nsp_error_estimator {
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
}
