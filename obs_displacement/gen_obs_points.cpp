#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

int main() {
    // read model domain
    std::ifstream ifs;
    std::string buf;
    std::string filename;
    filename = "data/modeldomain.dat";
    ifs.open(filename);
    if (!ifs) {
        std::cerr << "Error: file " << filename << " not found." << std::endl;
        std::exit(1);
    }
    getline(ifs, buf);  // xmin-plane
    getline(ifs, buf);
    double xmin = std::stod(buf);
    getline(ifs, buf);  // xmax-plane
    getline(ifs, buf);
    double xmax = std::stod(buf);
    getline(ifs, buf);  // ymin-plane
    getline(ifs, buf);
    double ymin = std::stod(buf);
    getline(ifs, buf);  // ymax-plane
    getline(ifs, buf);
    double ymax = std::stod(buf);
    getline(ifs, buf);  // zmin-plane
    getline(ifs, buf);
    double zmin = std::stod(buf);
    ifs.close();

    double zmax = 0.0;
    filename = "data/sur0001.dat";
    ifs.open(filename);
    if (!ifs) {
        std::cerr << "Error: file " << filename << " not found." << std::endl;
        std::exit(1);
    }
    getline(ifs, buf);  // nx, ny
    getline(ifs, buf);  // nx, ny
    getline(ifs, buf);  // DEM data
    while (getline(ifs, buf)) {
        zmax = std::max(zmax, std::stod(buf));
    }
    ifs.close();
    std::cout << "xmin: " << xmin << std::endl;
    std::cout << "xmax: " << xmax << std::endl;
    std::cout << "ymin: " << ymin << std::endl;
    std::cout << "ymax: " << ymax << std::endl;
    std::cout << "zmin: " << zmin << std::endl;
    std::cout << "zmax: " << zmax << std::endl;

    // generate observation points
    std::vector<double> x_obs = {6e4,  7e4,  8e4,  9e4,  10e4, 11e4, 12e4,
                                 13e4, 14e4, 15e4, 16e4, 17e4, 18e4};
    std::vector<double> y_obs = {5e4,  6e4,  7e4,  8e4,  9e4,  10e4,
                                 11e4, 12e4, 13e4, 14e4, 15e4, 16e4};
    std::vector<double> z_obs = {20e4, 21e4, 22e4, 23e4, 24e4, 25e4, 26e4};

    if (x_obs[0] < xmin || x_obs.back() > xmax) {
        std::cerr << "Error: x_obs is out of model domain." << std::endl;
        std::exit(1);
    }
    if (y_obs[0] < ymin || y_obs.back() > ymax) {
        std::cerr << "Error: y_obs is out of model domain." << std::endl;
        std::exit(1);
    }
    if (z_obs[0] < zmin || z_obs.back() > zmax) {
        std::cerr << "Error: z_obs is out of model domain." << std::endl;
        std::exit(1);
    }

    // write observation points
    filename = "data/obs_points.dat";
    std::ofstream ofs;
    ofs.open(filename);
    if (!ofs) {
        std::cerr << "Error: file " << filename << " could not be opened."
                  << std::endl;
        std::exit(1);
    }
    ofs << "number of observation points" << std::endl;
    ofs << x_obs.size() * y_obs.size() * z_obs.size() << std::endl;
    ofs << "x, y, z" << std::endl;
    for (int i = 0; i < x_obs.size(); i++) {
        for (int j = 0; j < y_obs.size(); j++) {
            for (int k = 0; k < z_obs.size(); k++) {
                ofs << x_obs[i] << " " << y_obs[j] << " " << z_obs[k]
                    << std::endl;
            }
        }
    }
    ofs.close();
}