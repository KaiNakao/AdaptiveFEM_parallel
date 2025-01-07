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
    std::vector<double> x_obs;
    std::vector<double> y_obs;
    std::vector<double> z_obs;

    double x = (xmin + xmax) / 2.0 - 60000.0;
    x = x - fmod(x, 10000.0);
    while (x < (xmin + xmax) / 2.0 + 60000.0) {
        x_obs.push_back(x);
        x += 10000.0;
    }
    double y = (ymin + ymax) / 2.0 - 60000.0;
    y = y - fmod(y, 10000.0);
    while (y < (ymin + ymax) / 2.0 + 60000.0) {
        y_obs.push_back(y);
        y += 10000.0;
    }
    double z = zmax - 60000.0;
    z = z - fmod(z, 10000.0);
    while (z < zmax) {
        z_obs.push_back(z);
        z += 10000.0;
    }
    std::cout << "x_obs.size(): " << x_obs.size() << std::endl;
    for (int i = 0; i < x_obs.size(); i++) {
        std::cout << x_obs[i] << std::endl;
    }
    std::cout << "y_obs.size(): " << y_obs.size() << std::endl;
    for (int i = 0; i < y_obs.size(); i++) {
        std::cout << y_obs[i] << std::endl;
    }
    std::cout << "z_obs.size(): " << z_obs.size() << std::endl;
    for (int i = 0; i < z_obs.size(); i++) {
        std::cout << z_obs[i] << std::endl;
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