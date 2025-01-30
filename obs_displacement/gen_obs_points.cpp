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
    // std::vector<double> x_obs = {6e4,  7e4,  8e4,  9e4,  10e4, 11e4, 12e4,
    //                              13e4, 14e4, 15e4, 16e4, 17e4, 18e4};
    // std::vector<double> y_obs = {5e4,  6e4,  7e4,  8e4,  9e4,  10e4,
    //                              11e4, 12e4, 13e4, 14e4, 15e4, 16e4};
    // std::vector<double> z_obs = {12e4, 13e4, 14e4, 15e4, 16e4, 17e4};
    std::vector<double> x_obs = {6e4, 8e4, 10e4, 12e4, 14e4, 16e4, 18e4};
    std::vector<double> y_obs = {6e4, 8e4, 10e4, 12e4, 14e4};
    std::vector<double> z_obs = {4e4, 6e4, 8e4, 10e4, 12e4, 14e4, 16e4, 18e4, 19e4, 19.5e4};

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

    std::vector<std::vector<double>> obs_points;
    for (int i = 0; i < x_obs.size(); i++) {
        for (int j = 0; j < y_obs.size(); j++) {
            for (int k = 0; k < z_obs.size(); k++) {
                obs_points.push_back({x_obs[i], y_obs[j], z_obs[k]});
            }
        }
    }

    // centroid
    std::vector<double> centroid(3, 0.0);
    filename = "data/target_centroid.dat";
    ifs.open(filename);
    if (!ifs) {
        std::cerr << "Error: file " << filename << " not found." << std::endl;
        std::exit(1);
    }
    getline(ifs, buf);  // centroid coordinate
    for (int idim = 0; idim < 3; idim++) {
        getline(ifs, buf);
        centroid[idim] = std::stod(buf);
    }
    ifs.close();

    double x1, x2, y1, y2, z1, z2;
    x1 = int(centroid[0] / 5e3) * 5e3;
    x2 = x1 + 5e3;
    y1 = int(centroid[1] / 5e3) * 5e3;
    y2 = y1 + 5e3;
    z1 = int(centroid[2] / 5e3) * 5e3;
    z2 = z1 + 5e3;

    obs_points.push_back({x1, y1, z1});
    obs_points.push_back({x2, y1, z1});
    obs_points.push_back({x2, y2, z1});
    obs_points.push_back({x1, y2, z1});
    obs_points.push_back({x1, y1, z2});
    obs_points.push_back({x2, y1, z2});
    obs_points.push_back({x2, y2, z2});
    obs_points.push_back({x1, y2, z2});

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
    ofs << obs_points.size() << std::endl;
    ofs << "x, y, z" << std::endl;
    for (int iobs = 0; iobs < obs_points.size(); iobs++) {
        ofs << obs_points[iobs][0] << " " << obs_points[iobs][1] << " "
            << obs_points[iobs][2] << std::endl;
    }
    ofs.close();
}