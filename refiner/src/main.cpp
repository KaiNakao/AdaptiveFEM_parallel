#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "inc/mesh_refiner.hpp"
#include "inc/posterior_error.hpp"

int main() {
    std::string data_dir("result/");
    ErrorEstimator error_estimator(data_dir);

    Refiner refiner = Refiner(data_dir);
    refiner.executeRefinement();
}
