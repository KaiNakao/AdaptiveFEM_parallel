#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "inc/mesh_refiner.hpp"
#include "inc/posterior_error.hpp"

int main() {
    std::string data_dir("result/");

    Refiner refiner = Refiner(data_dir);
    refiner.executeRefinement();
}
