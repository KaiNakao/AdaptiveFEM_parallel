#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "mesh_refiner.hpp"

int main() {
    std::string data_dir("result/");

    Refiner refiner = Refiner(data_dir);
    refiner.readData();
    refiner.executeRefinement();
    refiner.writeNewMesh();
}
