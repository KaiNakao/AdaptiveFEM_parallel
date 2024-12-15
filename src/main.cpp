#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "inc/mesh_refiner.hpp"

int main() 
{
    std::string data_dir("data/analysis_result_iburi_large/");

    Refiner refiner = Refiner(data_dir);
    refiner.executeRefinement();
}
