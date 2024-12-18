#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "inc/mesh_refiner.hpp"
#include "inc/posterior_error.hpp"

int main() 
{
    std::string data_dir("/data6/itou/AFEM/data/analysis_result_iburi_large/");
    ErrorEstimator error_estimator(data_dir);


    Refiner refiner = Refiner(data_dir);
    refiner.executeRefinement();
}
