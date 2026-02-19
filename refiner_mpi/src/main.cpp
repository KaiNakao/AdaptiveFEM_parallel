#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <mpi.h>

#include "mesh_refiner.hpp"

int main(int argc, char **argv) {
    int myid, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    Refiner refiner = Refiner(myid, numprocs);
    refiner.readData();
    refiner.executeRefinement();
    refiner.writeData();

    MPI_Finalize();
}
