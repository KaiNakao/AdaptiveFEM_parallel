#include <mpi.h>
#include <iostream>
#include <string>

#include "analysis_pointload/pointload_solver.hpp"
#include "error_estimator/error_estimator.hpp"
#include "obs_displacement/obs_displacement.hpp"

int main(int argc, char **argv) {
    int myid, numprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    bool write_displacement = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--write-displacement") {
            write_displacement = true;
        } else if (arg == "--no-write-displacement") {
            write_displacement = false;
        }
    }

    Pointload_Solver pointload_solver(myid, numprocs);
    pointload_solver.prep();
    MPI_Barrier(MPI_COMM_WORLD);

    Error_Estimator error_estimator(myid, numprocs);
    error_estimator.prep();
    MPI_Barrier(MPI_COMM_WORLD);

    Obs_Displacement obs_displacement(myid, numprocs);
    obs_displacement.prep();
    MPI_Barrier(MPI_COMM_WORLD);

    int nload = pointload_solver.nload;
    // int nload = 4; // for testing
    std::vector<double> uvec;
    MPI_Barrier(MPI_COMM_WORLD);

    for (int iload = 0; iload < nload; iload++) {
        if (myid == 0) {
            std::cout << "iload:" << iload << std::endl;
        }
        pointload_solver.exec(iload, write_displacement);
        MPI_Barrier(MPI_COMM_WORLD);

        uvec = pointload_solver.up;
        error_estimator.exec(iload, uvec);
        MPI_Barrier(MPI_COMM_WORLD);

        obs_displacement.exec(iload, uvec);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    error_estimator.load_elem = pointload_solver.load_elem_arr;
    error_estimator.post();

    MPI_Finalize();
    return 0;
}
