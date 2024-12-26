#!/bin/bash
#PBS -j oe
#PBS -q calc
##PBS -N hpl
#PBS -l select=1:ncpus=80:mpiprocs=20

#########################################################
export WORKDIR=work1
echo "WORKDIR: ${WORKDIR}"

# number of partitions
export NPART=16
echo "NPART: ${NPART}"
cd ${PBS_O_WORKDIR}

source /home/ap/intel_2021.1/oneapi/setvars.sh
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../hdf5-1.10.7_lib_intel/lib
ulimit -s unlimited

mkdir -p ../tmp
cd ../tmp
mv ../work/${WORKDIR}/* ./

# generate second order mesh
echo "generate second order mesh"
cd ../hdata_to_cdatamdata
./compile_modeling_x86-64.sh 
cd ../tmp
mpiexec -n ${NPART} ../hdata_to_cdatamdata/mpi_prepare_hybrid_model_refine.exe

# generate surface mesh
echo "generate surface mesh"
cd ../gen_surf_mesh
./compile_ibis.sh
cd ../tmp
mpiexec -n ${NPART} ../gen_surf_mesh/main

# finite element solver
export OMP_NUM_THREADS=4
echo "finite element solver"
cd ../src_analysis_pointload
make
cd ../tmp
mpiexec -n ${NPART} ../src_analysis_pointload/analysis.exe.npc1


cd ../work
cd ${WORKDIR}
mv ../../tmp/* ./
rmdir ../../tmp
