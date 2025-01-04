#!/bin/bash
#PBS -j oe
#PBS -q calc
##PBS -N hpl
#PBS -l select=1:ncpus=80:mpiprocs=20

#########################################################
export WORKDIR=work3
echo "WORKDIR: ${WORKDIR}"

# number of partitions
export NPART=16
echo "NPART: ${NPART}"
cd ${PBS_O_WORKDIR}

source /home/ap/intel_2021.1/oneapi/setvars.sh
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
export LD_LIBRARY_PATH=../hdf5-1.10.7_lib_intel/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=../metis_lib_intel_64bit/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=../metis-5.1.0_icc_32bit/lib:$LD_LIBRARY_PATH
ulimit -s unlimited

mkdir -p ../${WORKDIR}_tmp
cd ../${WORKDIR}_tmp
mv ../work/${WORKDIR}/* ./

# generate second order mesh
echo "generate second order mesh"
cd ../hdata_to_cdatamdata
./compile_modeling_x86-64.sh 
cd ../${WORKDIR}_tmp
mpiexec -n ${NPART} ../hdata_to_cdatamdata/mpi_prepare_hybrid_model_refine.exe > hdata_to_cdatamdata.log

# generate surface mesh
echo "generate surface mesh"
cd ../gen_surf_mesh
./compile_ibis.sh
cd ../${WORKDIR}_tmp
mpiexec -n ${NPART} ../gen_surf_mesh/main > gen_surf_mesh.log

# finite element solver
export OMP_NUM_THREADS=4
echo "finite element solver"
cd ../src_analysis_pointload
make
cd ../${WORKDIR}_tmp
mpiexec -n ${NPART} ../src_analysis_pointload/analysis.exe.npc1 > analysis.log

# error analysis
echo "error analysis"
cd ../error_estimator
./compile_ibis.sh
cd ../${WORKDIR}_tmp
mpiexec -n ${NPART} ../error_estimator/main > error_estimator.log

cd ../work
cd ${WORKDIR}
mv ../../${WORKDIR}_tmp/* ./
rmdir ../../${WORKDIR}_tmp
