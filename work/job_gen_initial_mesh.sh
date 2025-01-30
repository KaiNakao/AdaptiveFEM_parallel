#!/bin/bash
#PBS -j oe
#PBS -q calc
##PBS -N hpl
#PBS -l select=1:ncpus=80:mpiprocs=1

#########################################################
# check DNSQUASH in compile_modeling_x86-64_ibis.sh #
export WORKDIR=work_5600_crust_tmp
echo "WORKDIR: ${WORKDIR}"

# number of partitions
export NPART=320
echo "NPART: ${NPART}"
cd ${PBS_O_WORKDIR}

source /home/ap/intel_2021.1/oneapi/setvars.sh
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../hdf5-1.10.7_lib_intel/lib
ulimit -s unlimited

mkdir -p ../${WORKDIR}_tmp
cd ../${WORKDIR}_tmp

mkdir -p data hdata cdata mdata 2Doutput surf_mesh displacement result

# coordinate transformation of DEM data
echo "coordinate transformation"
python3 ../coord_trans/coord_trans.py > coord_trans.log

# generate mesh
echo "generate mesh"
cd ../dem_to_hdf5
./compile_modeling_x86-64_ibis.sh
cd ../${WORKDIR}_tmp
export OMP_NUM_THREADS=80
../dem_to_hdf5/modeling_hybrid_hdf.exe > dem_to_hdf5.log

# check surface
echo "check surface"
cd ../hdf5_to_surf
./compile.sh
cd ../${WORKDIR}_tmp
../hdf5_to_surf/hdf5_to_surf.exe -s > hdf5_to_surf.log
# check that only data/tri_surf.001.vtk was generated. Otherwise, mesh generation failed.

# partition mesh
echo "partition mesh"
echo -e "num of MPI process\n ${NPART}\nnum of OpenMP threads\n 12" > data/para_setting.dat
cd ../hdf5_to_hdata/
./compile_partition_2part_x86-64.sh 
cd ../${WORKDIR}_tmp
export OMP_NUM_THREADS=1
../hdf5_to_hdata/partition_hybrid_model_2part.exe > hdf5_to_hdata.log

# generate observation points
echo "generate observation points"
cd ../obs_displacement
icpx -O3 gen_obs_points.cpp -o gen_obs_points
cd ../${WORKDIR}_tmp
../obs_displacement/gen_obs_points > gen_obs_points.log

cd ../work
mkdir -p ${WORKDIR}
cd ${WORKDIR}
mv ../../${WORKDIR}_tmp/* ./
rmdir ../../${WORKDIR}_tmp
