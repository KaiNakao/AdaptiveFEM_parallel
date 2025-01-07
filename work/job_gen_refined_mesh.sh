#!/bin/bash
#PBS -j oe
#PBS -q calc
##PBS -N hpl
#PBS -l select=4:ncpus=80:mpiprocs=80

#########################################################
export WORKDIR_ORG=work3
export WORKDIR_NEW=work4
echo "WORKDIR_ORG: ${WORKDIR_ORG}"
echo "WORKDIR_NEW: ${WORKDIR_NEW}"

# number of partitions
export NPART=16
echo "NPART: ${NPART}"
cd ${PBS_O_WORKDIR}

source /home/ap/intel_2021.1/oneapi/setvars.sh
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../hdf5-1.10.7_lib_intel/lib
ulimit -s unlimited

mkdir -p ../tmp_org
cd ../tmp_org
mv ../work/${WORKDIR_ORG}/* ./

mkdir -p ../tmp_new
cd ../tmp_new
mkdir -p data hdata cdata mdata 2Doutput surf_mesh displacement result

# write hdf5 file for new mesh
echo "write hdf5 file for new mesh"
cd ../tmp_org
julia ../refiner/write_hdf5.jl > write_hdf5.log
mv result/tet4* ../tmp_new/data/
cp data/log_setting.dat ../tmp_new/data/
cp data/material.dat ../tmp_new/data/ 
cp data/modeldomain.dat ../tmp_new/data/
cp data/modeling_setting.dat ../tmp_new/data/
cp data/para_setting.dat ../tmp_new/data/
cp data/pointload.dat ../tmp_new/data/
cp result/refinement_edge.bin ../tmp_new/result/
cp result/marked_edge.bin ../tmp_new/result/

# reduce ds in modeling_setting.dat to half
cd ../tmp_new
file="data/modeling_setting.dat"
current_ds=$(awk '/^ds/{getline; print}' $file)
new_ds=$(echo "$current_ds / 2" | bc -l)
awk -v new_ds="$new_ds" '
/^ds/ {print; getline; print new_ds; next} 
{print}
' $file > temp_file && mv temp_file $file

# check surface
echo "check surface"
cd ../hdf5_to_surf
./compile.sh
cd ../tmp_new
../hdf5_to_surf/hdf5_to_surf.exe -s > hdf5_to_surf.log
# check that only data/tri_surf.001.vtk was generated. Otherwise, mesh generation failed.

# partition mesh
echo "partition mesh"
cd ../hdf5_to_hdata/
./compile_partition_2part_x86-64.sh
cd ../tmp_new
export OMP_NUM_THREADS=1
../hdf5_to_hdata/partition_hybrid_model_2part.exe > hdf5_to_hdata.log

cd ../work
mkdir -p ${WORKDIR_ORG}
cd ${WORKDIR_ORG}
mv ../../tmp_org/* ./
rmdir ../../tmp_org

cd ../
mkdir -p ${WORKDIR_NEW}
cd ${WORKDIR_NEW}
mv ../../tmp_new/* ./
rmdir ../../tmp_new
