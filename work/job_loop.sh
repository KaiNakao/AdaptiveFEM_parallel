#!/bin/bash
#PBS -j oe
#PBS -q calc
##PBS -N hpl
#PBS -l select=4:ncpus=80:mpiprocs=20

#########################################################
export WORKDIR=work1
echo "WORKDIR: ${WORKDIR}"

# number of partitions
export NPART=64
echo "NPART: ${NPART}"
cd ${PBS_O_WORKDIR}

source /home/ap/intel_2021.1/oneapi/setvars.sh
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=1G
export LD_LIBRARY_PATH=../hdf5-1.10.7_lib_intel/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=../metis_lib_intel_64bit/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=../metis-5.1.0_icc_32bit/lib:$LD_LIBRARY_PATH
ulimit -s unlimited

mkdir -p ${WORKDIR}
cd ${WORKDIR}
mkdir -p data hdata cdata mdata 2Doutput surf_mesh displacement result
cd ../
mkdir -p ../${WORKDIR}_tmp
cd ../${WORKDIR}_tmp
mv ../work/${WORKDIR}/* ./

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

for iter in {1..5}
do
    echo "iteration: ${iter}"
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

    # merge local results
    echo "merge local results"
    julia ../merge_local_result/to_AFEM.jl > merge_local_result.log

    # mesh refinement
    echo "mesh refinement"
    cd ../refiner
    make
    cd ../${WORKDIR}_tmp
    ../refiner/main > refiner.log

    mkdir -p result/vtu
    mkdir -p result/fig
    # cd ../refiner
    julia ../refiner/write_mesh.jl
    julia ../refiner/fig_out.jl

    echo "write hdf5 file for new mesh"
    julia ../refiner/write_hdf5.jl > write_hdf5.log

    # set new workdir
    export WORKDIR_ORG=${WORKDIR}
    WORKDIR_NUM=$(echo ${WORKDIR} | grep -o -E '[0-9]+$')
    WORKDIR_NEW=$(echo ${WORKDIR} | sed "s/${WORKDIR_NUM}$/$((${WORKDIR_NUM} + 1))/")
    export WORKDIR=${WORKDIR_NEW}

    echo "NEW WORKDIR: ${WORKDIR}"

    # write hdf5 file for new mesh
    cd ../work
    mkdir -p ${WORKDIR}
    cd ${WORKDIR}
    mkdir -p data hdata cdata mdata 2Doutput surf_mesh displacement result
    mv ../../${WORKDIR_ORG}_tmp/result/tet* ./data/
    cp ../../${WORKDIR_ORG}_tmp/data/log_setting.dat ./data/
    cp ../../${WORKDIR_ORG}_tmp/data/material.dat ./data/
    cp ../../${WORKDIR_ORG}_tmp/data/modeldomain.dat ./data/
    cp ../../${WORKDIR_ORG}_tmp/data/modeling_setting.dat ./data/
    cp ../../${WORKDIR_ORG}_tmp/data/para_setting.dat ./data/
    cp ../../${WORKDIR_ORG}_tmp/data/pointload.dat ./data/

    # reduce ds in modeling_setting.dat to half
    file="data/modeling_setting.dat"
    current_ds=$(awk '/^ds/{getline; print}' $file)
    new_ds=$(echo "$current_ds / 2" | bc -l)
    awk -v new_ds="$new_ds" '
    /^ds/ {print; getline; print new_ds; next} 
    {print}
    ' $file > temp_file && mv temp_file $file

    cd ../
    mv ../${WORKDIR_ORG}_tmp/* ${WORKDIR_ORG}/
    rmdir ../${WORKDIR_ORG}_tmp

    mkdir -p ../${WORKDIR}_tmp
    cd ../${WORKDIR}_tmp
    mv ../work/${WORKDIR}/* ./
done

cd ../work
cd ${WORKDIR}
mv ../../${WORKDIR}_tmp/* ./
rmdir ../../${WORKDIR}_tmp
