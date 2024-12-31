#!/bin/bash
#PBS -j oe
#PBS -q calc
##PBS -N hpl
#PBS -l select=1:ncpus=80:mpiprocs=80

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
ulimit -s unlimited

mkdir -p ../tmp
cd ../tmp
mv ../work/${WORKDIR}/* ./

# merge local results
echo "merge local results"
julia ../merge_local_result/to_AFEM.jl > merge_local_result.log

# mesh refinement
echo "mesh refinement"
cd ../refiner
make
cd ../tmp
../refiner/main > refiner.log

mkdir -p result/vtu
mkdir -p result/fig
# cd ../refiner
julia ../refiner/write_mesh.jl
julia ../refiner/fig_out.jl

cd ../work
cd ${WORKDIR}
mv -f ../../tmp/* ./
rm -rf   ../../tmp

