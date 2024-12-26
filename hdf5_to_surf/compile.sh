#!/bin/sh
# OMP=-fopenmp
HDFPATH="../hdf5-1.10.7_lib_intel/"
INCLHDF="-I${HDFPATH}/include/"
LINKHDF="-L${HDFPATH}/lib/ -lhdf5 -lhdf5_hl"
FLAG="-g -traceback -CB"
FLAG="-g -traceback -CB -integer-size 64 -mcmodel=large -heap-arrays -zero"
FLAG="-O3 -integer-size 64 -mcmodel=large -heap-arrays -zero"
#FLAG="-O3"
SRC=hdf5_to_surf
icc -c $INCLHDF -qopenmp  hdf_lib.c
ifort -c $FLAG ${SRC}.F90
ifort -qopenmp $FLAG  $LINKHDF hdf_lib.o ${SRC}.o  -o ${SRC}.exe
