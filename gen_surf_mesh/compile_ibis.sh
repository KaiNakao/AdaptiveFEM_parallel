CXX="mpiicpc -cxx=icpx"
CC=icc

HDFPATH=" ../hdf5-1.10.7_lib_intel/"
INCLHDF="-I${HDFPATH}/include"
LINKHDF="-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib"

OPTLEV="-O3 -fopenmp"
OPTCPP="-std=c++20"
OPTC=" "

$CC hdf_lib.c -c $INCLHDF $OPTC $OPTLEV
$CXX main.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -o main