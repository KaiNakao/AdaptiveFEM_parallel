CXX="mpiicpc -cxx=icpx"
CC=icc

HDFPATH=" ../hdf5-1.10.7_lib_intel/"
INCLHDF="-I${HDFPATH}/include"
LINKHDF="-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib"

OPTLEV="-O3 -fopenmp"
# OPTLEV="-g -fopenmp"
OPTCPP="-std=c++20"
OPTC=" "

$CXX reader.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
$CXX extend_mesh.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
$CXX posterior_error.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
$CXX mark_elem.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
$CXX writer.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
$CXX main.cpp reader.o extend_mesh.o posterior_error.o mark_elem.o writer.o $INCLHDF $LINKHDF $OPTLEV $OPTCPP -o main 