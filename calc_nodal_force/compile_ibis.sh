CXX="mpiicpc -cxx=icpx"
CC=icc
FC=ifx

HDFPATH=" ../hdf5-1.10.7_lib_intel/"
INCLHDF="-I${HDFPATH}/include"
LINKHDF="-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib"
LINKFORTRAN="-lifcore"

# OPTLEV="-O3 -fopenmp"
OPTLEV="-g -fopenmp"
OPTCPP="-std=c++20"
OPTC=" "
OPTF=" "

$CXX reader.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
$CXX find_nodes.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
$CXX writer.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
# $CXX extend_mesh.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
# $CXX posterior_error.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
# $CXX mark_elem.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
# $CXX writer.cpp $INCLHDF $LINKHDF $OPTLEV $OPTCPP -c
$FC sub.F90 $INCLHDF $LINKHDF $OPTLEV $OPTF -c
$CXX main.cpp sub.o reader.o find_nodes.o writer.o $INCLHDF $LINKHDF $LINKFORTRAN $OPTLEV $OPTCPP -o main 