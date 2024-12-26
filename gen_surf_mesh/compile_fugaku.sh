HDFPATH=${HOME}/local
INCLHDF="-I${HDFPATH}/include"
LINKHDF="-L${HDFPATH} ${HDFPATH}/lib/libhdf5_hl.a ${HDFPATH}/lib/libhdf5.a -lz -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib"

OPTLEV="-Kfast"
OPTC=" "
OPTCPP=" "
OMP="-Kopenmp"

mpifccpx hdf_lib.c -c $INCLHDF $OPTC $OPTLEV $OMP
mpiFCCpx main.cpp $INCLHDF $LINKHDF $OMP $OPTLEV $OPTCPP -o main