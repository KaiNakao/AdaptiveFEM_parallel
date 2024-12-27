# Compile on Kerria
# use 64 bit integers by default

##### DEBUG
OPTLEV="-g -DOUTPUT_LOG"
OPTF="-CcdII8"

##### RELEASE
OPTLEV="-Kfast -DOUTPUT_LOG"
OPTLEV="-Kfast"
#OPTLEV="-Kfast"
OPTF="-CcdII8"

##### path to libraries
HDFPATH="/opt/aics/netcdf/k-serial-noszip"
HDFPATH=${HOME}/local

OPTC=" "
OPTCPP=" "
OMP="-Kopenmp"
INCLHDF="-I${HDFPATH}/include"
LINKHDF="-L${HDFPATH} ${HDFPATH}/lib/libhdf5_hl.a ${HDFPATH}/lib/libhdf5.a -lz -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib"
# for linking fortran subroutines to c++ program
LINKFORT="--linkfortran"

##### hdf5 calling routine
mpifccpx hdf_lib.c -c $INCLHDF $OPTC $OPTLEV $OMP

##### prepare model (usually done on K/FX10)
mpiFCCpx mpi_prepare_hybrid_model_refine.cpp $INCLHDF $OMP $OPTLEV $OPTCPP -c
mpifrtpx MPI_mapping.f MPI_ABC2_mod_lib.f MPI_ABC4_mod.F 4to10_grid.F 4to10_grid_merge4and10.F read_tet10_geometry_hdf.f 4to10_confirming_mesh.F MPI_CRS_connect_block4.F MPI_CRS_connect_block.F MPI_CRS_nodes4.F MPI_CRS_nodes.F -c $OPTF $OPTLEV $OMP
mpiFCCpx read_tet10_geometry_hdf.o MPI_ABC2_mod_lib.o MPI_ABC4_mod.o MPI_mapping.o 4to10_grid.o 4to10_grid_merge4and10.o mpi_prepare_hybrid_model_refine.o 4to10_confirming_mesh.o MPI_CRS_connect_block4.o MPI_CRS_connect_block.o MPI_CRS_nodes4.o MPI_CRS_nodes.o $LINKFORT $INCLHDF $LINKHDF $OMP $OPTLEV $OPTCPP -o ./mpi_prepare_hybrid_model_refine.exe

