# Compile on x86-64
# use 64 bit integers by default

#CXX=g++
#CC=gcc
#FC=gfortran
#HDFPATH="/home/fujita/Softwares/hdf5-1.8.13_gcc/hdf5"
#HDFPATH="/home/fujita/Softwares/hdf5-1.8.13/hdf5"
#METISPATH="/home/fujita/Softwares/bin"

##### path to libraries
CXX=icpc
CC=icc
FC=ifort
#HDFPATH="/home/fujita/Softwares/hdf5-1.10.7/hdf5/"
#METISPATH="/home/fujita/Softwares/metis-5.1.0_icc_64bit/metis-5.1.0_icc_64bit/"
METISPATH="../metis_lib_intel_64bit/"
HDFPATH=" ../hdf5-1.10.7_lib_intel/"

##### DEBUG
OPTLEV="-g -fopenmp -DOUTPUT_LOG"
OPTF="-integer-size 64 -mcmodel=large -i-dynamic -heap-arrays -zero -check -warn all"

##### RELEASE
OPTLEV="-O3 -fopenmp"
OPTF="-integer-size 64 -mcmodel=large -heap-arrays -zero"

OPTC=" "
OPTCPP="-std=c++11"
INCLHDF="-I${HDFPATH}/include"
LINKHDF="-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib"
INCLMETIS="-I${METISPATH}/include"
LINKMETIS="-L${METISPATH}/lib/ -lmetis"
# for linking fortran subroutines to c++ program
LINKFORT="-limf -lifcore -lifport -lifcoremt"

##### hdf5 calling routine
#$CC hdf_lib.c -c $INCLHDF $OPTC $OPTLEV

##### modeling settings
#$FC geometry_setting.f $OPTF $OPTLEV $OMP -o ./geometry_setting.exe
#$FC modelingparasetting.f $OPTF $OPTLEV $OMP -o ./modelingparasetting.exe

##### modeling
# modify NSQUASH for changing height to width ratio of mesh (DEFAULT is -DNSQUASH=1)
#$FC modeling_main.F delauny.F OMP_modeling_sub1.F -c $OPTF $OPTLEV -DNSQUASH=1
#$FC modeling_main.o delauny.o OMP_modeling_sub1.o hdf_lib.o $OPTLEV $LINKHDF -o ./modeling_hybrid_hdf.exe

##### modify mesh
#$CXX modify_mesh.cpp $INCLHDF $OPTLEV $LINKHDF -limf -lifcore -o modify_mesh_infinite.exe

##### hdf5 calling routine
icc -lmpi hdf_lib.c -c $INCLHDF $OPTC $OPTLEV $OMP

##### prepare model (usually done on K/FX10)
icpc -lmpi mpi_prepare_hybrid_model_refine.cpp $INCLHDF $OMP $OPTLEV $OPTCPP -c
ifort -lmpi MPI_mapping.f MPI_ABC2_mod_lib.f MPI_ABC4_mod.F 4to10_grid.F 4to10_grid_merge4and10.F read_tet10_geometry_hdf.f 4to10_confirming_mesh.F MPI_CRS_connect_block4.F MPI_CRS_connect_block.F MPI_CRS_nodes4.F MPI_CRS_nodes.F -c $OPTF $OPTLEV $OMP
icpc -lmpi read_tet10_geometry_hdf.o MPI_ABC2_mod_lib.o MPI_ABC4_mod.o MPI_mapping.o 4to10_grid.o 4to10_grid_merge4and10.o mpi_prepare_hybrid_model_refine.o 4to10_confirming_mesh.o MPI_CRS_connect_block4.o MPI_CRS_connect_block.o MPI_CRS_nodes4.o MPI_CRS_nodes.o $LINKFORT $INCLHDF $LINKHDF $OMP $OPTLEV $OPTCPP -o ./mpi_prepare_hybrid_model_refine.exe

