# Compile on x86-64
# use 64 bit integers by default

#CXX=g++
#CC=gcc
#FC=gfortran
#HDFPATH="/home/fujita/Softwares/hdf5-1.8.13_gcc/hdf5"
#HDFPATH="/home/fujita/Softwares/hdf5-1.8.13/hdf5"
#METISPATH="/home/fujita/Softwares/bin"

##### path to libraries
CXX=icc
CC=icc
FC=ifort
METISPATH="../metis_lib_intel_64bit/"
HDFPATH=" ../hdf5-1.10.7_lib_intel/"

##### DEBUG
#OPTLEV="-g -fopenmp -DOUTPUT_LOG"
#OPTF="-integer-size 64 -mcmodel=large -i-dynamic -heap-arrays -zero -check -warn all"
OPTLEV="-g -fopenmp"
OPTF="-fdefault-integer-8"

##### RELEASE
#OPTLEV="-O3 -fopenmp"
OPTF="-integer-size 64 -mcmodel=large -heap-arrays -zero"
OPTLEV="-O3 -qopenmp"
# OPTLEV="-g -traceback -O3 -qopenmp -heap-arrays"
OPTF="-i8"

OPTC=" "
OPTCPP=" "
INCLHDF="-I${HDFPATH}/include"
LINKHDF="-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib"
#LINKHDF="-L${HDFPATH} ${HDFPATH}/lib/libhdf5_hl.a ${HDFPATH}/lib/libhdf5.a -lz -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib"
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
#$CXX modify_mesh.cpp $INCLHDF $OPTLEV $LINKHDF -o modify_mesh.exe

##### partitioning
$CXX partition_hybrid_model_2part.cpp $INCLMETIS $INCLHDF $LINKHDF $LINKMETIS $OPTLEV $OPTCPP -o ./partition_hybrid_model_2part.exe

##### postprocessing
#$CXX postprocess_hybrid_model.cpp $OMP $OPTEV $OPTCPP -o ./postprocess_hybrid_model.exe

