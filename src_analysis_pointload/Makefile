# Compile on x86-64
# use 64 bit integers by default

NPC = 1

ifeq ($(NPC), 1)
  NSTA = 0
  NKL = 0
else ifeq ($(NPC), 21)
  NSTA = 5
  NKL = 5
else ifeq ($(NPC), 288)
  NSTA = 5
  NKL = 10
endif

hostname = $(shell hostname)
ifeq ($(findstring ibis, $(hostname)), ibis)
  ##### path to libraries
  CXX = icpc
  CC = icc
  FC = mpiifort
  HDFPATH = ../hdf5-1.10.7_lib_intel
  METISPATH = ../metis-5.1.0_icc_32bit

  ##### DEBUG
  OPTLEV = -g -fopenmp -DOUTPUT_LOG -fpe0
  OPTF =  -mcmodel=large -i-dynamic -heap-arrays -zero -check -warn all
  OPTF =  -check -warn all
  OPTF =  -g -CB -traceback
  OPTF =  -g -CB -traceback -mcmodel=large -i-dynamic -heap-arrays -zero -check -warn all
  OPTF =  -g -CB -traceback -mcmodel=large -i-dynamic -heap-arrays -zero -check bound

  ##### RELEASE
#   OPTLEV = -O3 -fopenmp
  OPTLEV = -O3 -fopenmp -fpe0
  OPTF =  -mcmodel=large -heap-arrays -zero

  OPTC =
  OPTCPP =
  # for linking fortran subroutines to c++ program
  LINKFORT = -limf -lifcore -lifport -lifcoremt
endif

ifeq ($(findstring fn01sv, $(hostname)), fn01sv)
#ifeq ($(findstring bmt4-ln, $(hostname)), bmt4-ln)
  CC = fccpx
  FC = mpifrtpx
  HDFPATH = ${HOME}/local
  METISPATH = ${HOME}/local
  OPTLEV =
  OPTF =

#   # for release
  OPTLEV = -Kfast,openmp,parallel,ocl
  OPTF = -Kautoobjstack -Nlst=t

  # for debug
#   OPTLEV = -H aesofux -g -NRtrap
#   OPTF =
endif

INCLHDF = -I${HDFPATH}/include
LINKHDF = -L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib
INCLMETIS = -I${METISPATH}/include
LINKMETIS = -L${METISPATH}/lib/ -lmetis


OPTIONS = \
	-DREAL_4=real*4 -DREAL_4_SIZE=4 -DMPI_REAL_4=MPI_REAL4 \
	-DUSEIEPENTA -DNOOMPSAME -DTET4GRID -DMAXCOLOR=20 \
	-Dnpc=${NPC} -Dnkl=${NKL} -Dnsta=${NSTA} -DDIV_N=32 -DNOUSE_PA -DNONOUTPUTBLOCK=100 -DREAD_NODAL_FORCE \
#	-DREAL_4=real*8 -DREAL_4_SIZE=8 -DMPI_REAL_4=MPI_REAL8 \

.SUFFIXES: .o .F .F90 .cpp

SRC_SUB_F = \
	setcolor_1th.F \
	sync_parallel.F \
	compAmats.F \
	compAmatd.F \
	read_data.F \
	other_functions.F \
	IEpenta.F \
	IEpenta2nd.F \
	IEpenta2nd_KG.F \
	IEpenta2nd_s.F \
	IEpenta1st.F \
	IEpenta1st_KG.F \
	IEpenta1st_s.F \
	nontetra10_double.F \
	nontetra4_double.F \
	pickupfaultelement.F \
	pointsource.F \
	tet10ku_single.F \
	tet10ku_double.F \
	compAmat4s.F \
	tet4kud_single.F \
	compsig_npc288.F \
	compsig_npc21.F \
	set_CRS.F \
	IPCG.F \
	CGs.F \
	compdiag.F \
	output_obs.F \

SRC_SUB_F90 = \
	pointload.F90 \
	read_nodal_force.F90


SRC_C = hdf_lib.c scale_value.c

OBJS = crustFEM.o ${SRC_C:.c=.o} ${SRC_SUB_F:.F=.o} ${SRC_SUB_F90:.F90=.o}
OBJS_DEBUG = debug.o ${SRC_C:.c=.o} ${SRC_SUB_F:.F=.o} ${SRC_SUB_F90:.F90=.o}
OBJS_DEBUG2 = debug_split_crs.o ${SRC_C:.c=.o} ${SRC_SUB_F:.F=.o} ${SRC_SUB_F90:.F90=.o}

#######################################################################
default: FORCE
	make -j GAMERA

GAMERA: ${OBJS}
	${FC} ${OBJS} $(OPTLEV) ${LINKHDF} ${LINKMETIS} -o ./analysis.exe.npc${NPC}

debug: ${OBJS_DEBUG}
	${FC} ${OBJS_DEBUG} $(OPTLEV) ${LINKHDF} ${LINKMETIS} -o ./debug.exe

debug2: ${OBJS_DEBUG2}
	${FC} ${OBJS_DEBUG2} $(OPTLEV) ${LINKHDF} ${LINKMETIS} -o ./debug_split_crs.exe

.F.o:
	${FC} ${OPTF} ${OPTLEV} ${OPTIONS} -c $<

.F90.o:
	${FC} ${OPTF} ${OPTLEV} ${OPTIONS} -c $<

hdf_lib.o: hdf_lib.c
	${CC} hdf_lib.c -c ${INCLHDF} ${OPTC} ${OPTLEV} $<

scale_value.o: scale_value.c
	${CC} -c -O3 -fopenmp -restrict ${OPTIONS} $<


#######################################################################
clean: FORCE
	rm -f *.o *.exe *__genmod*

all: FORCE
	make clean
	make default

check: FORCE
	@echo FC: ${FC}

FORCE:
