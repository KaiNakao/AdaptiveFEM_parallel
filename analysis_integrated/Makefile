NPC = 1
CRSSIZE = 400
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

####### Compiler, tools and options
CXX           = mpiicpc -cxx=icpx
CC            = icx
FC            = mpiifort -fc=ifx

HDFPATH=../hdf5-1.10.7_lib_intel/
INCLHDF=-I${HDFPATH}/include
LINKHDF=-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib

METISPATH=../metis-5.1.0_icc_32bit/
INCLMETIS=-I${METISPATH}/include
LINKMETIS=-L${METISPATH}/lib/ -lmetis

LINKFORT=-limf -lifcore -lifport -lifcoremt

OPTLEV=-O3 -fopenmp -qmkl
OPTF=-mcmodel=large -heap-arrays -zero
# OPTF =  -g -CB -traceback -mcmodel=large -i-dynamic -heap-arrays -zero -check bound
OPTCXX=-Wall -W -std=c++20 -mcmodel=large
OPTC=-Wall -W -mcmodel=large -restrict

OPTIONS=\
	-DREAL_4=real*4 -DREAL_4_SIZE=4 -DMPI_REAL_4=MPI_REAL4 \
	-DUSEIEPENTA -DNOOMPSAME -DTET4GRID -DMAXCOLOR=20 \
	-Dnpc=${NPC} -Dnkl=${NKL} -Dnsta=${NSTA} -DDIV_N=32 -DNOUSE_PA -DNONOUTPUTBLOCK=100 -DGFLIB \
	-DCRSSIZE=${CRSSIZE}

CXXFLAGS= $(OPTLEV) $(OPTCXX) $(INCLHDF) $(INCLMETIS) $(OPTIONS)
CFLAGS= $(OPTLEV) $(OPTC) $(INCLHDF) $(INCLMETIS) $(OPTIONS)
FFLAGS= $(OPTLEV) $(OPTF) $(INCLHDF) $(INCLMETIS) $(OPTIONS)
LDFLAGS= $(LINKHDF) $(LINKMETIS) $(LINKFORT)

####### Directories

SOURCES_DIR = src
BUILD_DIR = build
SRC_DIRS = $(SOURCES_DIR) $(SOURCES_DIR)/analysis_pointload $(SOURCES_DIR)/common $(SOURCES_DIR)/error_estimator $(SOURCES_DIR)/obs_displacement

####### Files

SOURCES_CPP     = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.cpp))
SOURCES_C       = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.c))
SOURCES_F       = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.F))
SOURCES_F90     = $(foreach dir,$(SRC_DIRS),$(wildcard $(dir)/*.F90))

OBJECTS_CPP     = $(patsubst $(SOURCES_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES_CPP))
OBJECTS_C       = $(patsubst $(SOURCES_DIR)/%.c, $(BUILD_DIR)/%.o, $(SOURCES_C))
OBJECTS_F       = $(patsubst $(SOURCES_DIR)/%.F, $(BUILD_DIR)/%.o, $(SOURCES_F))
OBJECTS_F90     = $(patsubst $(SOURCES_DIR)/%.F90, $(BUILD_DIR)/%.o, $(SOURCES_F90))
OBJECTS         = $(OBJECTS_CPP) $(OBJECTS_C) $(OBJECTS_F) $(OBJECTS_F90)
TARGET = main

####### Build rules

all: $(BUILD_DIR) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.c | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.F | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.F90 | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

####### Cleaning Rules

clean:
	rm -r $(BUILD_DIR)
	rm *__genmod.f90 *__genmod.mod
