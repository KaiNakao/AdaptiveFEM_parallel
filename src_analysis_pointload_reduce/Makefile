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

####### Compiler, tools and options
CXX           = mpiicpc
CC            = mpiicc
FC            = mpiifort

HDFPATH=../hdf5-1.10.7_lib_intel/
INCLHDF=-I${HDFPATH}/include
LINKHDF=-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib
METISPATH=../metis-5.1.0_icc_32bit/
INCLMETIS=-I${METISPATH}/include
LINKMETIS=-L${METISPATH}/lib/ -lmetis
OPTLEV=-O3 -fopenmp -qmkl -fpe0
OPTF=-mcmodel=large -heap-arrays -zero
OPTCXX=-Wall -W -std=c++17 -mcmodel=large
OPTC=-Wall -W -mcmodel=large -restrict
LINKFORT=-limf -lifcore -lifport -lifcoremt

OPTIONS=\
	-DREAL_4=real*4 -DREAL_4_SIZE=4 -DMPI_REAL_4=MPI_REAL4 \
	-DUSEIEPENTA -DNOOMPSAME -DTET4GRID -DMAXCOLOR=20 \
	-Dnpc=${NPC} -Dnkl=${NKL} -Dnsta=${NSTA} -DDIV_N=32 -DNOUSE_PA -DNONOUTPUTBLOCK=100 -DGFLIB

CXXFLAGS= $(OPTLEV) $(OPTCXX) $(INCLHDF) $(INCLMETIS) $(OPTIONS)
CFLAGS= $(OPTLEV) $(OPTC) $(INCLHDF) $(INCLMETIS) $(OPTIONS)
FFLAGS= $(OPTLEV) $(OPTF) $(INCLHDF) $(INCLMETIS) $(OPTIONS)
LDFLAGS= $(LINKHDF) $(LINKMETIS) $(LINKFORT)

####### Directories

SOURCES_DIR = src/
BUILD_DIR = build/

####### Files

SOURCES_CPP     = $(wildcard $(SOURCES_DIR)/*.cpp)
SOURCES_C       = $(wildcard $(SOURCES_DIR)/*.c)
SOURCES_F       = $(wildcard $(SOURCES_DIR)/*.F)
SOURCES_F90     = $(wildcard $(SOURCES_DIR)/*.F90)

OBJECTS_CPP     = $(patsubst $(SOURCES_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES_CPP))
OBJECTS_C       = $(patsubst $(SOURCES_DIR)/%.c, $(BUILD_DIR)/%.o, $(SOURCES_C))
OBJECTS_F       = $(patsubst $(SOURCES_DIR)/%.F, $(BUILD_DIR)/%.o, $(SOURCES_F))
OBJECTS_F90     = $(patsubst $(SOURCES_DIR)/%.F90, $(BUILD_DIR)/%.o, $(SOURCES_F90))
OBJECTS         = $(OBJECTS_CPP) $(OBJECTS_C) $(OBJECTS_F) $(OBJECTS_F90)
TARGET = main

####### Build rules

all: $(BUILD_DIR) $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.F | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.F90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

####### Cleaning Rules

clean:
	rm -r $(BUILD_DIR) $(TARGET)
