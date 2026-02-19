####### Compiler, tools and options

CXX           = mpiicpc -cxx=icpx
CC            = mpiicc

HDFPATH=../hdf5-1.10.7_lib_intel/
INCLHDF=-I${HDFPATH}/include
LINKHDF=-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib

METISPATH=../metis-5.1.0_icc_32bit/
INCLMETIS=-I${METISPATH}/include
LINKMETIS=-L${METISPATH}/lib/ -lmetis

PARMETISPATH=../parmetis_icc_32bit/
INCLPARMETIS=-I${PARMETISPATH}/include
LINKPARMETIS=-L${PARMETISPATH}/lib -lparmetis

GKLIBPATH=../gklib_icc/
INCLGKLIB=-I${GKLIBPATH}/include
LINKGKLIB=-L${GKLIBPATH}/lib -lGKlib

OPTLEV=-O3 -fopenmp -mkl
OPTCXX=-Wall -W -std=c++17 -mcmodel=large

CXXFLAGS	  = $(OPTLEV) $(OPTCXX) $(INCLHDF) $(INCLMETIS) $(INCLPARMETIS) $(INCLGKLIB)
CFLAGS        = $(OPTLEV) -mcmodel=large $(INCLHDF)
LDFLAGS       = $(LINKHDF) $(LINKPARMETIS) $(LINKMETIS) $(LINKGKLIB)

####### Directories

SOURCES_DIR = src/
BUILD_DIR = build/

####### Files

SOURCES_CPP     = $(wildcard $(SOURCES_DIR)/*.cpp)
SOURCES_C       = $(wildcard $(SOURCES_DIR)/*.c)

OBJECTS_CPP     = $(patsubst $(SOURCES_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES_CPP))
OBJECTS_C       = $(patsubst $(SOURCES_DIR)/%.c, $(BUILD_DIR)/%.o, $(SOURCES_C))
OBJECTS         = $(OBJECTS_CPP) $(OBJECTS_C)
TARGET = repartition

####### Build rules

all: $(BUILD_DIR) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

####### Cleaning Rules

clean:
	rm -r $(BUILD_DIR) $(TARGET)
