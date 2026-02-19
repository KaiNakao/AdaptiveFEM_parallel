
####### Compiler, tools and options

CXX           = mpiicpc -cxx=icpx
CC            = mpiicc

HDFPATH=../hdf5-1.10.7_lib_intel/
INCLHDF=-I${HDFPATH}/include
LINKHDF=-L${HDFPATH}/lib/ -lhdf5_hl -lhdf5 -lrt -ldl -lm -Wl,-rpath -Wl,${HDFPATH}/lib
OPTLEV=-O3 -fopenmp -mkl
OPTCXX=-Wall -W -std=c++17 -mcmodel=large

#CXXFLAGS      = -Wall -W -mcmodel=large -O3 -std=c++17 -qopenmp -mkl -DNO_SMOOTHING
# CXXFLAGS      = -Wall -W -mcmodel=large -O3 -std=c++17 -qopenmp -mkl 
CXXFLAGS	  = $(OPTLEV) $(OPTCXX) $(INCLHDF) $(LINKHDF)
CFLAGS        = $(OPTLEV) -mcmodel=large $(INCLHDF)
# CXXFLAGS      = -Wall -W -mcmodel=large -g -std=c++17 -qopenmp -mkl

TAR           = tar -cf
COMPRESS      = gzip -9f


####### Directories

SOURCES_DIR = src/
BUILD_DIR = build/

####### Files

SOURCES_CPP     = $(wildcard $(SOURCES_DIR)/*.cpp)
SOURCES_C       = $(wildcard $(SOURCES_DIR)/*.c)

OBJECTS_CPP     = $(patsubst $(SOURCES_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES_CPP))
OBJECTS_C       = $(patsubst $(SOURCES_DIR)/%.c, $(BUILD_DIR)/%.o, $(SOURCES_C))
OBJECTS         = $(OBJECTS_CPP) $(OBJECTS_C)
TARGET = main

####### Build rules

all: $(BUILD_DIR) $(TARGET)

$(TARGET): $(OBJECTS) $(REFINER_OBJECTS) $(UTIL_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

####### Cleaning Rules

clean:
	rm -r $(BUILD_DIR) $(TARGET)
