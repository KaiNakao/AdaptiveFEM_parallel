
####### Compiler, tools and options

CXX           = icpx
#CXXFLAGS      = -Wall -W -mcmodel=large -O3 -std=c++17 -qopenmp -mkl -DNO_SMOOTHING
CXXFLAGS      = -Wall -W -mcmodel=large -O3 -std=c++17 -qopenmp -mkl
# CXXFLAGS      = -Wall -W -mcmodel=large -g -std=c++17 -qopenmp -mkl
TAR           = tar -cf
COMPRESS      = gzip -9f


####### Directories

SOURCES_DIR = src/
BUILD_DIR = build/

####### Files

SOURCES         = $(wildcard $(SOURCES_DIR)/*.cpp)

OBJECTS         = $(patsubst $(SOURCES_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))
TARGET = main

####### Build rules

all: $(BUILD_DIR) $(TARGET)

$(TARGET): $(OBJECTS) $(REFINER_OBJECTS) $(UTIL_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

####### Cleaning Rules

clean:
	rm -r $(BUILD_DIR) $(TARGET)
