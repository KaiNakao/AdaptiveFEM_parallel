
####### Compiler, tools and options

CXX           = g++
CXXFLAGS      = -Wall -W
INCPATH       = -I./src/inc
TAR           = tar -cf
COMPRESS      = gzip -9f


####### Output directly

SOURCES_DIR = src/
ESTIMATOR_DIR = src/estimator/
MARKER_DIR = src/marker/
REFINER_DIR = src/refiner/
UTIL_DIR = src/util/
BUILD_DIR = build/
BUILD_ESTIMATOR_DIR = build/estimator/
BUILD_MARKER_DIR = build/marker/
BUILD_REFINER_DIR = build/refiner/
BUILD_UTIL_DIR = build/util/

####### Files

SOURCES           = $(wildcard $(SOURCES_DIR)/*.cpp)
ESTIMATOR_SOURCES = $(wildcard $(ESTIMATOR_DIR)/*.cpp)
MARKER_SOURCES    = $(wildcard $(MARKER_DIR)/*.cpp)
REFINER_SOURCES   = $(wildcard $(REFINER_DIR)/*.cpp)
UTIL_SOURCES      = $(wildcard $(UTIL_DIR)/*.cpp)

OBJECTS           = $(patsubst $(SOURCES_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SOURCES))
ESTIMATOR_OBJECTS = $(patsubst $(ESTIMATOR_DIR)/%.cpp, $(BUILD_ESTIMATOR_DIR)/%.o, $(ESTIMATOR_SOURCES))
MARKER_OBJECTS    = $(patsubst $(MARKER_DIR)/%.cpp, $(BUILD_MARKER_DIR)/%.o, $(MARKER_SOURCES))
REFINER_OBJECTS   = $(patsubst $(REFINER_DIR)/%.cpp, $(BUILD_REFINER_DIR)/%.o, $(REFINER_SOURCES))
UTIL_OBJECTS      = $(patsubst $(UTIL_DIR)/%.cpp, $(BUILD_UTIL_DIR)/%.o, $(UTIL_SOURCES))
TARGET = main

####### Build rules

all: $(BUILD_DIR) $(TARGET)

$(TARGET): $(OBJECTS) $(ESTIMATOR_OBJECTS) $(MARKER_OBJECTS) $(REFINER_OBJECTS) $(UTIL_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SOURCES_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(BUILD_ESTIMATOR_DIR)/%.o: $(ESTIMATOR_DIR)/%.cpp | $(BUILD_ESTIMATOR_DIR)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(BUILD_MARKER_DIR)/%.o: $(MARKER_DIR)/%.cpp | $(BUILD_MARKER_DIR)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(BUILD_REFINER_DIR)/%.o: $(REFINER_DIR)/%.cpp | $(BUILD_REFINER_DIR)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(BUILD_UTIL_DIR)/%.o: $(UTIL_DIR)/%.cpp | $(BUILD_UTIL_DIR)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR) $(BUILD_ESTIMATOR_DIR) $(BUILD_MARKER_DIR) $(BUILD_REFINER_DIR) $(BUILD_UTIL_DIR)

####### Cleaning Rules

clean:
	rm -r $(BUILD_DIR) $(TARGET)