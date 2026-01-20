# Makefile for CASort project
PROJECT_NAME := CASort

# Directories
SRC_DIR  := src
OBJ_DIR  := obj
BIN_DIR  := bin
#LIB_DIR  := lib
INC_DIR := include

# Compiler and flags
CXX      := g++
CXXFLAGS := `root-config --cflags` -I./include
LDFLAGS  := `root-config --libs`
DEBUGFLAGS := -g -O0

# Target executable names
TARGET   := $(BIN_DIR)/$(PROJECT_NAME)
TARGET_GAINMATCH := $(BIN_DIR)/GainMatch

# Find all source files in the src directory
SOURCES  := $(wildcard $(SRC_DIR)/*.cpp)

# Separate common sources (shared between binaries) from main files
COMMON_SOURCES := $(filter-out $(SRC_DIR)/CASort.cpp $(SRC_DIR)/GainMatch.cpp,$(SOURCES))
COMMON_OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(COMMON_SOURCES))

# Object collections per binary
CASORT_OBJECTS := $(COMMON_OBJECTS) $(OBJ_DIR)/CASort.o
GAINMATCH_OBJECTS := $(COMMON_OBJECTS) $(OBJ_DIR)/GainMatch.o

# Default target
all: $(TARGET) $(TARGET_GAINMATCH)

# Build only CASort binary
casort: $(TARGET)

# Build only GainMatch binary
gainmatch: $(TARGET_GAINMATCH)

# Debug target
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: clean all

# Link object files into the executables
$(TARGET): $(CASORT_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(TARGET_GAINMATCH): $(GAINMATCH_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS) -lSpectrum

# Compile source files into object files
# This rule is used for all .o files, which will be PIC due to -fPIC in CXXFLAGS
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

install: $(TARGET) $(TARGET_GAINMATCH)
	@mkdir -p ~/.local/bin
	@cp $(TARGET) ~/.local/bin/$(PROJECT_NAME)
	@cp $(TARGET_GAINMATCH) ~/.local/bin/GainMatch
	@mkdir -p ~/.local/include/${PROJECT_NAME}
	@cp $(INC_DIR)/*.hpp ~/.local/include/${PROJECT_NAME}

uninstall:
	@rm -f ~/.local/bin/$(PROJECT_NAME)
	@rm -f ~/.local/bin/GainMatch
	@rm -rf ~/.local/lib/$(PROJECT_NAME)
	@rm -rf ~/.local/include/${PROJECT_NAME}

# Clean up build artifacts
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

test: $(TARGET)
	$(TARGET) ~/TUNL/Data/NRF/70Ge/energy_calibration/calibrations . examples 1 out.001.root

.PHONY: all clean debug casort gainmatch
