# Makefile for CASort project
PROJECT_NAME := CASort

# Directories
SRC_DIR  := src
OBJ_DIR  := obj
BIN_DIR  := bin
INC_DIR  := include

# Compiler and flags
CXX       := g++
CXXFLAGS  := `root-config --cflags` -I./include
LDFLAGS   := `root-config --libs`
DEBUGFLAGS := -g -O0

# Target executable name
TARGET := $(BIN_DIR)/$(PROJECT_NAME)

# Source and object files
SOURCES  := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS  := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# Default target
all: $(TARGET)

# Debug target
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: clean all

# Link object files into the executable
$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

install: $(TARGET)
	@mkdir -p ~/.local/bin
	@cp $(TARGET) ~/.local/bin/$(PROJECT_NAME)
	@mkdir -p ~/.local/include/${PROJECT_NAME}
	@cp $(INC_DIR)/*.hpp ~/.local/include/${PROJECT_NAME}

uninstall:
	@rm -f ~/.local/bin/$(PROJECT_NAME)
	@rm -rf ~/.local/lib/$(PROJECT_NAME)
	@rm -rf ~/.local/include/${PROJECT_NAME}

# Clean up build artifacts
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

test: $(TARGET)
	$(TARGET) ~/TUNL/Data/NRF/70Ge/energy_calibration/calibrations . examples 1 out.001.root

.PHONY: all clean debug
