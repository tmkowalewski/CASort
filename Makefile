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

# Target executable name
TARGET   := $(BIN_DIR)/$(PROJECT_NAME)

# Find all source files in the src directory
SOURCES  := $(wildcard $(SRC_DIR)/*.cpp)
# All object files
OBJECTS  := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# Default target
all: $(TARGET)

# Debug target
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: clean $(TARGET)

# Link object files into the executable
$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Compile source files into object files
# This rule is used for all .o files, which will be PIC due to -fPIC in CXXFLAGS
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
	$(TARGET) examples/root_data_70Ge_run001.mvmelst.bin_tree.root out.root

.PHONY: all clean library debug
