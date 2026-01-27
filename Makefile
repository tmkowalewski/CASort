# Makefile for CASort project
PROJECT_NAME := CASort

# Directories
SRC_DIR  := src
OBJ_DIR  := obj
LIB_DIR  := lib
BIN_DIR  := bin
INC_DIR  := include

# Compiler and flags
CXX       := g++
CXXFLAGS  := `root-config --cflags` -I./include -fPIC
LDFLAGS   := `root-config --libs` -shared
DEBUGFLAGS := -g -O0

# Target shared library name
TARGET := $(LIB_DIR)/lib$(PROJECT_NAME).so

# Source and object files
SOURCES  := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS  := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

# Default target
all: $(TARGET)

# Debug target
debug: CXXFLAGS += $(DEBUGFLAGS)
debug: clean all

# Link object files into the shared library
$(TARGET): $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

install: $(TARGET)
	@mkdir -p ~/.local/lib
	@cp $(TARGET) ~/.local/lib
	@mkdir -p ~/.local/include/${PROJECT_NAME}
	@cp $(INC_DIR)/*.hpp ~/.local/include/${PROJECT_NAME}

uninstall:
	@rm -f ~/.local/lib/lib$(PROJECT_NAME).so
	@rm -rf ~/.local/include/${PROJECT_NAME}

# Clean up build artifacts
clean:
	rm -rf $(OBJ_DIR) $(LIB_DIR) $(BIN_DIR)

.PHONY: all clean debug install uninstall

