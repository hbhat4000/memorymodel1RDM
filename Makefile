# Compiler (AOCC Clang C++)
CXX      = clang++

# AOCL Installation Prefix
AOCL_DIR = /global/cfs/cdirs/m5214/hbhat512

# Compiler flags
# -march=native replaces Intel's -xHost
# -fopenmp replaces Intel's -qopenmp
CXXFLAGS = -O3 -march=native -m64 -std=c++20 -fopenmp \
           -I $(AOCL_DIR)/include \
           -I $(AOCL_DIR)/include/eigen3

# Linker flags
# Link order matters: High-level math (flame) -> Low-level math (blis) -> Utilities (aoclutils)
LDFLAGS  = -L$(AOCL_DIR)/lib \
           -lcnpy -lflame -lblis -laoclutils -fopenmp

# Target
TARGET   = memoryFO
SRC      = memoryFO.cpp

# Default rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Cleanup
clean:
	rm -f $(TARGET)
