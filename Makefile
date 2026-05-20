# Compiler
CXX      = nvcc 

# AOCL Installation Prefix
AOCL_DIR = /global/cfs/cdirs/m5214/hbhat512

# CUDA 12.9 directory
NVIDIA_DIR = /opt/nvidia/hpc_sdk/Linux_x86_64/25.5/math_libs/12.9/targets/x86_64-linux

# Compiler flags
# -arch=sm_90 is the critical flag for the H200 (Hopper architecture)
# -Xcompiler forwards CPU-specific flags to the underlying host compiler
CXXFLAGS = -O3 -m64 -std=c++20 -arch=sm_80 \
           -Xcompiler "-march=native -fopenmp" \
           -I $(AOCL_DIR)/include \
           -I $(AOCL_DIR)/include/eigen3 \
           -I $(NVIDIA_DIR)/include

# Linker flags
# Link order matters: High-level math (flame) -> Low-level math (blis) -> Utilities (aoclutils)
# -lcusolver handles the GPU math, and -Xcompiler "-fopenmp" ensures the CPU threads link properly
LDFLAGS  = -L$(AOCL_DIR)/lib \
           -L$(AOCL_DIR)/lib64 \
           -L$(NVIDIA_DIR)/lib \
           -lcnpy -lflame -lblis-mt -laoclutils -lcusolver \
           -Xcompiler "-fopenmp"

# Target
TARGET   = memoryFF
SRC      = memoryFFgpu.cpp

# Default rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Cleanup
clean:
	rm -f $(TARGET)
