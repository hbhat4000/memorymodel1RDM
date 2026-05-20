# underlying C++ compiler
HOST_COMPILER = /opt/spack-v1.1.1/opt/spack/linux-cascadelake/intel-oneapi-compilers-2025.3.1-hc2qvxmvb44sbrcv4ialqzent27vxc2v/compiler/2025.3/bin/icpx

# Compiler
CXX      = nvcc 

# base directory for Eigen, cnpy and cxxopts
BASE_DIR = /home/hbhat

# CUDA 12.9 directory
NVIDIA_DIR = /opt/spack-v1.1.1/opt/spack/linux-cascadelake/nvhpc-25.7-kfehydowyjvjmze7hn7po2jkarxfml77/Linux_x86_64/2025/math_libs/12.9

# Compiler flags
# -ccbin tells nvcc to use AOCC for the host code
CXXFLAGS = -x cu -O3 -m64 -std=c++20 -arch=sm_80 -allow-unsupported-compiler -ccbin $(HOST_COMPILER) \
           -Xcompiler "-xHost -qopenmp -qmkl=parallel" \
           -I $(BASE_DIR)/include \
           -I $(BASE_DIR)/include/eigen3 \
           -I $(NVIDIA_DIR)/include

# Linker flags
LDFLAGS  = -L$(BASE_DIR)/lib \
           -L$(BASE_DIR)/lib64 \
           -L$(NVIDIA_DIR)/lib64 \
           -lcnpy -lcusolver -lcublas \
           -Xcompiler "-qopenmp -qmkl=parallel"

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
