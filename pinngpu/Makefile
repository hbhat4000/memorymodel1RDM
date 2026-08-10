# DELTA (INTEL)
# HOST_COMPILER = /sw/rh9.4/spack/v1.0.0/sw/linux-x86_64_v2/intel-oneapi-compilers-2025.2.0-a2xli3d/compiler/2025.2/bin/icpx

# PINNACLES (INTEL)
HOST_COMPILER = /opt/spack-v1.1.1/opt/spack/linux-cascadelake/intel-oneapi-compilers-2025.3.1-hc2qvxmvb44sbrcv4ialqzent27vxc2v/compiler/2025.3/bin/icpx

# Compiler
CXX      = nvcc 

# DELTA CUDA directory
# NVIDIA_DIR = /opt/nvidia/hpc_sdk/Linux_x86_64/25.3/math_libs
# PINNACLES CUDA directory
NVIDIA_DIR = /opt/spack-v1.1.1/opt/spack/linux-cascadelake/nvhpc-25.7-kfehydowyjvjmze7hn7po2jkarxfml77/Linux_x86_64/25.7/math_libs/12.9/

# DELTA home directory
# HOME_DIR = /u/hbhat

# PINNACLES home directory
HOME_DIR = /home/hbhat

# Compiler flags
# -Xcompiler forwards CPU-specific flags to the underlying host compiler
CXXFLAGS = -x cu -O3 -m64 -std=c++20 -arch=sm_80 \
           -allow-unsupported-compiler --expt-relaxed-constexpr \
           -ccbin $(HOST_COMPILER) \
           -Xcompiler "-xHost -qopenmp -qmkl=parallel" \
           -I $(HOME_DIR)/include -I $(HOME_DIR)/include/eigen3 \

# Linker flags
LDFLAGS  = -L$(HOME_DIR)/lib -L$(HOME_DIR)/lib64 \
           -L$(NVIDIA_DIR)/lib64 \
           -lcnpy -lcusolver -lcublas \
           -Xcompiler "-qopenmp -qmkl=parallel"

# Target
TARGET   = memoryFF
SRC      = memoryFFgpuV2.cpp

# Default rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Cleanup
clean:
	rm -f $(TARGET)
