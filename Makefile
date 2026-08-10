# Compiler
CXX      = icpx

# base directory for Eigen, cnpy and cxxopts
BASE_DIR = /home/hbhat
# BASE_DIR = /u/hbhat

# Compiler flags
CXXFLAGS = -O3 -m64 -std=c++20 \
           -xHost -qopenmp -qmkl=parallel \
           -I $(BASE_DIR)/include \
           -I $(BASE_DIR)/include/eigen3 \

# Linker flags
LDFLAGS  = -L$(BASE_DIR)/lib \
           -L$(BASE_DIR)/lib64 \
           -lcnpy -qopenmp -qmkl=parallel

# Target
TARGET   = memoryFF
SRC      = memoryFFV2.cpp

# Default rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Cleanup
clean:
	rm -f $(TARGET)
