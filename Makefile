# Compiler
CXX      = icpx

# Compiler flags
CXXFLAGS = -O3 -xHost -m64 \
           -I/home/hbhat/data/cnpy \
           -I/home/hbhat/include \
           -I/home/hbhat/include/eigen3

# Linker flags
LDFLAGS  = -L/home/hbhat/data/cnpy -lcnpy -qmkl

# Target
TARGET   = memoryFF
SRC      = memoryFF.cpp

# Default rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

# Cleanup
clean:
	rm -f $(TARGET)

