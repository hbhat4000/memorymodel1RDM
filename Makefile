# Compiler
CXX      = icpx

# Compiler flags
CXXFLAGS = -O3 -xHost -m64 -std=c++20 -qopenmp \
           -I /home/hbhat/include/eigen3/ \
           -I /home/hbhat/include/

# Linker flags
LDFLAGS  = -L/home/hbhat/lib -lcnpy -qmkl -qopenmp

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

