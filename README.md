At the moment, the Makefile specifically uses the Intel compiler and MKL available.

To get this to run on the Pinnacles cluster at UC Merced, after logging in, load the following modules at the command line:

module load intel-oneapi-compilers/2025.3.1
module load intel-oneapi-mkl/2025.3.1
module load intel-oneapi-tbb/2022.3.0
module load intel-oneapi-vtune/2025.3.0

Three other prerequisites:
- Installing the latest version of Eigen: first check out the repo

git clone https://gitlab.com/libeigen/eigen.git

Then follow the instructions using cmake and make install into your home directory.

- Build and install cnpy: https://github.com/rogersce/cnpy

- Build and install cxxopts: https://github.com/jarro2783/cxxopts

Check the paths in the Makefile.  Then use "make" to compile the code.

To run the code, the command-line syntax with required arguments is:

./memoryFF --dt [float] --delay [integer] --infile [string]

There are two optional arguments, --verbose to print out lots of information and --tol [float] to set the SVD tolerance.
