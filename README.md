# Current status

The current branch is tested with continuous integration on travis-ci.org:

    http://https://travis-ci.org/phys-tools/pi-qmc

    [![Build Status](https://travis-ci.org/phys-tools/pi-qmc.png)](https://travis-ci.org/phys-tools/pi-qmc)

# Obtaining the source code

*   Version control through gitHUB:

        git clone git@github.com:phys-tools/pi-qmc.git
        
*   Download ZIP File: https://github.com/phys-tools/pi-qmc/zipball/master
*   Download TAR Ball: https://github.com/phys-tools/pi-qmc/tarball/master


# Building with cmake (recommended)

Running "cmake ." in the source directory configures the program.
After the program is configured you can always edit settings in CMakeCache.txt.

You can set the intel compilers with:

    env CC=icc CXX=icpc cmake .

To generate an eclipse project with gcc-4.7 tools from mac ports:

    env CC=gcc-mp-4.7 CXX=g++-mp-4.7 cmake -G"Eclipse CDT4 - Unix Makefiles" .
    

# Building with legacy autotools (supported, but not  recommended)

Important configure flags:

*   Compile the code for 1, 2, 3 or 4 dimensional
    simulations. (Four dimensional simulations are mostly
    used for testing spin algorithms or other mathematical
    models.) 

        --with-ndim=NDIM
 
*   Enable MPI for parallel simulation.

        --enable-mpi

*   Use the SPRNG library for random numbers.

        --enable-sprng


*   If configure has trouble locating the proper blas
    and lapack libraries:

        --with-blas=<lib>
        --with-lapack=<lib>


