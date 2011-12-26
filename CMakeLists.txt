cmake_minimum_required (VERSION 2.6)

project (PI_QMC)
set (PI_QMC_VERSION_MAJOR 0)
set (PI_QMC_VERSION_MINOR 1)
enable_language(CXX)
enable_language(Fortran)

find_library(LAPACK_LIB
NAMES
    lapack
PATHS
    ENV LAPACK
    /usr/lib
    /usr/local/lib
    c:/libraries/lib
    ENV LAPACK
)

if(NOT EXISTS ${LAPACK_LIB})
    message("** Error finding PI_QMC required library:")
    message("**   Try setting the LAPACK environtment variable")
    message("**   to the directory containing lapack")
    message(FATAL_ERROR "** PI_QMC requires the LAPACK libary")
else()
    message("Found LAPACK: " ${LAPACK_LIB})
endif(NOT EXISTS ${LAPACK_LIB})

find_library(BLAS_LIB
NAMES
    blas
PATHS
    /usr/lib
    /usr/local/lib
    c:/libraries/lib
    ENV BLAS
PATHS_SUFFIXES
    lib
)

if(NOT EXISTS ${BLAS_LIB})
    message("** Error finding PI_QMC required library:")
    message("**   Try setting the BLAS environtment variable")
    message("**   to the directory containing blas")
    message(FATAL_ERROR "** PI_QMC requires the BLAS libary")
else()
    message("Found BLAS: " ${BLAS_LIB})    
endif(NOT EXISTS ${BLAS_LIB})

find_library(LIBF2C_LIB
NAMES
    libf2c
PATHS
    /usr/lib
    /usr/local/lib
    c:/libraries/lib
    ENV BLAS
PATHS_SUFFIXES
    lib
)

FIND_PACKAGE(LibXml2 REQUIRED)
INCLUDE_DIRECTORIES (SYSTEM "${LIBXML2_INCLUDE_DIR}")

find_library(HDF5_LIB
NAMES
    hdf5
PATHS
    /usr/lib
    /usr/local/lib
    c:/libraries/lib
    ENV HDF5
PATHS_SUFFIXES
    lib
)

if(NOT EXISTS ${HDF5_LIB})
    message("** Error finding PI_QMC required library:")
    message("**   Try setting the HDF5 environtment variable")
    message("**   to the directory containing hdf5")
    message(FATAL_ERROR "** PI_QMC requires the hdf5 libary")
else()
    message("Found HDF5: " ${HDF5_LIB})    
endif(NOT EXISTS ${HDF5_LIB})


find_library(GSL_LIB
NAMES
    gsl
PATHS
    /usr/lib
    /usr/local/lib
    c:/libraries/lib
    ENV GSL
PATHS_SUFFIXES
    lib
)

if(NOT EXISTS ${GSL_LIB})
    message("** Error finding PI_QMC required library:")
    message("**   Try setting the GSL environtment variable")
    message("**   to the directory containing gsl")
    message(FATAL_ERROR "** PI_QMC requires the GSL libary")
else()
    message("Found GSL: " ${GSL_LIB})    
endif(NOT EXISTS ${GSL_LIB})

find_library(FFTW3_LIB
NAMES
    fftw3
PATHS
    /usr/lib
    /usr/local/lib
    c:/libraries/lib
    ENV FFTW3
PATHS_SUFFIXES
    lib
)

if(NOT EXISTS ${FFTW3_LIB})
    message("** Error finding PI_QMC required library:")
    message("**   Try setting the FFTW3 environtment variable")
    message("**   to the directory containing fftw3")
    message(FATAL_ERROR "** PI_QMC requires the FFTW3 libary")
else()
    message("Found FFTW3: " ${FFTW3_LIB})    
endif(NOT EXISTS ${FFTW3_LIB})


#find_package(FFTW3 REQUIRED)

configure_file (
  "${PI_QMC_SOURCE_DIR}/config_cmake.h.in"
  "${PI_QMC_BINARY_DIR}/config.h"
)
INCLUDE_DIRECTORIES(${PI_QMC_BINARY_DIR} ${PI_QMC_SOURCE_DIR}/lib ${PI_QMC_SOURCE_DIR}/src)

add_definitions(-DHAVE_CONFIG_H)


set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PI_QMC_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PI_QMC_BINARY_DIR}/lib)


add_subdirectory(lib)
add_subdirectory(src)
