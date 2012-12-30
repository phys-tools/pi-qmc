cmake_minimum_required (VERSION 2.6)

project (PI_QMC)
set (PI_QMC_VERSION_MAJOR 0)
set (PI_QMC_VERSION_MINOR 1)
enable_language(CXX)

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

set(BLITZ_INCLUDE_DIR ${PI_QMC_SOURCE_DIR}/contrib/blitz-0.9)

include(CheckIncludeFiles)
CHECK_INCLUDE_FILES(getopt.h HAVE_GETOPT_H)

configure_file (
  "${PI_QMC_SOURCE_DIR}/config_cmake.h.in"
  "${PI_QMC_BINARY_DIR}/config.h"
)

INCLUDE_DIRECTORIES(
    ${PI_QMC_BINARY_DIR}
    ${PI_QMC_SOURCE_DIR}/lib
    ${PI_QMC_SOURCE_DIR}/src
    ${BLITZ_INCLUDE_DIR})

add_definitions(-DHAVE_CONFIG_H)


set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PI_QMC_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PI_QMC_BINARY_DIR}/lib)


add_custom_target(unittest)
set (GTEST_SOURCE_DIR ${PI_QMC_SOURCE_DIR}/contrib/gtest-1.6.0)
add_subdirectory(${GTEST_SOURCE_DIR})

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/etc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/etc/Doxyfile @ONLY)
add_custom_target(doxygen
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/etc/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc/etc
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

add_subdirectory(src)
add_subdirectory(unit-tests)
add_subdirectory(test)

