add_subdirectory(demo)
add_subdirectory(algorithm)
add_subdirectory(action)
add_subdirectory(parser)
add_subdirectory(sampler)
add_subdirectory(emarate)
add_subdirectory(estimator)
add_subdirectory(fixednode)
add_subdirectory(spin)
add_subdirectory(stats)
add_subdirectory(util)

set (sources
    BeadFactory.cc
    DoubleParallelPaths.cc
    EnumeratedModelState.cc
    ParallelPaths.cc
    Paths.cc
    SHOPhase.cc
    SerialPaths.cc
    SimInfoWriter.cc
    SimulationInfo.cc
    Species.cc
    Spin4DPhase.cc
    SpinModelState.cc
    XMLWriter.cc
)

add_executable(pi ${sources} main.cc)

TARGET_LINK_LIBRARIES(pi ${LIBXML2_LIBRARIES} )
target_link_libraries(pi ${BLAS_LIB})
target_link_libraries(pi ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
  target_link_libraries(pi ${LIBF2C_LIB})
endif(EXISTS ${LIBF2C_LIB})
target_link_libraries(pi ${FFTW3_LIB})
target_link_libraries(pi ${HDF5_LIB})
target_link_libraries(pi ${GSL_LIB})

set_target_properties(pi
    PROPERTIES
    LINKER_LANGUAGE CXX)

target_link_libraries(pi demo parser sampler estimator emarate fixednode 
    spin stats action algorithm util)

#add_subdirectory(test)
INCLUDE_DIRECTORIES(${PI_BINARY_DIR} ${PI_SOURCE_DIR}/src)
