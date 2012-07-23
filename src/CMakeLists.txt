add_subdirectory(action)
add_subdirectory(advancer)
add_subdirectory(algorithm)
add_subdirectory(base)
add_subdirectory(demo)
add_subdirectory(emarate)
add_subdirectory(estimator)
add_subdirectory(fixednode)
add_subdirectory(parser)
add_subdirectory(spin)
add_subdirectory(stats)
add_subdirectory(util)

add_executable(pi main.cc)

target_link_libraries(pi parser demo advancer estimator emarate 
    fixednode spin stats action algorithm base util)

target_link_libraries(pi ${LIBXML2_LIBRARIES} )
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
