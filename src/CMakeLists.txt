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

target_link_libraries(pi-qmc parser advancer estimator emarate 
    fixednode spin stats action algorithm base util demo)

target_link_libraries(pi-qmc ${LIBXML2_LIBRARIES} )
target_link_libraries(pi-qmc ${BLAS_LIB})
target_link_libraries(pi-qmc ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
  target_link_libraries(pi-qmc ${LIBF2C_LIB})
endif(EXISTS ${LIBF2C_LIB})
target_link_libraries(pi-qmc ${FFTW3_LIB})
target_link_libraries(pi-qmc ${HDF5_LIB})
target_link_libraries(pi-qmc ${GSL_LIB})

set_target_properties(pi-qmc
    PROPERTIES
    LINKER_LANGUAGE CXX)
