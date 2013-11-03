add_subdirectory(action)
#add_subdirectory(algorithm)
add_subdirectory(advancer)
add_subdirectory(base)
#add_subdirectory(demo)
add_subdirectory(emarate)
#add_subdirectory(estimator)
add_subdirectory(fixednode)
add_subdirectory(parser)
#add_subdirectory(spin)
add_subdirectory(stats)
add_subdirectory(util)

set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
set(TEST_DIR ${PI_QMC_SOURCE_DIR}/unit-tests)

include_directories(${SRC_DIR} ${TEST_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(unittest_pi EXCLUDE_FROM_ALL 
    unittest_main.cc
    ${sources}
)

target_link_libraries(unittest_pi 
    parser emarate action advancer algorithm demo estimator fixednode 
    spin stats base util)

target_link_libraries(unittest_pi gtest -lpthread)
target_link_libraries(unittest_pi ${BLAS_LIB})
target_link_libraries(unittest_pi ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(unittest_pi ${LIBF2C_LIB})
endif()
target_link_libraries(unittest_pi ${GSL_LIB})

target_link_libraries(unittest_pi ${FFTW3_LIB})
    
add_dependencies(unittest unittest_pi)
