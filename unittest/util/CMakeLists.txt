set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/testbin)
set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
set(TEST_DIR ${PI_QMC_SOURCE_DIR}/unittest)
include_directories(${SRC_DIR} ${TEST_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(utilTests EXCLUDE_FROM_ALL
    PermutationTest.cc
    ${SRC_DIR}/util/Permutation.cc
    HungarianTest.cc
    ${SRC_DIR}/util/Hungarian.cc
    AperiodicGaussianTest.cc
    ${SRC_DIR}/util/AperiodicGaussian.cc
    PeriodicGaussianTest.cc
    ${SRC_DIR}/util/PeriodicGaussian.cc
    SuperCellTest.cc
    ${SRC_DIR}/util/SuperCell.cc
)

target_link_libraries(utilTests gtest gtest_main -lpthread)  
target_link_libraries(utilTests ${BLAS_LIB})
target_link_libraries(utilTests ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(utilTests ${LIBF2C_LIB})
endif()
add_dependencies(unittest utilTests)
