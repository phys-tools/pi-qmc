set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/test)
set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
include_directories(${SRC_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(utilTests EXCLUDE_FROM_ALL
    PermutationTest.cc
    ../Permutation.cc
    HungarianTest.cc
    ../Hungarian.cc
    AperiodicGaussianTest.cc
    ../AperiodicGaussian.cc
    PeriodicGaussianTest.cc
    ../PeriodicGaussian.cc
    SuperCellTest.cc
    ../SuperCell.cc
)

target_link_libraries(utilTests gtest gtest_main -lpthread)  
target_link_libraries(utilTests ${BLAS_LIB})
target_link_libraries(utilTests ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(utilTests ${LIBF2C_LIB})
endif()
add_dependencies(unittest utilTests)
