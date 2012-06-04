set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/test)
set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
include_directories(${SRC_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(fixedNodeTests EXCLUDE_FROM_ALL
    AugmentedNodesTest.cc
    Atomic1sDMTest.cc
    ${SRC_DIR}/fixednode/AtomicOrbitalDM.cc
    ${SRC_DIR}/fixednode/Atomic1sDM.cc
    Atomic2spDMTest.cc
    ${SRC_DIR}/fixednode/Atomic2spDM.cc
)
target_link_libraries(fixedNodeTests gtest gtest_main -lpthread)
target_link_libraries(fixedNodeTests ${BLAS_LIB})
target_link_libraries(fixedNodeTests ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(fixedNodeTests ${LIBF2C_LIB})
endif()

add_dependencies(unittest fixedNodeTests)

