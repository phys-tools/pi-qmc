set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/test)
set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
include_directories(${SRC_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(emaRateTests EXCLUDE_FROM_ALL
    EMARateActionTest.cc
    EMARateMoverTest.cc
    EMARateEstimatorTest.cc
    EMARateTestBeadPositioner.cc
    ${SRC_DIR}/sampler/test/MultiLevelSamplerFake.cc
    ${SRC_DIR}/emarate/EMARateAction.cc
    ${SRC_DIR}/emarate/EMARateEstimator.cc
    ${SRC_DIR}/emarate/EMARateMover.cc
    ${SRC_DIR}/Beads.h
    ${SRC_DIR}/BeadFactory.cc
    ${SRC_DIR}/SimulationInfo.cc
    ${SRC_DIR}/Species.cc
    ${SRC_DIR}/Paths.cc
    ${SRC_DIR}/SerialPaths.cc
    ${SRC_DIR}/util/SuperCell.cc
    ${SRC_DIR}/util/Permutation.cc
    ${SRC_DIR}/util/RandomNumGenerator.cc
    ${SRC_DIR}/stats/ScalarEstimator.cc
    ${SRC_DIR}/sampler/PermutationChooser.cc
    ${SRC_DIR}/sampler/ParticleChooser.cc
)
target_link_libraries(emaRateTests gtest gtest_main -lpthread)
target_link_libraries(emaRateTests ${BLAS_LIB})
target_link_libraries(emaRateTests ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(emaRateTests ${LIBF2C_LIB})
endif()
add_dependencies(unittest emaRateTests)
