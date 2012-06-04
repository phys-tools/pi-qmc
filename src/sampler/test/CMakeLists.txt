set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/test)
set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
include_directories(${SRC_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(samplerTests EXCLUDE_FROM_ALL
    MultiLevelSamplerTest.cc
    CollectiveSectionMoverTest.cc
    CollectiveSectionSamplerTest.cc
    ${SRC_DIR}/sampler/SectionChooser.cc
    ${SRC_DIR}/sampler/DoubleSectionChooser.cc
    ${SRC_DIR}/sampler/CollectiveSectionMover.cc
    ${SRC_DIR}/sampler/CollectiveSectionSampler.cc
    ${SRC_DIR}/sampler/DoubleCollectiveSectionSampler.cc
    ${SRC_DIR}/BeadFactory.cc
    ${SRC_DIR}/util/SuperCell.cc
    ${SRC_DIR}/util/RandomNumGenerator.cc
    ${SRC_DIR}/util/Permutation.cc
    ${SRC_DIR}/stats/AccRejEstimator.cc
)
target_link_libraries(samplerTests gtest gtest_main -lpthread)
target_link_libraries(samplerTests ${BLAS_LIB})
target_link_libraries(samplerTests ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(samplerTests ${LIBF2C_LIB})
endif()
target_link_libraries(samplerTests ${GSL_LIB})

add_dependencies(unittest samplerTests)
