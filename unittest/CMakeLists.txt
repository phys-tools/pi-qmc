#add_subdirectory(demo)
#add_subdirectory(algorithm)
add_subdirectory(action)
add_subdirectory(parser)
add_subdirectory(sampler)
add_subdirectory(emarate)
#add_subdirectory(estimator)
add_subdirectory(fixednode)
#add_subdirectory(spin)
#add_subdirectory(stats)
add_subdirectory(util)

set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
set(TEST_DIR ${PI_QMC_SOURCE_DIR}/unittest)

include_directories(${SRC_DIR} ${TEST_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(unittest_pi EXCLUDE_FROM_ALL 
    unittest_main.cc
    action/CoulombActionTest.cc
    action/PrimSHOActionTest.cc
    action/GaussianActionTest.cc
    action/SHOActionTest.cc
    action/GrapheneActionTest.cc
    action/DotGeomActionTest.cc
    action/GateActionTest.cc
    action/PrimAnisSHOActionTest.cc
    action/GaussianDotActionTest.cc
    action/PrimCosineActionTest.cc
    action/EFieldActionTest.cc
    action/PrimColloidalActionTest.cc
    emarate/EMARateActionTest.cc
    emarate/EMARateMoverTest.cc
    emarate/EMARateEstimatorTest.cc
    emarate/EMARateTestBeadPositioner.cc
    fixednode/AugmentedNodesTest.cc
    fixednode/Atomic1sDMTest.cc
    fixednode/Atomic2spDMTest.cc
    parser/EstimatorParserTest.cc
    sampler/MultiLevelSamplerFake.cc
    sampler/MultiLevelSamplerTest.cc
    sampler/CollectiveSectionMoverTest.cc
    sampler/CollectiveSectionSamplerTest.cc
    util/AperiodicGaussianTest.cc
    util/HungarianTest.cc
    util/PeriodicGaussianTest.cc
    util/PermutationTest.cc
    util/SuperCellTest.cc

    ${SRC_DIR}/BeadFactory.cc
    ${SRC_DIR}/Paths.cc
    ${SRC_DIR}/SerialPaths.cc
    ${SRC_DIR}/SimulationInfo.cc
    ${SRC_DIR}/Species.cc
)

target_link_libraries(unittest_pi gtest -lpthread)
target_link_libraries(unittest_pi ${BLAS_LIB})
target_link_libraries(unittest_pi ${LAPACK_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(unittest_pi ${LIBF2C_LIB})
endif()
target_link_libraries(unittest_pi ${GSL_LIB})

target_link_libraries(unittest_pi ${FFTW3_LIB})

target_link_libraries(unittest_pi 
    demo parser sampler estimator emarate fixednode 
    spin stats action algorithm util)
    
add_dependencies(unittest unittest_pi)
