set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/testbin)
set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
set(TEST_DIR ${PI_QMC_SOURCE_DIR}/unittest)
include_directories(${SRC_DIR} ${TEST_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(actionTests EXCLUDE_FROM_ALL
    CoulombActionTest.cc
    ${SRC_DIR}/action/CoulombAction.cc
    ${SRC_DIR}/action/PairAction.cc
    ${SRC_DIR}/action/PairIntegrator.cc
    ${SRC_DIR}/action/ImagePairAction.cc
    ${SRC_DIR}/action/EwaldImagePairAction.cc
    ${SRC_DIR}/util/EwaldSum.cc
    ${SRC_DIR}/util/TradEwaldSum.cc
    ${SRC_DIR}/util/OptEwaldSum.cc
    PrimSHOActionTest.cc
    ${SRC_DIR}/action/PrimSHOAction.cc
    GaussianActionTest.cc
    ${SRC_DIR}/action/GaussianAction.cc
    SHOActionTest.cc
    ${SRC_DIR}/action/SHOAction.cc
    GrapheneActionTest.cc
    ${SRC_DIR}/action/GrapheneAction.cc
    DotGeomActionTest.cc
    ${SRC_DIR}/action/DotGeomAction.cc
    GateActionTest.cc
    ${SRC_DIR}/action/GateAction.cc
    PrimAnisSHOActionTest.cc
    ${SRC_DIR}/action/PrimAnisSHOAction.cc
    GaussianDotActionTest.cc
    ${SRC_DIR}/action/GaussianDotAction.cc
    PrimCosineActionTest.cc
    ${SRC_DIR}/action/PrimCosineAction.cc
    EFieldActionTest.cc
    ${SRC_DIR}/action/EFieldAction.cc
    PrimColloidalActionTest.cc
    ${SRC_DIR}/action/PrimColloidalAction.cc
    ${SRC_DIR}/Species.cc
    ${SRC_DIR}/SimulationInfo.cc
    ${SRC_DIR}/util/SuperCell.cc
    ${SRC_DIR}/BeadFactory.cc
    ${TEST_DIR}/sampler/MultiLevelSamplerFake.cc
)
target_link_libraries(actionTests gtest gtest_main -lpthread)
target_link_libraries(actionTests ${BLAS_LIB})
target_link_libraries(actionTests ${LAPACK_LIB})
target_link_libraries(actionTests ${FFTW3_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(actionTests ${LIBF2C_LIB})
endif()
add_dependencies(unittest actionTests)
