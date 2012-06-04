set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/test)
set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
include_directories(${SRC_DIR} ${GTEST_SOURCE_DIR}/include)

add_executable(actionTests EXCLUDE_FROM_ALL
    CoulombActionTest.cc
    ../CoulombAction.cc
    ../PairAction.cc
    ../PairIntegrator.cc
    ../ImagePairAction.cc
    ../EwaldImagePairAction.cc
    ${SRC_DIR}/util/EwaldSum.cc
    ${SRC_DIR}/util/TradEwaldSum.cc
    ${SRC_DIR}/util/OptEwaldSum.cc
    PrimSHOActionTest.cc
    ../PrimSHOAction.cc
    GaussianActionTest.cc
    ../GaussianAction.cc
    SHOActionTest.cc
    ../SHOAction.cc
    GrapheneActionTest.cc
    ../GrapheneAction.cc
    DotGeomActionTest.cc
    ../DotGeomAction.cc
    GateActionTest.cc
    ../GateAction.cc
    PrimAnisSHOActionTest.cc
    ../PrimAnisSHOAction.cc
    GaussianDotActionTest.cc
    ../GaussianDotAction.cc
    PrimCosineActionTest.cc
    ../PrimCosineAction.cc
    EFieldActionTest.cc
    ../EFieldAction.cc
    PrimColloidalActionTest.cc
    ../PrimColloidalAction.cc
    ${SRC_DIR}/Species.cc
    ${SRC_DIR}/SimulationInfo.cc
    ${SRC_DIR}/util/SuperCell.cc
    ${SRC_DIR}/BeadFactory.cc
    ${SRC_DIR}/sampler/test/MultiLevelSamplerFake.cc
)
target_link_libraries(actionTests gtest gtest_main -lpthread)
target_link_libraries(actionTests ${BLAS_LIB})
target_link_libraries(actionTests ${LAPACK_LIB})
target_link_libraries(actionTests ${FFTW3_LIB})
if (EXISTS ${LIBF2C_LIB})
    target_link_libraries(actionTests ${LIBF2C_LIB})
endif()
add_dependencies(unittest actionTests)
