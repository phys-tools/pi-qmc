set(GTEST_FOUND "")
find_package(GTest QUIET)

if(GTEST_FOUND)
    message("GTEST found, unit testing is enabled") 
    set(EXECUTABLE_OUTPUT_PATH ${PI_QMC_BINARY_DIR}/test)
    include_directories(${GTEST_INCLUDE_DIRS})

    FIND_LIBRARY(GMOCK_LIBRARIES gmock
        HINTS
            $ENV{GMOCK_ROOT}
            ${GMOCK_ROOT}
            $ENV{GMOCK_ROOT}/make
            ${GMOCK_ROOT}/make
    )
    if(NOT DEFINED GMOCK_LIBRARIES)
        MESSAGE(SEND_ERROR "Google mock library path not Set to GMOCK_LIBRARIES")
    endif()

    find_path(GMOCK_INCLUDE_DIR gmock/gmock.h
        HINTS
            $ENV{GMOCK_ROOT}/include
            ${GMOCK_ROOT}/include
    )
    if(NOT DEFINED GMOCK_INCLUDE_DIR)
        MESSAGE(SEND_ERROR "Google mock include path not set to GMOCK_INCLUDE_DIR") 
    endif() 

    set(SRC_DIR ${PI_QMC_SOURCE_DIR}/src)
    include_directories(${SRC_DIR})

    add_executable(actionTests
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
        ${SRC_DIR}/Species.cc
        ${SRC_DIR}/SimulationInfo.cc
        ${SRC_DIR}/util/SuperCell.cc
        ${SRC_DIR}/BeadFactory.cc
        ${SRC_DIR}/sampler/test/MultiLevelSamplerFake.cc
    )
    target_link_libraries(actionTests ${GTEST_BOTH_LIBRARIES} -lpthread)
  
    target_link_libraries(actionTests ${BLAS_LIB})
    target_link_libraries(actionTests ${LAPACK_LIB})
    if (EXISTS ${LIBF2C_LIB})
        target_link_libraries(actionTests ${LIBF2C_LIB})
    endif()

    else()
    message(AUTHOR_WARNING "** GTEST not found, unit testing not enabled") 
endif(GTEST_FOUND)
