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

    add_executable(samplerTests
        MultiLevelSamplerTest.cc
        CollectiveSectionMoverTest.cc
        ${SRC_DIR}/sampler/CollectiveSectionMover.cc
    )
    
    target_link_libraries(samplerTests ${GTEST_BOTH_LIBRARIES} -lpthread)
  
    target_link_libraries(samplerTests ${BLAS_LIB})
    target_link_libraries(samplerTests ${LAPACK_LIB})
    if (EXISTS ${LIBF2C_LIB})
        target_link_libraries(samplerTests ${LIBF2C_LIB})
    endif()

    else()
    message(AUTHOR_WARNING "** GTEST not found, unit testing not enabled") 
endif(GTEST_FOUND)
