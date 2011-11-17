cmake_minimum_required(VERSION 2.6.2)
project(PITEST)

add_definitions(-DNDIM=3)

include_directories(${PITEST_SOURCE_DIR}/../../src)
set(SOURCE_DIR ${PITEST_SOURCE_DIR}/../../src)

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
set(GMOCK_BOTH_LIBRARIES gmock gmock_main)

add_executable(pitest 
    EMARateActionTest
    ${SOURCE_DIR}/EMARateAction
    ${SOURCE_DIR}/Species
    ${SOURCE_DIR}/SuperCell
    ${SOURCE_DIR}/SimulationInfo
    ${SOURCE_DIR}/BeadFactory
)
target_link_libraries(pitest ${GMOCK_BOTH_LIBRARIES})
target_link_libraries(pitest ${GTEST_BOTH_LIBRARIES})
target_link_libraries(pitest pthread)
add_test(pitest pitest)
