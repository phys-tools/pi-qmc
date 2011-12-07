include_directories(${PI_SOURCE_DIR}/src)
set(SOURCE_DIR ${PI_SOURCE_DIR}/src)

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
