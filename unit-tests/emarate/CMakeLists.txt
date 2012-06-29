get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources}
    ${dir}/EMARateActionTest.cc
    ${dir}/EMARateMoverTest.cc
    ${dir}/EMARateEstimatorTest.cc
    ${dir}/EMARateTestBeadPositioner.cc
    PARENT_SCOPE
)
