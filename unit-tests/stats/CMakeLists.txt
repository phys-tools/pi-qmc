get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources}
    ${dir}/ScalarAccumulatorTest.cpp
    ${dir}/ScalarEstimatorTest.cpp
    ${dir}/UnitsTest.cpp
    PARENT_SCOPE
)
