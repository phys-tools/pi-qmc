get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources}
    ${dir}/ScalarEstimatorTest.cpp
    ${dir}/SimpleScalarAccumulatorTest.cpp
    ${dir}/UnitsTest.cpp
    PARENT_SCOPE
)
