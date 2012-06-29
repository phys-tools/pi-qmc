get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources}
    ${dir}/EstimatorParserTest.cc
    PARENT_SCOPE
)
