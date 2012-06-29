get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources}
    ${dir}/MultiLevelSamplerTest.cc
    ${dir}/CollectiveSectionMoverTest.cc
    ${dir}/CollectiveSectionSamplerTest.cc
    ${dir}/MultiLevelSamplerFake.cc
    PARENT_SCOPE
)
