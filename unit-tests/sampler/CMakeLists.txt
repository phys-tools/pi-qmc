get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set_property(TARGET unittest APPEND PROPERTY sources
    ${dir}/MultiLevelSamplerTest.cc
    ${dir}/CollectiveSectionMoverTest.cc
    ${dir}/CollectiveSectionSamplerTest.cc
    ${dir}/MultiLevelSamplerFake.cc
)
