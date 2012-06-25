set (dir sampler)
set_property(TARGET unittest APPEND PROPERTY sources
    ${dir}/MultiLevelSamplerTest.cc
    ${dir}/CollectiveSectionMoverTest.cc
    ${dir}/CollectiveSectionSamplerTest.cc
    ${dir}/MultiLevelSamplerFake.cc
)
