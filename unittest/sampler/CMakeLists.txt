set (sources
    MultiLevelSamplerTest.cc
    CollectiveSectionMoverTest.cc
    CollectiveSectionSamplerTest.cc
)

add_library(sampler_tests STATIC EXCLUDE_FROM_ALL ${sources})
