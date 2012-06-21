set (sources
    EMARateActionTest.cc
    EMARateMoverTest.cc
    EMARateEstimatorTest.cc
    EMARateTestBeadPositioner.cc
)

add_library(emarate_tests STATIC EXCLUDE_FROM_ALL ${sources})
