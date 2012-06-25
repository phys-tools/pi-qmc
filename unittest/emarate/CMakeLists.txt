set (dir emarate)
set_property(TARGET unittest APPEND PROPERTY sources
    ${dir}/EMARateActionTest.cc
    ${dir}/EMARateMoverTest.cc
    ${dir}/EMARateEstimatorTest.cc
    ${dir}/EMARateTestBeadPositioner.cc
)
