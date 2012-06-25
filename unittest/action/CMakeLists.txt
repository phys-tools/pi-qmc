set (dir action)
set_property(TARGET unittest APPEND PROPERTY sources
    ${test_sources}
    ${dir}/CoulombActionTest.cc
    ${dir}/PrimSHOActionTest.cc
    ${dir}/GaussianActionTest.cc
    ${dir}/SHOActionTest.cc
    ${dir}/GrapheneActionTest.cc
    ${dir}/DotGeomActionTest.cc
    ${dir}/GateActionTest.cc
    ${dir}/PrimAnisSHOActionTest.cc
    ${dir}/GaussianDotActionTest.cc
    ${dir}/PrimCosineActionTest.cc
    ${dir}/EFieldActionTest.cc
    ${dir}/PrimColloidalActionTest.cc
)
