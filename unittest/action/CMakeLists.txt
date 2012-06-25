get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)
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
