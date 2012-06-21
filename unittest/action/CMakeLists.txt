set (sources
    CoulombActionTest.cc
    PrimSHOActionTest.cc
    GaussianActionTest.cc
    SHOActionTest.cc
    GrapheneActionTest.cc
    DotGeomActionTest.cc
    GateActionTest.cc
    PrimAnisSHOActionTest.cc
    GaussianDotActionTest.cc
    PrimCosineActionTest.cc
    EFieldActionTest.cc
    PrimColloidalActionTest.cc
)

add_library(action_tests STATIC EXCLUDE_FROM_ALL ${sources})
