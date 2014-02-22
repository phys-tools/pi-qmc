get_filename_component(dir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

set(sources
    ${sources}
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
    ${dir}/coulomb/Coulomb1DLinkActionTest.cpp
    ${dir}/coulomb/Coulomb3DLinkActionTest.cpp
    ${dir}/coulomb/CoulombLinkActionTest.cpp
    ${dir}/interaction/InverseCosh2PotentialTest.cpp
    ${dir}/interaction/PairIntegratorTest.cpp
    ${dir}/interaction/SHOInteractionTest.cpp
    PARENT_SCOPE
)
