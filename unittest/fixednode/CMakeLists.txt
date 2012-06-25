set (dir fixednode)
set_property(TARGET unittest APPEND PROPERTY sources
    ${dir}/AugmentedNodesTest.cc
    ${dir}/Atomic1sDMTest.cc
    ${dir}/Atomic2spDMTest.cc
)
