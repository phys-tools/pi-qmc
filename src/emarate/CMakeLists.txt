set (sources
    EMARateAction.cc
    EMARateEstimator.cc
    EMARateMover.cc
    EMARateWeight.cpp
)

add_library(emarate STATIC ${sources})
