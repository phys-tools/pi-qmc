set (sources
    EwaldSum.cc
    Hungarian.cc
    OptEwaldSum.cc
    PairDistance.cc
    TradEwaldSum.cc
    WireEwald.cc
)

add_library(util STATIC ${sources})
