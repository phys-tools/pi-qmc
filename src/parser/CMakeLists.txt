set (sources
    ActionParser.cc
    EstimatorParser.cc
    MainParser.cc
    PIMCParser.cc
    SimInfoParser.cc
    XMLParser.cc
    XMLUnitParser.cc
)

add_library(parser STATIC ${sources})

add_subdirectory(test)
