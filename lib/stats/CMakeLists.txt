set (sources
    AccRejEstimator.cc
    ArrayBlockedEstimator.cc
    ArrayEstimator.cc
    AsciiReportBuilder.cc
    EstimatorManager.cc
    H5ReportBuilder.cc
    MPIManager.cc
    ScalarEstimator.cc
    StdoutReportBuilder.cc
    Units.cc
)

add_library(stats STATIC ${sources})
