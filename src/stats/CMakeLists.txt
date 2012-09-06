set (sources
    AccRejEstimator.cc
    ArrayEstimator.cc
    AsciiReportBuilder.cc
    EstimatorReportBuilder.cpp
    EstimatorManager.cc
    H5ReportBuilder.cc
    H5ScalarReportWriter.cpp
    MPIManager.cc
    ReportWriters.cpp
    ScalarEstimator.cc
    ScalarReportWriter.cpp
    StdoutReportBuilder.cc
    Units.cc
)

add_library(stats STATIC ${sources})
