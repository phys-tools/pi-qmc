SUBDIRS = action demo emarate estimator fixednode parser \
	sampler spin stats util algorithm

bin_PROGRAMS = pi
pi_DEPENDENCIES = \
    sampler/libsampler.la \
    algorithm/libalgorithm.la \
    emarate/libemarate.la \
    estimator/libestimator.la \
    fixednode/libfixednode.la \
    action/libaction.la \
    demo/libdemo.la \
    spin/libspin.la \
    stats/libstats.la \
    parser/libparser.la \
    util/libutil.la
pi_LDADD = -Lparser -lparser -Laction -laction -Lalgorithm -lalgorithm \
	-Lsampler -lsampler -Lemarate -lemarate -Lfixednode -lfixednode \
	-Ldemo -ldemo -Lspin -lspin -Lstats -lstats -Lestimator -lestimator \
        -Lutil -lutil
pi_SOURCES = \
	BeadFactory.cc \
	DoubleParallelPaths.cc \
	EnumeratedModelState.cc \
	ParallelPaths.cc \
	Paths.cc \
	SHOPhase.cc \
	SerialPaths.cc \
	SimInfoWriter.cc \
	SimulationInfo.cc \
	Species.cc \
	Spin4DPhase.cc \
	SpinModelState.cc \
	XMLWriter.cc \
	main.cc
noinst_HEADERS = \
	BeadFactory.h \
	Beads.h \
	DoubleParallelPaths.h \
	EnumeratedModelState.h \
	Help.h \
	LinkSummable.h \
	ModelState.h \
	ParallelPaths.h \
	Paths.h \
	PhaseModel.h \
	SHOPhase.h \
	SerialPaths.h \
	SimInfoWriter.h \
	SimulationInfo.h \
	Species.h \
	Spin4DPhase.h \
	SpinModelState.h \
	XMLWriter.h
