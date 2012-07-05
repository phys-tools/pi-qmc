SUBDIRS = action advancer demo emarate estimator fixednode parser \
	spin stats util algorithm

bin_PROGRAMS = pi
pi_DEPENDENCIES = \
    advancer/libadvancer.la \
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
	-Ladvancer -ladvancer -Lemarate -lemarate -Lfixednode -lfixednode \
	-Ldemo -ldemo -Lspin -lspin -Lstats -lstats -Lestimator -lestimator \
        -Lutil -lutil
pi_SOURCES = \
	BeadFactory.cc \
	DoubleParallelPaths.cc \
	EnumeratedModelState.cc \
	ParallelPaths.cc \
	Paths.cc \
	SerialPaths.cc \
	SimInfoWriter.cc \
	SimulationInfo.cc \
	Species.cc \
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
	SerialPaths.h \
	SimInfoWriter.h \
	SimulationInfo.h \
	Species.h \
	SpinModelState.h \
	XMLWriter.h
