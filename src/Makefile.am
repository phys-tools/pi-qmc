SUBDIRS = action demo emarate estimator fixednode parser \
	sampler spin stats util

bin_PROGRAMS = pi
pi_DEPENDENCIES = \
    sampler/libsampler.la \
    emarate/libemarate.la \
    estimator/libestimator.la \
    fixednode/libfixednode.la \
    action/libaction.la \
    demo/libdemo.la \
    spin/libspin.la \
    stats/libstats.la \
    parser/libparser.la \
    util/libutil.la
pi_LDADD = -Lparser -lparser -Laction -laction -Lutil -lutil\
	-Lsampler -lsampler -Lemarate -lemarate -Lfixednode -lfixednode \
	-Ldemo -ldemo -Lspin -lspin -Lstats -lstats -Lestimator -lestimator
pi_SOURCES = \
	AperiodicGaussian.cc \
	BeadFactory.cc \
	BinProbDensity.cc \
	Collect.cc \
	ConditionalDensityGrid.cc \
	CubicLattice.cc \
	DoubleParallelPaths.cc \
	EnumeratedModelState.cc \
	Measure.cc \
	ParallelPaths.cc \
	PathReader.cc \
	Paths.cc \
	PeriodicGaussian.cc \
	Permutation.cc \
	ProbDensityGrid.cc \
	RandomNumGenerator.cc \
	RingLattice.cc \
	SHOPhase.cc \
	SeedRandom.cc \
	SerialPaths.cc \
	SimInfoWriter.cc \
	SimulationInfo.cc \
	Species.cc \
	Spin4DPhase.cc \
	SpinModelState.cc \
	StructReader.cc \
	SuperCell.cc \
	WorkerShifter.cc \
	WritePaths.cc \
	WriteProbDensity.cc \
	XMLWriter.cc \
	main.cc
noinst_HEADERS = \
	Algorithm.h \
	AperiodicGaussian.h \
	BeadFactory.h \
	Beads.h \
	BinProbDensity.h \
	Collect.h \
	CompositeAlgorithm.h \
	ConditionalDensityGrid.h \
	CubicLattice.h \
	DoubleParallelPaths.h \
	EnumeratedModelState.h \
	Help.h \
	LinkSummable.h \
	Loop.h \
	Measure.h \
	ModelState.h \
	ParallelPaths.h \
	PathReader.h \
	Paths.h \
	PeriodicGaussian.h \
	Permutation.h \
	PhaseModel.h \
	Positioner.h \
	ProbDensityGrid.h \
	RandomNumGenerator.h \
	RingLattice.h \
	SHOPhase.h \
	SeedRandom.h \
	SerialPaths.h \
	SimInfoWriter.h \
	SimulationInfo.h \
	Species.h \
	Spin4DPhase.h \
	SpinModelState.h \
	StructReader.h \
	SuperCell.h \
	WorkerShifter.h \
	WritePaths.h \
	WriteProbDensity.h \
	XMLWriter.h