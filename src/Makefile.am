SUBDIRS = action advancer base demo emarate estimator fixednode parser \
	spin stats util algorithm

pi_CXXFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_srcdir)/contrib/blitz-0.9


bin_PROGRAMS = pi
pi_DEPENDENCIES = \
    action/libaction.la \
    advancer/libadvancer.la \
    algorithm/libalgorithm.la \
    base/libbase.la \
    demo/libdemo.la \
    emarate/libemarate.la \
    estimator/libestimator.la \
    fixednode/libfixednode.la \
    parser/libparser.la \
    spin/libspin.la \
    stats/libstats.la \
    util/libutil.la
pi_LDADD = -Lparser -lparser -Lemarate -lemarate -Laction -laction \
    -Lalgorithm -lalgorithm -Ladvancer -ladvancer -Ldemo -ldemo \
    -Lestimator -lestimator -Lfixednode -lfixednode \
    -Lspin -lspin -Lstats -lstats -Lbase -lbase -Lutil -lutil
pi_SOURCES = \
    main.cc
