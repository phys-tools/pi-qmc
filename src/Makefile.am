SUBDIRS = action advancer base demo emarate estimator fixednode parser \
	spin stats util algorithm

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
pi_LDADD = -Laction -laction -Lalgorithm -lalgorithm -Ladvancer -ladvancer \
    -Lbase -lbase -Ldemo -ldemo -Lemarate -lemarate -Lestimator -lestimator \
    -Lfixednode -lfixednode -Lspin -lspin -Lstats -lstats -Lutil -lutil
pi_SOURCES = \
    main.cc
