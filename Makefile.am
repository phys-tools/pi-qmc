SUBDIRS = src doc
ACLOCAL_AMFLAGS = -I m4
docs:
	(cd doc; make docs)
.PHONY: test

test:
	bin/system-test-pi-qmc
