bin_PROGRAMS = risearch2 risearch2_dbg

AM_CPPFLAGS = @OPENMP_CFLAGS@
AM_CFLAGS = @OPENMP_CFLAGS@
AM_LDFLAGS = @OPENMP_CFLAGS@

risearch2_SOURCES = dsm.c dsm.h dsm-others.c fasta.c fasta.h lists.c lists.h main.c main.h sa.c sa.h search.c search.h seed.c seed.h util.c weights.c weights.h
risearch2_CPPFLAGS = -DRISEARCH2 -DRISVERSION=2 -I../libdivsufsort-2.0.1/build/include
risearch2_CFLAGS = $(AM_CFLAGS) -O3 -Wall
risearch2_LDADD = ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a -lz -lpcre -lm

risearch2_dbg_SOURCES = $(risearch2_SOURCES)
risearch2_dbg_CPPFLAGS = $(risearch2_CPPFLAGS) -DDEBUG -DVERBOSE=2
risearch2_dbg_CFLAGS = $(risearch2_CFLAGS) -ggdb -g3 -O0 -Wall
risearch2_dbg_LDADD = $(risearch2_LDADD)

CLEANFILES = risearch2 risearch2_dbg

# Install binaries to bin directory
install-exec-hook:
	mkdir -p $(DESTDIR)/bin
	cp risearch2 risearch2_dbg $(DESTDIR)/bin
