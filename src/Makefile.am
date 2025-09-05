bin_PROGRAMS = RIsearch1 RIsearch1_dbg RIsearch2 RIsearch2_dbg

RIsearch1_SOURCES = dsm.c dsm.h fasta.c fasta.h main.h risearch.c weights.c util.c
RIsearch1_CFLAGS = -O3 -Wall -DMATPATH='"$(pkgdatadir)/tables"'
RIsearch1_LDADD = -lm -lz

RIsearch1_dbg_SOURCES = $(RIsearch1_SOURCES)
RIsearch1_dbg_CFLAGS = $(RIsearch1_CFLAGS) -ggdb -O0 -DDEBUG -DVERBOSE=2 
RIsearch1_dbg_LDADD = $(RIsearch1_LDADD)

RIsearch2_SOURCES = dsm.c dsm.h fasta.c fasta.h lists.c lists.h main.c main.h sa.c sa.h search.c search.h seed.c seed.h util.c weights.c weights.h risearch2.h
RIsearch2_CFLAGS = @OPENMP_CFLAGS@ -O3 -Wall -I../libdivsufsort-2.0.1/build/include -DMATPATH='"$(pkgdatadir)/tables"'
RIsearch2_LDADD = @OPENMP_CFLAGS@ ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a -lz -lpcre -lm

RIsearch2_dbg_SOURCES = $(RIsearch2_SOURCES)
RIsearch2_dbg_CFLAGS = $(RIsearch2_CFLAGS) -ggdb -g3 -O0 -Wall -DDEBUG -DVERBOSE=2
RIsearch2_dbg_LDADD = $(RIsearch2_LDADD)


CLEANFILES = RIsearch RIsearch_dbg RIsearch2 RIsearch2_dbg
