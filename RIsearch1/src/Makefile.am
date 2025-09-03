bin_PROGRAMS = RIsearch RIsearch_dbg

RIsearch_SOURCES = dsm.c dsm.h fasta.c fasta.h main.h risearch.c weights.c
RIsearch_CPPFLAGS = -DRISVERSION=1
RIsearch_CFLAGS = -O3 -Wall -pedantic
RIsearch_LDADD = -lm

RIsearch_dbg_SOURCES = $(RIsearch_SOURCES)
RIsearch_dbg_CPPFLAGS = $(RIsearch_CPPFLAGS) -DDEBUG -DVERBOSE=2
RIsearch_dbg_CFLAGS = -g -O0 -Wall -pedantic
RIsearch_dbg_LDADD = -lm

CLEANFILES = RIsearch RIsearch_dbg

# Install binaries to bin directory
install-exec-hook:
	mkdir -p $(DESTDIR)/bin
	cp RIsearch RIsearch_dbg $(DESTDIR)/bin
