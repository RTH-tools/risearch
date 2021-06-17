CFLAGS += -Wall -pedantic -ansi

all: RISEARCH # debug

RISEARCH: risearch.c fasta.c dsm.c weights.c
	$(CC) $(CFLAGS) -O3 $^ -o RIsearch

debug: risearch.c fasta.c dsm.c weights.c
	$(CC) $(CFLAGS) -g -O0 -DDEBUG -DVERBOSE=2 $^ -o RIs-dbg


.PHONY: clean

clean:
	rm -f RIsearch RIs-dbg *.o
