#PATH to where RIsearch is to be found
MATPATH = \"/home/projects/users/anthon/projects/2024_risearch/temp_par/risearch/tables\"

CFLAGS += -Wall -pedantic -DMATPATH=$(MATPATH) -DRISVERSION=1

all: ../bin/RIsearch ../bin/RIsearch.dbg

../bin/RIsearch: *.c *.h
	$(CC) $(CFLAGS) -O3 $^ -lm -o $@

../bin/RIsearch.dbg: *.c *.h
	$(CC) $(CFLAGS) -g -O0 -DDEBUG -DVERBOSE=2 $^ -lm -o $@


.PHONY: clean

clean:
	rm -f ../bin/RIsearch ../bin/RIsearch.dbg *.o
