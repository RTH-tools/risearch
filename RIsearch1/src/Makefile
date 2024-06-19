CFLAGS += -Wall -pedantic

all: ../bin/RIsearch ../bin/RIsearch.dbg

../bin/RIsearch: *.c *.h
	$(CC) $(CFLAGS) -O3 $^ -DRISVERSION=1 -lm -o $@

../bin/RIsearch.dbg: *.c *.h
	$(CC) $(CFLAGS) -g -O0 -DRISVERSION=1 -DDEBUG -DVERBOSE=2 $^ -lm -o $@


.PHONY: clean

clean:
	rm -f ../bin/RIsearch ../bin/RIsearch.dbg *.o
