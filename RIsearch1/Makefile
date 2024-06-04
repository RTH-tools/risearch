CFLAGS += -Wall -pedantic

all: RIsearch RIsearch.dbg

RIsearch: *.c *.h
	$(CC) $(CFLAGS) -O3 $^ -DRISVERSION=1 -lm -o $@

RIsearch.dbg: *.c *.h
	$(CC) $(CFLAGS) -g -O0 -DRISVERSION=1 -DDEBUG -DVERBOSE=2 $^ -lm -o $@


.PHONY: clean

clean:
	rm -f RIsearch RIs-dbg *.o
