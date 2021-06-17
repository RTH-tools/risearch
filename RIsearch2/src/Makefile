CFLAGS += -fopenmp -Wall -I../libdivsufsort-2.0.1/build/include -DRISEARCH2

all: ../bin/risearch2.x
# ../bin/risearch2.dbg.x

../bin/risearch2.x: *.c
	$(CC) $(CFLAGS) -O3 $^ ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a -lz -o $@

../bin/risearch2.dbg.x: *.c
#	$(CC) $(CFLAGS) -g -O0 -DDEBUG -DVERBOSE=2 $^ ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a -lz -o $@
	$(CC) $(CFLAGS) -ggdb -g3 -O0 -DDEBUG $^ ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a -lz -o $@

.PHONY: clean

conv.x: dsm*
	$(CC) -O3 $^ -o $@

clean:
	rm -f ../bin/risearch2.x ../bin/risearch2.dbg.x conv.x
