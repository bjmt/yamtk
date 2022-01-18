CFLAGS=-std=gnu99 -g -O3 -Ikseq -Wall -Wextra -Wno-sign-compare
LDLIBS=-lz -lm -pthread

all: minimotif clean

minimotif: src/minimotif.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

clean:
	mkdir -p bin ; mv minimotif bin/minimotif
