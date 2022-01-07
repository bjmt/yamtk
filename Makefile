all: prep minimotif

prep:
	mkdir -p bin

minimotif:
	cc -lm -std=gnu99 -pthread -O3 src/minimotif.c -o bin/minimotif
