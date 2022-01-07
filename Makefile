all: build

build:
	mkdir -p bin;\
	cc -lm -std=gnu99 -pthread -O3 src/minimotif.c -o bin/minimotif
