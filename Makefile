all: build

build:
	mkdir -p bin;\
	cc -pthread -O3 src/minimotif.c -o bin/minimotif
