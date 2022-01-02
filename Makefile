all: build

build:
	mkdir -p bin;\
	cc -O2 src/minimotif.c -o bin/minimotif
