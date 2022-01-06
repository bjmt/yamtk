all: build

build:
	mkdir -p bin;\
	cc -pthread -O2 src/minimotif.c -o bin/minimotif
