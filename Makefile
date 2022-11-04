CFLAGS=-std=gnu99
LDLIBS=-lz -lm -pthread

ifneq ($(shell uname -s),Darwin)
	CFLAGS+=-march=native
endif

all: CFLAGS+=-O3
all: minidedup minimotif

debug: CFLAGS+=-g -Wall -Wextra -Wno-sign-compare
debug: minidedup minimotif

minimotif: src/minimotif.c
	mkdir -p bin ;\
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o bin/$@

minidedup: src/minidedup.c
	mkdir -p bin ;\
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o bin/$@

