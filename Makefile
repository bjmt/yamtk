CFLAGS=-std=gnu99
LDLIBS=-lz -lm -pthread

ifneq ($(shell uname -s),Darwin)
	CFLAGS+=-march=native
endif

all: CFLAGS+=-O3
all: yamdedup yamscan

debug: CFLAGS+=-g -Og -Wall -Wextra -Wno-sign-compare -fsanitize=address,undefined -fno-omit-frame-pointer -DDEBUG
debug: yamdedup yamscan

yamscan: src/yamscan.c
	mkdir -p bin ;\
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o bin/$@

yamdedup: src/yamdedup.c
	mkdir -p bin ;\
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o bin/$@

