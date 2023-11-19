CFLAGS=-std=gnu99
LDLIBS=-lz -lm -pthread

ifneq ($(shell uname -s),Darwin)
	CFLAGS+=-march=native
endif

release: CFLAGS+=-O3
release: yamdedup yamscan yamshuf

debug: CFLAGS+=-g -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare -fsanitize=address,undefined -fno-omit-frame-pointer -DDEBUG -Wcast-qual
debug: yamdedup yamscan yamshuf

yamdedup: src/yamdedup.c
	mkdir -p bin ;\
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o bin/$@

yamscan: src/yamscan.c
	mkdir -p bin ;\
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o bin/$@

yamshuf: src/yamshuf.c
	mkdir -p bin ;\
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o bin/$@

