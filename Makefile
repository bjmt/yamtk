CFLAGS=-std=gnu99
LDLIBS=-lz -lm -pthread

ifneq ($(shell uname -s),Darwin)
	CFLAGS+=-march=native
endif

release: CFLAGS+=-O3
release: yamdedup yamscan yamshuf

debug: CFLAGS+=-g -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare -fsanitize=address,undefined -fno-omit-frame-pointer -DDEBUG -Wcast-qual
debug: yamdedup yamscan yamshuf

clean:
	-rm -f ./src/*.o
	-rm -f ./yamdedup ./yamscan ./yamshuf

yamdedup: src/yamdedup.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

yamscan: src/yamscan.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

yamshuf: src/yamshuf.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

