CFLAGS+=-std=gnu99 -O3 -Ikseq
LDLIBS+=-lz -lm -pthread

all: minimotif clean

minimotif: src/minimotif_kseq.c
	$(CC) $(CFLAGS) $(LDLIBS) $(LDFLAGS) $^ -o $@

clean:
	mkdir -p bin ; mv minimotif src/minimotif
