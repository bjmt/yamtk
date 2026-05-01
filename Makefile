#    yamtk: Yet Another Motif Toolkit
#    Copyright (C) 2025  Benjamin Jean-Marie Tremblay
#  
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#  
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#  
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#  

CC ?=cc
CFLAGS=-std=gnu99
LDFLAGS=-std=gnu99
LDLIBS=-lz -lm -pthread
PREFIX  ?=/usr/local
BINDIR  ?=bin

release: CFLAGS+=-O3
release: LDFLAGS+=-O3
release: yamtk

debug: CFLAGS+=-g -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer -DDEBUG -Wcast-qual
debug: LDFLAGS+=-g -Og -fsanitize=address,undefined
debug: yamtk

clean:
	-rm -f ./src/*.o
	-rm -f ./yamtk

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $^ -o $@

objects := $(patsubst %.c,%.o,$(wildcard src/*.c))

yamtk: $(objects)
	$(CC) $(LDFLAGS) $(objects) -o $@ $(LDLIBS)

install: yamtk
	install -p ./yamtk $(PREFIX)/$(BINDIR)

uninstall:
	-rm -f $(PREFIX)/$(BINDIR)/yamtk

check: yamtk
	@bash test/run_tests.sh

check-debug:
	$(MAKE) clean
	$(MAKE) debug
	@bash test/run_tests.sh

.PHONY: check check-debug

