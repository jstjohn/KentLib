CC=gcc
LDFLAGS=
CFLAGS=-c -Wall -Iinc
LIBDIR=lib
INCDIR=inc
SOURCES=$(shell find $(LIBDIR) -type f -name '*.c')
HEADERS=$(shell find $(INCDIR) -type f -name '*.h')
OBJECTS=$(patsubst %.c,%.o,$(SOURCES))
LIBOUT=jkweb.a

all: $(SOURCES) $(LIBOUT) $(HEADERS)

$(LIBOUT): $(OBJECTS)
	ar rcus $(LIBOUT) $(OBJECTS) 

$(OBJECTS): $(HEADERS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	-rm $(OBJECTS) $(LIBOUT)

