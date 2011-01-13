C=gcc
MACHTYPE=x86_64
LDFLAGS=
CFLAGS=-c -Wall -Iinc -DMACHTYPE_$(MACHTYPE)
LIBDIR=lib
INCDIR=inc
SOURCES=$(shell find $(LIBDIR) -type f -name '*.c')
HEADERS=$(shell find $(INCDIR) -type f -name '*.h')
OBJECTS=$(patsubst %.c,%.o,$(SOURCES))
LIBOUT=libkent.a
LEGACYOUT=jkweb.a

all: $(SOURCES) $(LIBOUT) $(HEADERS)

$(LIBOUT): $(OBJECTS)
	ar rcus $(LIBOUT) $(OBJECTS)
	ln -s ${LIBOUT} ${LEGACYOUT}

$(OBJECTS): $(HEADERS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	-rm $(OBJECTS) $(LIBOUT) ${LEGACYOUT}

