CC=gcc
MACHTYPE=x86_64
LDFLAGS=
CFLAGS=-c -fPIC -Wall -Iinc -DMACHTYPE_$(MACHTYPE)
LIBDIR=lib
INCDIR=inc
SOURCES=$(shell find $(LIBDIR) -type f -name '*.c')
HEADERS=$(shell find $(INCDIR) -type f -name '*.h')
OBJECTS=$(patsubst %.c,%.o,$(SOURCES))
LIBOUT=libkent.a
LEGACYOUT=jkweb.a
SHAREDOUT=libkent.so

all: $(SOURCES) $(LIBOUT) $(HEADERS)

$(LIBOUT): $(OBJECTS)
	ar rcus $(LIBOUT) $(OBJECTS)
	ln -sf ${LIBOUT} ${LEGACYOUT}

$(OBJECTS): $(HEADERS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	-rm $(OBJECTS) $(LIBOUT) ${LEGACYOUT}

#not working yet
shared: ${LIBOUT}
	${CC} -shared -L. -Iinc -o ${SHAREDOUT} ${OBJECTS} -lz -lm -lc
shared-darwin: ${LIBOUT}
	${CC} -shared -L. -Iinc -Wl,-install_name,${SHAREDOUT} -o ${SHAREDOUT} ${OBJECTS} -lz -lm -lc -lkent


