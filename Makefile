CC=gcc
MACHTYPE=x86_64
LDFLAGS=-Lthirdparty/samtools -lbam -pthread
CFLAGS=-c -fPIC -Wall -Iinc -Ithirdparty/samtools -DUSE_BAM=1 -DMACHTYPE_$(MACHTYPE)
LIBDIR=lib
INCDIR=inc
SOURCES=$(shell find $(LIBDIR) -type f -name '*.c')
HEADERS=$(shell find $(INCDIR) -type f -name '*.h')
OBJECTS=$(patsubst %.c,%.o,$(SOURCES))
LIBOUT=libkent.a
LEGACYOUT=jkweb.a
SHAREDOUT=libkent.so

all: $(SOURCES) $(LIBOUT) $(HEADERS) THIRDPARTY

THIRDPARTY: thirdparty/samtools/libbam.a

thirdparty/samtools/libbam.a:
	cd thirdparty/samtools && make && cd ../..

$(LIBOUT): $(OBJECTS)
	ar rcus $(LIBOUT) $(OBJECTS)
	ln -sf ${LIBOUT} ${LEGACYOUT}

$(OBJECTS): $(HEADERS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	-rm $(OBJECTS) $(LIBOUT) ${LEGACYOUT}
	cd thirdparty/samtools && make clean && cd ../..

#not working yet
shared: ${LIBOUT}
	${CC} -shared -L. -Iinc -o ${SHAREDOUT} ${OBJECTS} -lz -lm -lc
shared-darwin: ${LIBOUT}
	${CC} -shared -L. -Iinc -Wl,-install_name,${SHAREDOUT} -o ${SHAREDOUT} ${OBJECTS} -lz -lm -lc -lkent


