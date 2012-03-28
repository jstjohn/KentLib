include ../include.mk
DESTDIR=$(HOME)/
BINDIR=bin
CC=gcc
A=faNStatsOverNoCovBamRegions
O=$(patsubst %.c,%.o,$(wildcard *.c))
SOURCES=$(wildcard *.c)

all: ${A} ${O} ${SOURCES}

install: ${O} ${MYLIBS} ${SOURCES}
	${CC} ${USEROPTS} -o ${DESTDIR}${BINDIR}/${A} ${O} ${LDFLAGS}

${A}: ${O} ${MYLIBS} ${SOURCES}
	${CC} ${USEROPTS} -o ${A} ${O} ${LDFLAGS}

clean::
	rm -f ${A} ${O}

.c.o:
	$(CC) ${CFLAGS} ${USEROPTS} -c  $< -o $@

check-syntax:
	$(CC) ${CFLAGS} ${USEROPTS} -c -o .nul -S ${CHK_SOURCES}
