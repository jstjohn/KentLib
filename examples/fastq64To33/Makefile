L=-lm -lz
KENTLIBDIR=../../
KENTLIBS=${KENTLIBDIR}/libkent.a
KENTINC=${KENTLIBDIR}/inc
DESTDIR=$(HOME)/
BINDIR=bin
COPTS=-I${KENTINC}
CC=gcc
A=fastq64To33
USEROPTS=
O=$(patsubst %.c,%.o,$(wildcard *.c))
SOURCES=$(wildcard *.c)

all: ${A} ${O} ${SOURCES}

install: ${O} ${MYLIBS} ${SOURCES}
	${CC} ${USEROPTS} -o ${DESTDIR}${BINDIR}/${A} ${O} ${KENTLIBS} ${L}

${A}: ${O} ${MYLIBS} ${SOURCES}
	${CC} ${USEROPTS} -o ${A} ${O} ${KENTLIBS} ${L}

clean::
	rm -f ${A} ${O}

.c.o:
	$(CC) ${COPTS} ${USEROPTS} -c  $< -o $@

check-syntax:
	$(CC) ${COPTS} ${USEROPTS} -c -o .nul -S ${CHK_SOURCES}
