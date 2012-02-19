L=-lm -lz
KENTLIBDIR=../../
KENTLIBS=${KENTLIBDIR}/libkent.a
KENTINC=${KENTLIBDIR}/inc
SAMTOOLSDIR=$(HOME)/src/samtools
SAMTOOLS=${SAMTOOLSDIR}/libbam.a
DESTDIR=$(HOME)/
BINDIR=bin
COPTS=-I${KENTINC}
CC=gcc
A=readIDsFromBam
USEROPTS=-DHASH_FUNCTION=HASH_MUR -DHASH_USING_NO_STRICT_ALIASING -fno-strict-aliasing -DHASH_BLOOM=27
O=$(patsubst %.c,%.o,$(wildcard *.c))
SOURCES=$(wildcard *.c)

all: ${A} ${O} ${SOURCES}

install: ${O} ${MYLIBS} ${SOURCES}
	${CC} ${USEROPTS} -o ${DESTDIR}${BINDIR}/${A} ${O} ${KENTLIBS} ${SAMTOOLS} ${L}

${A}: ${O} ${MYLIBS} ${SOURCES}
	${CC} ${USEROPTS} -o ${A} ${O} ${KENTLIBS} ${L}

clean::
	rm -f ${A} ${O}

.c.o:
	$(CC) ${COPTS} ${USEROPTS} -c  $< -o $@

check-syntax:
	$(CC) ${COPTS} ${USEROPTS} -c -o .nul -S ${CHK_SOURCES}
