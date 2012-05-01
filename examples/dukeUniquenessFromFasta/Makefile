L=-lm -lz
KENTLIBDIR=../../
KENTLIBS=${KENTLIBDIR}/libkent.a
KENTINC=${KENTLIBDIR}/inc
DESTDIR=$(HOME)/
BINDIR=bin
COPTS=-I${KENTINC}
CXX=g++
A=dukeUniquenessFromFasta
USEROPTS=-g
O=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
SOURCES=$(wildcard *.c)

all: ${A} ${O} ${SOURCES}

install: ${O} ${MYLIBS} ${SOURCES}
	${CXX} ${USEROPTS} -o ${DESTDIR}${BINDIR}/${A} ${O} ${KENTLIBS} ${L}

${A}: ${O} ${MYLIBS} ${SOURCES}
	${CXX} ${USEROPTS} -o ${A} ${O} ${KENTLIBS} ${L}

clean::
	rm -f ${A} ${O}

.cpp.o:
	$(CXX) ${COPTS} ${USEROPTS} -c  $< -o $@

check-syntax:
	$(CXX) ${COPTS} ${USEROPTS} -c -o .nul -S ${CHK_SOURCES}
