DESTDIR=$(HOME)/
BINDIR=bin
CXX=g++
A=parsimonyFromMaf
USEROPTS=-g
O=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
SOURCES=$(wildcard *.c)
include ../include.mk

all: ${A} ${O} ${SOURCES}

install: ${O} ${MYLIBS} ${SOURCES}
	${CXX} ${USEROPTS}  ${CFLAGS} -o ${DESTDIR}${BINDIR}/${A} ${O} ${LDFLAGS}

${A}: ${O} ${MYLIBS} ${SOURCES}
	${CXX} ${USEROPTS} ${CFLAGS} -o ${A} ${O} ${KENTLIBS} ${L} ${LDFLAGS}

clean::
	rm -f ${A} ${O}

.cpp.o:
	$(CXX) ${CFLAGS} ${USEROPTS} -c  $< -o $@

check-syntax:
	$(CXX) ${CFLAGS} ${USEROPTS} -c -o .nul -S ${CHK_SOURCES}
