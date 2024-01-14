# This makefile can be used to build the logCompensatedOctaveSmooth
# program.  To build it, use the command 'make logCompensatedOctaveSmooth'.

VERSION=10.0
MLINKDIR = /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64/
SYS = MacOSX-x86-64
CADDSDIR = ${MLINKDIR}/CompilerAdditions
EXTRA_CFLAGS=-stdlib=libstdc++

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

MPREP = ${CADDSDIR}/mprep
RM = rm

CC = /usr/bin/clang
CXX = /usr/bin/clang++

all : logCompensatedOctaveSmooth clean

logCompensatedOctaveSmooth : logCompensatedOctaveSmoothtm.o logCompensatedOctaveSmooth.o
	${CXX} ${EXTRA_CFLAGS} -I${INCDIR} logCompensatedOctaveSmoothtm.o logCompensatedOctaveSmooth.o -L${LIBDIR} -lMLi4 -lstdc++ -framework Foundation -o $@

.c.o :
	${CC} -c ${EXTRA_CFLAGS} -I${INCDIR} $<

logCompensatedOctaveSmoothtm.c : logCompensatedOctaveSmooth.tm
	${MPREP} $? -o $@

clean :
	@ ${RM} -rf *.o *tm.c
