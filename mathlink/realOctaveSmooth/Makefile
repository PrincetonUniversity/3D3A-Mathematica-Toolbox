# This makefile can be used to build the realOctaveSmooth
# program.  To build it, use the command 'make realOctaveSmooth'.

VERSION=10.0
MLINKDIR = /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64
SYS = MacOSX-x86-64
CADDSDIR = ${MLINKDIR}/CompilerAdditions
EXTRA_CFLAGS=-stdlib=libstdc++ -stdlib=libc++

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

MPREP = ${CADDSDIR}/mprep
RM = rm

all : realOctaveSmooth clean

realOctaveSmooth : realOctaveSmoothtm.o realOctaveSmooth.o
	${CXX} ${EXTRA_CFLAGS} -I${INCDIR} realOctaveSmoothtm.o realOctaveSmooth.o -L${LIBDIR} -lMLi4 -lc++ -framework Foundation -o $@

.c.o :
	${CC} -c ${EXTRA_CFLAGS} -I${INCDIR} $<

realOctaveSmoothtm.c : realOctaveSmooth.tm
	${MPREP} $? -o $@

clean :
	@ ${RM} -rf *.o *tm.c
