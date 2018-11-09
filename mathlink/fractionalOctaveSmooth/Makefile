# This makefile can be used to build the fractionalOctaveSmooth
# program.  To build it, use the command 'make fractionalOctaveSmooth'.

VERSION=10.0
MLINKDIR = /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64
SYS = MacOSX-x86-64
CADDSDIR = ${MLINKDIR}/CompilerAdditions
EXTRA_CFLAGS=-stdlib=libstdc++ -stdlib=libc++

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

MPREP = ${CADDSDIR}/mprep
RM = rm

all : fractionalOctaveSmooth clean

fractionalOctaveSmooth : fractionalOctaveSmoothtm.o fractionalOctaveSmooth.o
	${CXX} ${EXTRA_CFLAGS} -I${INCDIR} fractionalOctaveSmoothtm.o fractionalOctaveSmooth.o -L${LIBDIR} -lMLi4 -lc++ -framework Foundation -o $@

.c.o :
	${CC} -c ${EXTRA_CFLAGS} -I${INCDIR} $<

fractionalOctaveSmoothtm.c : fractionalOctaveSmooth.tm
	${MPREP} $? -o $@

clean :
	@ ${RM} -rf *.o *tm.c
