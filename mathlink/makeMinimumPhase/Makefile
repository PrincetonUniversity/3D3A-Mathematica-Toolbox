# This makefile can be used to build the makeMinimumPhase
# program.  To build it, use the command 'make makeMinimumPhase'.

VERSION=10.0
MLINKDIR = /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64
SYS = MacOSX-x86-64
CADDSDIR = ${MLINKDIR}/CompilerAdditions
EXTRA_CFLAGS=-stdlib=libstdc++ -stdlib=libc++
EXTRA_LIBS = -lstdc++ -lfftw3 -lm -framework Foundation

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

MPREP = ${CADDSDIR}/mprep
RM = rm

all : makeMinimumPhase clean

makeMinimumPhase : makeMinimumPhasetm.o makeMinimumPhase.o
	${CXX} ${EXTRA_CFLAGS} -I${INCDIR} makeMinimumPhasetm.o makeMinimumPhase.o -L${LIBDIR} -lMLi4 -lc++ ${EXTRA_LIBS} -o $@

.c.o :
	${CC} -c ${EXTRA_CFLAGS} -I${INCDIR} $<

makeMinimumPhasetm.c : makeMinimumPhase.tm
	${MPREP} $? -o $@

clean :
	@ ${RM} -rf *.o *tm.c
