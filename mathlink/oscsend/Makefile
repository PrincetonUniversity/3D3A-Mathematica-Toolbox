# This makefile can be used to build the oscsend
# program.  To build it, use the command 'make oscsend'.

VERSION=10.0
MLINKDIR = /Applications/Mathematica.app/Contents/SystemFiles/Links/MathLink/DeveloperKit/MacOSX-x86-64
SYS = MacOSX-x86-64
CADDSDIR = ${MLINKDIR}/CompilerAdditions
EXTRA_CFLAGS=-stdlib=libstdc++ -stdlib=libc++
EXTRA_CFILES=htmsocket.c OSC-client.c OSC-timetag.c

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

MPREP = ${CADDSDIR}/mprep
MCC = ${CADDSDIR}/mcc
RM = rm

all : oscsend clean

oscsend :
	${MCC} ${EXTRA_CFLAGS} -I${INCDIR} oscsend.tm oscsend.c ${EXTRA_CFILES} -o $@

clean :
	@ ${RM} -rf *.o *tm.c
