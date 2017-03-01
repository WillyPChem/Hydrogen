LIB        = -L. -lm 
INCLUDE    = -I.
CFLAGS     = -O2
EXEC       = Radial.exe
CXX        = g++

${EXEC}: EXAM4.c Radial.o
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} EXAM4.c Radial.o -o ${EXEC}

blas.o: Radial.c Radial.h
	${CXX} ${LIB} -c Radial.c ${CFLAGS}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<
