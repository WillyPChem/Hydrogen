LIB        = -L. -lm 
INCLUDE    = -I.
CFLAGS     = -O2
EXEC       = rad.exe
CXX        = g++

${EXEC}: Exam4.c Radial.o
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} Exam4.c Radial.o -o ${EXEC}

Radial.o: Radial.c Radial.h
	${CXX} ${LIB} -c Radial.c ${CFLAGS}

clean:
	rm -f *.o
