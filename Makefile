LIB        = -L. -lm 
INCLUDE    = -I.
CFLAGS     = -O2
EXEC       = Hydrogen.exe
CXX        = g++

${EXEC}: Hydrogen_Simulator.c Radial.o
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} Hydrogen_Simulator.c Radial.o -o ${EXEC}

Radial.o: Radial.c Radial.h
	${CXX} ${LIB} -c Radial.c ${CFLAGS}

clean:
	rm -f *.o
