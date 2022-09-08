FC=gfortran
EXE=test
LFLAGS= -lblas -llapack -lunittest -L ../src
CFLAGS=-cpp -ffree-line-length-512

all: compile run

debug: CFLAGS+=-DDEBUG -g
debug: all

compile: main.f90 logger.o face.o
	${FC} ${CFLAGS} $^ -o ${EXE} ${LFLAGS}

run: ${EXE}
	./$<

%.o: %.f90
	${FC} ${CFLAGS} -c $<


# dependencies
logger.o: face.o
face.o: face.F90
	${FC} ${CFLAGS} -c $<

clean:
	rm -rdf ${EXE} *.mod *.o