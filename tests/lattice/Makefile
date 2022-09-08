FC=gfortran
EXE=test
LFLAGS= -lblas -llapack -lunittest -L ../src
OPTFLAGS=-cpp -ffree-line-length-512

all: compile run

compile: control.o lattice.o main.o
	${FC} ${OPTFLAGS} $^ -o ${EXE} ${LFLAGS}

run: ${EXE}
	./$<

%.o: %.f90
	${FC} ${OPTFLAGS} -c $<


# dependencies
lattice.o: control.o
main.o: lattice.o

clean:
	rm -rdf ${EXE} *.mod *.o