FC=gfortran
EXE=test
LFLAGS= -lblas -llapack -lunittest -L ../src
OPTFLAGS=-cpp -ffree-line-length-512

all: compile run

compile: control.o lattice.o self.o calculation.o main.o
	${FC} ${OPTFLAGS} $^ -o ${EXE} ${LFLAGS}

run: ${EXE}
	./$<

%.o: %.f90
	${FC} ${OPTFLAGS} -c $<


# dependencies
calculation.o: self.o lattice.o control.o
self.o: control.o
lattice.o: control.o
main.o: calculation.o

clean:
	rm -rdf ${EXE} *.mod *.o