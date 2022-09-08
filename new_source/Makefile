FC=gfortran
EXE=rs.x
ifeq ($(FC),ifort)
LFLAGS= -Nmkl -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
else
ifeq ($(if $(MKLROOT),1,0),1)
LFLAGS= -Nmkl -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
else
LFLAGS= -llapack -lblas
endif
endif
OPTFLAGS=-g -fbacktrace -fcheck=all -cpp -Ofast -ffree-line-length-0 -fpic 
OBJ=safe_alloc.o logger.o face.o precision.o math.o string.o control.o lattice.o potential.o self.o charge.o calculation.o element.o globals.o symbolic_atom.o hamiltonian.o recursion.o os.o green.o density_of_states.o bands.o mix.o namelist_generator.o timer.o
CPPFLAGS=-DDEBUG -DCOLOR
FLAGS=${OPTFLAGS} ${CPPFLAGS}

all: compile

librslmto.so: ${OBJ}
	${FC} -shared ${FLAGS} $^  -o $@ ${LFLAGS}

compile: main.o librslmto.so
	${FC} ${FLAGS} -o ${EXE} $< ${LFLAGS} -lrslmto -L$(CURDIR) -Wl,-rpath=$(CURDIR)

compile-no-lib: main.o ${OBJ}
	${FC} ${FLAGS} -o ${EXE} $^ ${LFLAGS}

run: ${EXE}
	./$<

%.o: %.f90
	${FC} ${FLAGS} -c $<

%.o: %.F90
	${FC} ${FLAGS} -c $<

# dependencies
safe_alloc.o: string.o
logger.o: string.o face.o
string.o: precision.o
math.o: precision.o
calculation.o: self.o lattice.o control.o charge.o symbolic_atom.o hamiltonian.o recursion.o string.o green.o density_of_states.o bands.o
charge.o: lattice.o precision.o math.o string.o self.o
self.o: math.o lattice.o control.o precision.o
lattice.o: string.o math.o control.o precision.o
# potential.o: lattice.o precision.o string.o
main.o: calculation.o precision.o os.o logger.o
element.o: globals.o precision.o
potential.o: globals.o precision.o
symbolic_atom.o: element.o potential.o lattice.o
hamiltonian.o: symbolic_atom.o element.o potential.o charge.o lattice.o precision.o math.o
recursion.o: hamiltonian.o lattice.o
green.o: self.o control.o lattice.o symbolic_atom.o recursion.o density_of_states.o precision.o math.o
density_of_states.o: control.o lattice.o self.o symbolic_atom.o recursion.o precision.o math.o
bands.o: control.o self.o green.o lattice.o symbolic_atom.o density_of_states.o precision.o math.o
mix.o: precision.o
namelist_generator.o: precision.o string.o
timer.o: precision.o string.o globals.o

clean:
	rm -rdf ${EXE} *.mod *.o *.so *.x
