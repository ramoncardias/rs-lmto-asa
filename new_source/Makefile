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
OPTFLAGS=-fallow-argument-mismatch -O2 -ffree-line-length-0 -g -fbacktrace -cpp #-g -fbacktrace -cpp -Ofast -ffree-line-length-0 -fpic -fcheck=all 
OBJ=safe_alloc.o logger.o face.o precision.o math.o string.o control.o lattice.o potential.o energy.o self.o charge.o calculation.o element.o globals.o symbolic_atom.o hamiltonian.o recursion.o os.o green.o density_of_states.o bands.o mix.o namelist_generator.o timer.o report.o xc.o
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
control.o: precision.o string.o namelist_generator.o logger.o
safe_alloc.o: string.o report.o
logger.o: string.o face.o report.o
string.o: precision.o
math.o: precision.o
calculation.o: self.o lattice.o control.o charge.o symbolic_atom.o hamiltonian.o recursion.o string.o green.o density_of_states.o bands.o mix.o
charge.o: lattice.o precision.o math.o string.o
self.o: math.o lattice.o control.o charge.o precision.o symbolic_atom.o logger.o xc.o recursion.o density_of_states.o green.o bands.o energy.o mix.o timer.o
energy.o: lattice.o
lattice.o: string.o math.o control.o precision.o symbolic_atom.o namelist_generator.o
# potential.o: lattice.o precision.o string.o
main.o: calculation.o precision.o os.o logger.o
element.o: globals.o precision.o logger.o string.o
potential.o: globals.o precision.o logger.o string.o
symbolic_atom.o: math.o element.o potential.o
hamiltonian.o: symbolic_atom.o element.o potential.o charge.o lattice.o precision.o math.o
recursion.o: hamiltonian.o lattice.o energy.o
green.o: energy.o control.o lattice.o symbolic_atom.o recursion.o density_of_states.o precision.o math.o
density_of_states.o: control.o lattice.o energy.o symbolic_atom.o recursion.o precision.o math.o
bands.o: control.o energy.o green.o lattice.o symbolic_atom.o density_of_states.o precision.o math.o
mix.o: precision.o string.o control.o lattice.o symbolic_atom.o
namelist_generator.o: precision.o string.o
timer.o: precision.o string.o globals.o
report.o: precision.o string.o
xc.o: control.o logger.o string.o precision.o math.o

globals.o: globals.f90
	#@(awk -i inplace -f .makefile_scripts/global_root_folder_filler.awk globals.f90)
	${FC} ${FLAGS} -c $<

clean:
	rm -rdf ${EXE} *.mod *.o *.so *.x
