tests=$(subst /,, $(subst src/,,$(wildcard */)))
compile-tests=$(patsubst %,compile-%,${tests})
run-tests=$(patsubst %,run-%,${tests})

FC=gfortran
EXE=test
LFLAGS= -lblas -llapack -lunittest -L ../src -lrslmto -L../${SRC_DIR} -Wl,-rpath=../${SRC_DIR} -DDEBUG
OPTFLAGS=-cpp -ffree-line-length-512
SRC=calculation charge control lattice logger math precision self string face element globals potential symbolic_atom
SRC_DIR=../new_source

help:
	@echo "To use this make file execute the options below:"
	@echo
	@echo " 1. To create a new test with the name 'test01':"
	@echo "  make create test01"
	@echo
	@echo " 2. To compile and run a test with the name 'test01':"
	@echo "  make test01"
	@echo
	@echo " 3. To remove a test with the name 'test01':"
	@echo "  rm -rdf test01"
	@echo
	@echo " 4. To compile all tests:"
	@echo "  make all"
	@echo
	@echo " 5. To clean all tests:"
	@echo "  make clean"
	@echo
	@echo " 6. To show this message:"
	@echo "  make help"
	@echo

*:
	@(if [ -f .create ] ; then \
		echo "Test '$@' already exists! Remove it if you want to create a new one!";\
		rm .create; \
		exit; \
	else \
		# echo 11; \
		${MAKE} -C ${SRC_DIR} librslmto.so; \
		cd $@; \
		ln -sf ../${SRC_DIR}/*.mod ../src/*.mod .; \
		${FC} ${OPTFLAGS} -o ${EXE} main.f90 ${LFLAGS} && ./${EXE}; \
		for link in `ls *.mod`; do unlink $$link; done; \
	fi)

.PHONY: clean show ${tests} help all

create:
	@touch .create

clean:
	@(for folder in ${tests};do rm -f $${folder}/test;done)

all:
	@(for folder in ${tests};do ${MAKE} $${folder};done)

.DEFAULT:
	@(if [ -e .create ] && [ ! -e $@ ]; then \
		echo "Creating test $@"; \
		rm .create; \
		mkdir -p $@; \
		test_title=$(shell echo '$@' | tr '[:lower:]' '[:upper:]' | sed "s/\_/' '/g"); \
		sed "s/TEST TITLE/$$test_title/1" src/example.f90 > $@/main.f90; \
	else \
		echo "Test '$@' does not exist!";\
	fi)

show:
	@echo ${SRC}