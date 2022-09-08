# assert_mod.F90  assert_test.F90  test_case_mod.F90  test_common_mod.F90  test_suite_mod.F90  unit_test.F90


FC=gfortran
FLAGS=-cpp -ffree-line-length-512

all: unit_test test

unit_test: unit_test.F90 unit_test.o assert_mod.o test_case_mod.o  test_common_mod.o  test_suite_mod.o
	$(FC) $(FLAGS) -c $<
	ar rU libunittest.a $^
	ranlib libunittest.a

test: assert_test.F90 assert_test.o
	$(FC) -o $@ $(FLAGS) $< -lunittest -L .

assert_mod.o: test_common_mod.o test_case_mod.o
assert_test.o: unit_test.o
test_case_mod.o: test_common_mod.o 
test_suite_mod.o: test_common_mod.o
unit_test.o: assert_mod.o test_case_mod.o test_suite_mod.o

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -rdf *.mod *.o *.a test