program Main
    use lattice_mod
    use control_mod
    use self_mod
    use unit_test

    type(test_suite_type) :: test_suite
    type(lattice) :: lattice_obj
    type(control), target :: control_obj
    type(self) :: self_obj

    type(control), pointer :: control_pointer

    control_obj = control('input.nml')
    lattice_obj = lattice('input.nml',control_obj)
    call lattice_obj%build_data()
    call lattice_obj%bravais()
    call lattice_obj%newclu()
    call lattice_obj%structb()
    self_obj = self('input.nml',lattice_obj)

    control_pointer => control_obj

    ! (1) Create a test
    call test_suite_init('Lattice Control dependency', test_suite)

    ! (2) Create as many cases as you need
    call test_case_create('CASE TITLE')

    ! (3) Test your variables:
    
    select type(self_obj%lattice%control)
        type is (self_obj)
            write(*,*) 'type is self'
        type is (control_obj)
            write(*,*) 'type is control'
        type is (lattice_obj)
            write(*,*) 'type is lattice'
    end select
    ! Number, character and vector test examples
    ! - call assert_equal(var1, var2, __FILE__, __LINE__,test_suite)
    ! - call assert_great_than(var1, var2, __FILE__, __LINE__,test_suite)
    ! - call assert_approximate(var1, var2, eps=tolerance, __FILE__, __LINE__,test_suite)

    ! Logical test example
    ! - call assert_false(var, __FILE__, __LINE__,test_suite)
    ! - call assert_true(var, __FILE__, __LINE__,test_suite)

    ! More details open `src/assert_test.F90`

    ! (4) Finish your test with the two following lines
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
