program Main
    use unit_test
    use namelist_generator_mod

    type(test_suite_type) :: test_suite
    type(namelist_generator) :: namelist_gen
    real(8), dimension(2) :: array
    
    ! (1) Create a test
    call test_suite_init('TEST NAMELIST GENERATOR 01', test_suite)

    ! (2) Create as many cases as you need
    call test_case_create('CREATING NEW NAMELIST', test_suite)

    ! (3) Test your variables:

    namelist_gen = namelist_generator('my_new_namelist')
    call namelist_gen%add_variable('a',2)
    call namelist_gen%add_variable('b',4)
    call namelist_gen%add_variable('c',3.14152122)
    call namelist_gen%add_variable('d','testando...')
    call namelist_gen%add_variable('e(1)',1)
    call namelist_gen%add_variable('e(2)',2)
    call namelist_gen%add_variable('e(3)',3)
    call namelist_gen%add_variable('e(4)',4)
    call namelist_gen%add_variable('e(5)',5)
    call namelist_gen%add_variable('f',.True.)
    call namelist_gen%add_variable('g',.False.)
    array = (/1.2,2.0/)
    call namelist_gen%add_variable('h',array)
    call namelist_gen%generate_namelist('my_new_namelist_test')


    ! Number, character and vector test examples
    ! - call assert_equal(var1, var2, __FILE__, __LINE__, test_suite)
    ! - call assert_great_than(var1, var2, __FILE__, __LINE__, test_suite)
    ! - call assert_approximate(var1, var2, eps=tolerance, __FILE__, __LINE__, test_suite)

    ! Logical test example
    ! - call assert_false(var, __FILE__, __LINE__, test_suite)
    ! - call assert_true(var, __FILE__, __LINE__, test_suite)

    ! More details open `src/assert_test.F90`

    ! (4) Finish your test with the two following lines
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
