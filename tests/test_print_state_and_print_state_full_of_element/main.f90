program Main
    use unit_test
    use element_mod

    type(test_suite_type) :: test_suite
    type(element) :: Fe
    
    ! (1) Create a test
    ! call test_suite_init('TEST PRINT STATE AND PRINT STATE FULL OF ELEMENT', test_suite)

    ! (2) Create as many cases as you need
    ! call test_case_create('CASE TITLE', test_suite)

    ! (3) Test your variables:

    Fe = element('Fe')
    call Fe%print_state()
    call Fe%print_state_full()
    ! Number, character and vector test examples
    ! - call assert_equal(var1, var2, __FILE__, __LINE__, test_suite)
    ! - call assert_great_than(var1, var2, __FILE__, __LINE__, test_suite)
    ! - call assert_approximate(var1, var2, eps=tolerance, __FILE__, __LINE__, test_suite)

    ! Logical test example
    ! - call assert_false(var, __FILE__, __LINE__, test_suite)
    ! - call assert_true(var, __FILE__, __LINE__, test_suite)

    ! More details open `src/assert_test.F90`

    ! (4) Finish your test with the two following lines
    ! call test_suite_report(test_suite)
    ! call test_suite_final(test_suite)

end program Main
