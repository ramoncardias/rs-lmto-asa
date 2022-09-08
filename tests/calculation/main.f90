program Main
    use calculation_mod
    use unit_test
    implicit none

    type(test_suite_type) :: test_suite
    type(calculation) :: calculation_obj
    

    call test_suite_init('CALCULATION INITALIZATION TEST', test_suite)

    call test_case_create('Restore to Default')

    call calculation_obj%restore_to_default()

    call assert_equal(calculation_obj%pre_processing, 'none',file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(calculation_obj%processing, 'none',file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(calculation_obj%post_processing, 'none',file_name=__FILE__, line_number=__LINE__, suite=test_suite)

    call test_case_create('Build from Constructor')

    calculation_obj = calculation('input.nml')

    call assert_equal(calculation_obj%pre_processing, 'bravais',file_name=__FILE__, line_number=__LINE__, suite=test_suite)

  
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)
end program Main
