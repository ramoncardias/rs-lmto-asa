program Main
    use control_mod
    use unit_test

    type(test_suite_type) :: test_suite
    type(control) :: control_obj
    
    call test_suite_init('CONTROL MODULE TEST', test_suite)

    call test_case_create('Restore to Default')

    call control_obj%restore_to_default()

    call assert_equal(control_obj%llsp,       21,    suite=test_suite)
    call assert_equal(control_obj%lld,        21,    suite=test_suite)
    call assert_equal(control_obj%idos,       0,     suite=test_suite)
    call assert_equal(control_obj%mext,       0,     suite=test_suite)
    call assert_equal(control_obj%txc,        1,     suite=test_suite)
    call assert_equal(control_obj%partype,    1,     suite=test_suite)
    call assert_equal(control_obj%terminator, 5,     suite=test_suite)
    call assert_equal(control_obj%conca,      0.0d0, suite=test_suite)
    call assert_equal(control_obj%concb,      0.0d0, suite=test_suite)
    call assert_equal(control_obj%ruban,      0.0d0, suite=test_suite)
    call assert_equal(control_obj%njij,       0,     suite=test_suite)

    call assert_false(control_obj%lrot,     suite=test_suite)
    call assert_false(control_obj%incorb,   suite=test_suite)
    call assert_false(control_obj%svac,     suite=test_suite)
    call assert_false(control_obj%blockrec, suite=test_suite)
    call assert_false(control_obj%do_asd,   suite=test_suite)
    call assert_false(control_obj%do_cochg, suite=test_suite)
    call assert_false(control_obj%asd_jij,  suite=test_suite)
    call assert_false(control_obj%do_comom, suite=test_suite)

    
    call test_case_create('Values From Constructor')

    control_obj = control('input.nml')

    call assert_equal(control_obj%llsp,             456,           __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%lld,              789,           __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%idos,             123,           __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%mext,             456,           __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%txc,              123,           __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%partype,          123,           __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%terminator,       1,             __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%conca,            1.0d0,         __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%concb,            2.0d0,         __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%ruban,            3.0d0,         __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%njij,             1,             __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%calctype,         "S",           __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%nsp,              111,           __FILE__, __LINE__, suite=test_suite)
    call assert_equal(control_obj%nlim,             222,           __FILE__, __LINE__, suite=test_suite)
    
    call assert_true(control_obj%lrot,     __FILE__, __LINE__, suite=test_suite)
    call assert_true(control_obj%incorb,   __FILE__, __LINE__, suite=test_suite)
    call assert_true(control_obj%svac,     __FILE__, __LINE__, suite=test_suite)
    call assert_true(control_obj%blockrec, __FILE__, __LINE__, suite=test_suite)
    call assert_true(control_obj%do_asd,   __FILE__, __LINE__, suite=test_suite)
    call assert_true(control_obj%do_cochg, __FILE__, __LINE__, suite=test_suite)
    call assert_true(control_obj%asd_jij,  __FILE__, __LINE__, suite=test_suite)
    call assert_true(control_obj%do_comom, __FILE__, __LINE__, suite=test_suite)


    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
