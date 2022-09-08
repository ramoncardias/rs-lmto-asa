program Main
    use unit_test
    ! use globals_mod
    use safe_alloc_mod, only: safe_alloc, g_safe_alloc, free

    type(test_suite_type) :: test_suite
    integer, dimension(:), allocatable :: v
    
    g_safe_alloc = safe_alloc()
    
    ! (1) Create a test
    call test_suite_init('SAFE ALLOC 01', test_suite)

    call test_case_create('CASE TITLE', test_suite)

    call g_safe_alloc%alloc('v',v,10)
    call assert_equal(size(v), 10, __FILE__, __LINE__, test_suite)

    call g_safe_alloc%dealloc('v',v)

    call g_safe_alloc%alloc('v',v,100)
    call assert_equal(size(v), 100, __FILE__, __LINE__, test_suite)

    call g_safe_alloc%report()

    call free(g_safe_alloc)

    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
