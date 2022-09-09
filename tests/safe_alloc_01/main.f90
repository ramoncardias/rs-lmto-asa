program Main
    use unit_test
    ! use globals_mod
    use string_mod, only: sl
    use safe_alloc_mod, only: safe_alloc, g_safe_alloc, free

    type(test_suite_type) :: test_suite
    integer, dimension(:), allocatable :: v
    integer :: i, j, k, vec_size

    character(len=sl), dimension(5) :: mod_labels
    character(len=sl), dimension(3) :: calc_labels
    character(len=sl), dimension(4) :: array_labels
    mod_labels(1) = 'selfcon'
    mod_labels(2) = 'mix'
    mod_labels(3) = 'hamiltonian'
    mod_labels(4) = 'green'
    mod_labels(5) = 'charge'

    calc_labels(1) = 'calc_01'
    calc_labels(2) = 'lapack'
    calc_labels(3) = 'matmul'

    array_labels(1) = 'array'
    array_labels(2) = 'vector'
    array_labels(3) = 'electrons'
    array_labels(4) = 'atoms'
    
    g_safe_alloc = safe_alloc()
    
    ! (1) Create a test
    call test_suite_init('SAFE ALLOC 01', test_suite)

    call test_case_create('ALLOATE/DEALLOCATE', test_suite)

    do i=1,size(mod_labels)
        do j=1,size(calc_labels)
            do k=1,size(array_labels)
                vec_size = 10*i*j*k
                call g_safe_alloc%allocate('calculation.'//trim(mod_labels(i))//'.'//trim(calc_labels(j))//'.'//trim(array_labels(k)),v,vec_size)
                call assert_equal(size(v), vec_size, __FILE__, __LINE__, test_suite)
                call g_safe_alloc%deallocate('calculation.'//trim(mod_labels(i))//'.'//trim(calc_labels(j))//'.'//trim(array_labels(k)),v)
            enddo
        enddo
    enddo
    
    call g_safe_alloc%print_report()

    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
