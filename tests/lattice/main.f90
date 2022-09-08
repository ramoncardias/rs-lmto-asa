program Main
    use lattice_mod
    use control_mod
    use unit_test

    implicit none

    type(test_suite_type) :: test_suite
    type(lattice) :: lattice_obj
    type(control) :: control_obj
    call test_suite_init('LATTICE TEST', test_suite)

    call test_case_create('Restore to Default',test_suite)

    call lattice_obj%restore_to_default()

    call assert_equal(lattice_obj%ndim, 9900000, file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%npe, 49, file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(size(lattice_obj%izp), lattice_obj%ndim,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(size(lattice_obj%no), lattice_obj%ndim,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(size(lattice_obj%crd), 3*lattice_obj%ndim,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(size(lattice_obj%inclu), 3*lattice_obj%nclu,file_name=__FILE__, line_number=__LINE__, suite=test_suite)

    call test_case_create('Values From Constructor',test_suite)

    lattice_obj = lattice('input.nml')
    
    ! General intialization
    call assert_equal(lattice_obj%iu(:), [1],file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%ib(:), [2],file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%irec(:), [3],file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    
    ! Bulk initialization
    call assert_equal(lattice_obj%ndim, 9,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%npe, 123,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%rc, real(100,8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%crystal_sym, 'bcc',file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%a(1,:), real([1,2,3],8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%a(2,:), real([4,5,6],8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%a(3,:), real([7,8,9],8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    
    call assert_equal(size(lattice_obj%izp), lattice_obj%ndim,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(size(lattice_obj%no), lattice_obj%ndim,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(size(lattice_obj%crd), 3 * lattice_obj%ndim,file_name=__FILE__, line_number=__LINE__, suite=test_suite)

    call assert_equal(lattice_obj%izp(:), [11,12,13,14,15,16,17,18,19],file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%no(:), [21,22,23,24,25,26,27,28,29],file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_approximate(lattice_obj%crd(1,:), real([31,32,33,34,35,36,37,38,39],8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_approximate(lattice_obj%crd(2,:), real([41,42,43,44,45,46,47,48,49],8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_approximate(lattice_obj%crd(3,:), real([51,52,53,54,55,56,57,58,59],8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    
    ! ! Impurity initialization
    call assert_equal(lattice_obj%nclu, 2,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(size(lattice_obj%inclu), lattice_obj%nclu*3,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%inclu(1,:), real([1,2,3],8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%inclu(2,:), real([4,5,6],8),file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    
    ! ! Surface initialization
    call assert_equal(lattice_obj%surftype, "abc",file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(lattice_obj%nlay, 1,file_name=__FILE__, line_number=__LINE__, suite=test_suite)

    call test_suite_report(test_suite)
    call test_suite_final(test_suite)
end program Main