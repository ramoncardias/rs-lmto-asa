program Main
    use lattice_mod
    use unit_test
    implicit none

    type(lattice) :: obj
    type(test_suite_type) :: test_suite

    call test_suite_init('TEST: "calculation_nbas"', test_suite)
    call test_case_create('CASE: "reduced_nbas" and "reduced_acr"',test_suite)

    obj%nbas = 10
    allocate(obj%iz(obj%nbas))
    
    obj%iz(:) = [10,11,2,2,2,3,3,3,4,4]
    call obj%calculate_nbas()
    call assert_equal(obj%reduced_nbas, 5,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(obj%reduced_acr(:), [10,11,2,3,4],file_name=__FILE__, line_number=__LINE__, suite=test_suite)

    obj%iz(:) = [2,11,2,2,2,3,3,3,4,4]
    call obj%calculate_nbas()
    call assert_equal(obj%reduced_nbas, 5,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(obj%reduced_acr(:), [2,11,2,3,4],file_name=__FILE__, line_number=__LINE__, suite=test_suite)


    obj%iz(:) = [2,3,2,2,2,3,3,3,4,4]
    call obj%calculate_nbas()
    call assert_equal(obj%reduced_nbas, 5,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(obj%reduced_acr(:), [2,3,2,3,4],file_name=__FILE__, line_number=__LINE__, suite=test_suite)

    obj%iz(:) = [3,3,2,2,2,3,3,3,4,4]
    call obj%calculate_nbas()
    call assert_equal(obj%reduced_nbas, 4,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(obj%reduced_acr(:), [3,2,3,4],file_name=__FILE__, line_number=__LINE__, suite=test_suite)

    
    obj%nbas = 5
    deallocate(obj%iz)
    allocate(obj%iz(obj%nbas))
    obj%iz(:) = [1,2,3,1,1]
    call obj%calculate_nbas()
    call assert_equal(obj%reduced_nbas, 4,file_name=__FILE__, line_number=__LINE__, suite=test_suite)
    call assert_equal(obj%reduced_acr(:), [1,2,3,1],file_name=__FILE__, line_number=__LINE__, suite=test_suite)


    call test_suite_report(test_suite)
    call test_suite_final(test_suite)
end program Main
