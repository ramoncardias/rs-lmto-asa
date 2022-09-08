program Main
    use unit_test
    use symbolic_atom_mod

    type(test_suite_type) :: test_suite
    type(symbolic_atom), dimension(3) :: array, array_read
    character(len=300) :: content_of_result_with_file, content_of_result_with_unit
    integer :: unit1,unit2
    
    ! (1) Create a test
    call test_suite_init('TEST SYMBOLIC ATOM PRINT STATE', test_suite)

    ! (2) Create as many cases as you need
    call test_case_create('RESULTS OF PRINT_STATE WITH UNIT AND FILE PARAMETERS', test_suite)

    ! (3) Test your variables:
    array = array_of_symbolic_atoms((/'Fe','Co','Cu'/),'./')

    open(newunit=unit1,file='result_with_unit.nml',action='write')
    call print_state_full(array,unit=unit1)
    close(unit1)
    
    call print_state_full(array,file='result_with_file.nml')

    open(newunit=unit1,file='result_with_file.nml',action='read')
    open(newunit=unit2,file='result_with_unit.nml',action='read')
    do
        read(unit1,*,end=10) content_of_result_with_file
        read(unit2,*,end=10) content_of_result_with_unit
        call assert_equal(content_of_result_with_file,content_of_result_with_unit,__FILE__, __LINE__, test_suite)
    enddo
    10 continue
    close(unit1)
    close(unit2)

    
    call test_case_create('RESULTS OF PRINT_STATE WITH SUFFIX PARAMETER', test_suite)
    
    call print_state(array,suffix='_out')
    
    array_read = array_of_symbolic_atoms((/'Fe_out','Co_out','Cu_out'/),'./')

    do i=1,size(array)
        call assert_equal(array(i)%element%symbol,array_read(i)%element%symbol,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%atomic_number,array_read(i)%element%atomic_number,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%core,array_read(i)%element%core,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%valence,array_read(i)%element%valence,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%f_core,array_read(i)%element%f_core,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%num_quant_s,array_read(i)%element%num_quant_s,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%num_quant_p,array_read(i)%element%num_quant_p,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%num_quant_d,array_read(i)%element%num_quant_d,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%q0,array_read(i)%element%q0,__FILE__, __LINE__, test_suite)
    enddo

    call test_case_create('RESULTS OF PRINT_STATE_FULL WITH SUFFIX PARAMETER', test_suite)
    
    call print_state_full(array,suffix='_out')
    
    array_read = array_of_symbolic_atoms((/'Fe_out','Co_out','Cu_out'/),'./')

    do i=1,size(array)
        call assert_equal(array(i)%element%symbol,array_read(i)%element%symbol,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%atomic_number,array_read(i)%element%atomic_number,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%core,array_read(i)%element%core,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%valence,array_read(i)%element%valence,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%f_core,array_read(i)%element%f_core,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%num_quant_s,array_read(i)%element%num_quant_s,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%num_quant_p,array_read(i)%element%num_quant_p,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%num_quant_d,array_read(i)%element%num_quant_d,__FILE__, __LINE__, test_suite)
        call assert_equal(array(i)%element%q0,array_read(i)%element%q0,__FILE__, __LINE__, test_suite)
    enddo

    ! (4) Finish your test with the two following lines
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
