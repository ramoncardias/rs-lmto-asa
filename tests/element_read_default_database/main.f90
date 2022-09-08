program Main
    use element_mod
    use unit_test

    type(element) :: el
    type(element), dimension(4) :: ar
    integer :: unit,iostatus
    character(len=6) :: symbol
    integer :: atomic_number, core, valence, f_core, num_quant_s, num_quant_p, num_quant_d
    type(test_suite_type) :: test_suite
    
    ! (1) Create a test
    call test_suite_init('Element Read Default Database', test_suite)

    ! (2) Create as many cases as you need
    call test_case_create('Testing if all namelists are correctly readed',test_suite)

    ! (3) Test your variables:

    open(newunit=unit,file='elements-per-table')
    ! jump one line
    read(unit,'(A)',iostat=iostatus) 
    do while (.TRUE.)
        read(unit,'(A,I6,I15,I8,I16,I12,I21,I15)',iostat=iostatus) symbol,atomic_number,core,valence,f_core,num_quant_s,num_quant_p,num_quant_d
        if(iostatus /= 0) exit
        el = element(symbol)
        call assert_equal(el%symbol, symbol, __FILE__, __LINE__, test_suite)
        call assert_equal(el%atomic_number, atomic_number, __FILE__, __LINE__, test_suite)
        call assert_equal(el%core, core, __FILE__, __LINE__, test_suite)
        call assert_equal(el%valence, valence, __FILE__, __LINE__, test_suite)
        call assert_equal(el%f_core, f_core, __FILE__, __LINE__, test_suite)
        call assert_equal(el%num_quant_s, num_quant_s, __FILE__, __LINE__, test_suite)
        call assert_equal(el%num_quant_p, num_quant_p, __FILE__, __LINE__, test_suite)
        call assert_equal(el%num_quant_d, num_quant_d, __FILE__, __LINE__, test_suite)
    enddo
    close(unit)

    call test_case_create('Loading elements with function "array_of_elements"',test_suite)


    ar = array_of_elements((/'Fe','O ','C ','H '/))
    
    call assert_equal(ar(1)%symbol, 'Fe', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%atomic_number,  26, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%valence, 8, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call assert_equal(ar(2)%symbol, 'O', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%atomic_number, 8, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%core, 2, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%valence, 6, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%num_quant_s, 2, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%num_quant_p, 2, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call assert_equal(ar(3)%symbol, 'C', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%atomic_number, 6 , __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%core, 2, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%valence, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%num_quant_s, 2, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%num_quant_p, 2, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%num_quant_d, 3, __FILE__, __LINE__, test_suite)
    
    call assert_equal(ar(4)%symbol, 'H', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%atomic_number, 1 , __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%valence, 1, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%num_quant_s, 1, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%num_quant_p, 2, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call test_case_create('Loading elements outside default database', test_suite)

    el = element('NewElem','./')

    call assert_equal(el%symbol, 'NewElem', __FILE__, __LINE__, test_suite)
    call assert_equal(el%atomic_number, 1 , __FILE__, __LINE__, test_suite)
    call assert_equal(el%core, 2, __FILE__, __LINE__, test_suite)
    call assert_equal(el%valence, 3, __FILE__, __LINE__, test_suite)
    call assert_equal(el%f_core, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(el%num_quant_s, 5, __FILE__, __LINE__, test_suite)
    call assert_equal(el%num_quant_p, 6, __FILE__, __LINE__, test_suite)
    call assert_equal(el%num_quant_d, 7, __FILE__, __LINE__, test_suite)

    ! (4) Finish your test with the two following lines
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
