program Main
    use precision_mod, only: rp
    use string_mod, only: sl
    use report_mod
    use unit_test

    type(test_suite_type) :: test_suite
    type(report) :: rep
    real(rp), dimension(:), allocatable :: values
    character(len=sl), dimension(:), allocatable :: labels
    
    ! (1) Create a test
    call test_suite_init('TEST REPORT 02', test_suite)

    
    rep = report('total')

    call rep%add_value('total.calc_01.subcalc_01',9.1)
    call rep%add_value('total.calc_01.subcalc_01',11.12)
    call rep%add_value('total.calc_01.subcalc_01',13.14)
    
    call rep%add_value('total.calc_01.subcalc_02',3)
    call rep%add_value('total.calc_01.subcalc_02',5)
    call rep%add_value('total.calc_01.subcalc_02',4)
    call rep%add_value('total.calc_01.subcalc_02',5.5)

    
    call rep%add_value('total.calc_02.subcalc_01',1.2)
    call rep%add_value('total.calc_02.subcalc_01',3.4)
    call rep%add_value('total.calc_02.subcalc_01',5.6)
    
    call test_case_create('GET_LABELS', test_suite)
    
    call rep%get_labels(labels)
    
    call assert_equal(labels(1),'total',__FILE__, __LINE__, test_suite)
    call assert_equal(labels(2),'total.calc_01',__FILE__, __LINE__, test_suite)
    call assert_equal(labels(3),'total.calc_01.subcalc_01',__FILE__, __LINE__, test_suite)
    call assert_equal(labels(4),'total.calc_01.subcalc_02',__FILE__, __LINE__, test_suite)
    call assert_equal(labels(5),'total.calc_02',__FILE__, __LINE__, test_suite)
    call assert_equal(labels(6),'total.calc_02.subcalc_01',__FILE__, __LINE__, test_suite)
    
    call test_case_create('REP%GET_NUM_CHILDREN', test_suite)

    call assert_equal(rep%get_num_children(),2,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total'),2,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_01'),2,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_01.subcalc_01'),0,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_01.subcalc_02'),0,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_02'),1,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_02.subcalc_01'),0,__FILE__, __LINE__, test_suite)

    call test_case_create('REP%GET_SIZE_VALUES', test_suite)

    call assert_equal(rep%get_size_values(),0,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total'),0,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_01'),0,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_01.subcalc_01'),3,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_01.subcalc_02'),4,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_02'),0,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_02.subcalc_01'),3,__FILE__, __LINE__, test_suite)

    call test_case_create('VALUES', test_suite)

    call rep%get_child_values(values)

    call assert_equal(real(values(1),4),real(9.1,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(2),4),real(11.12,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(3),4),real(13.14,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(4),4),real(3,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(5),4),real(5,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(6),4),real(4,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(7),4),real(5.5,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(8),4),real(1.2,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(9),4),real(3.4,4),__FILE__,__LINE__,test_suite)
    call assert_equal(real(values(10),4),real(5.6,4),__FILE__,__LINE__,test_suite)

    call test_case_create('REPORT EXAMPLE', test_suite)

    call rep%print_report('REPORT EXAMPLE 02')

    call assert_true(.True.,__FILE__, __LINE__, test_suite)




    ! (4) Finish your test with the two following lines
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
