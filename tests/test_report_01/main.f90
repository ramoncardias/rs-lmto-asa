    program Main
        use precision_mod, only: rp
    use string_mod, only: sl
    use report_mod
    use unit_test

    type(test_suite_type) :: test_suite
    type(report) :: rep
    type(report), pointer :: child
    real(rp), dimension(:), allocatable :: values
    character(len=sl), dimension(:), allocatable :: labels
    logical :: found
    
    ! (1) Create a test
    call test_suite_init('TEST REPORT 01', test_suite)

    ! (2) Create as many cases as you need

    call test_case_create('REP%ADD_VALUE', test_suite)

    rep = report('total')
    call rep%add_value('total',1.2)
    call rep%add_value('total.calc_01',3.4)
    call rep%add_value('total.calc_01',5.6)
    call rep%add_value('total.calc_01',7.8)

    call rep%add_value('total.calc_01.subcalc_01',9.1)
    call rep%add_value('total.calc_01.subcalc_01',11.12)
    call rep%add_value('total.calc_01.subcalc_01',13.14)

    call rep%add_value('total.calc_01.subcalc_01.subsub_01',1)
    call rep%add_value('total.calc_01.subcalc_01.subsub_01',2)
    call rep%add_value('total.calc_01.subcalc_01.subsub_01',3)

    call rep%add_value('total.calc_01.subcalc_01.subsub_01.subsubsub_01',4)
    call rep%add_value('total.calc_01.subcalc_01.subsub_01.subsubsub_01',5)
    call rep%add_value('total.calc_01.subcalc_01.subsub_01.subsubsub_01',6)

    call rep%add_value('total.calc_02',11.2)
    call rep%add_value('total.calc_02',23.4)
    call rep%add_value('total.calc_02',35.6)
    
    call assert_true(.True.,__FILE__, __LINE__, test_suite)

    call test_case_create('REP%GET_ABELS', test_suite)
    call rep%get_labels(labels)
    call assert_equal(labels(1),'total', __FILE__, __LINE__, test_suite)
    call assert_equal(labels(2),'total.calc_01', __FILE__, __LINE__, test_suite)
    call assert_equal(labels(3),'total.calc_01.subcalc_01', __FILE__, __LINE__, test_suite)
    call assert_equal(labels(4),'total.calc_01.subcalc_01.subsub_01', __FILE__, __LINE__, test_suite)
    call assert_equal(labels(5),'total.calc_01.subcalc_01.subsub_01.subsubsub_01', __FILE__, __LINE__, test_suite)
    call assert_equal(labels(6),'total.calc_02', __FILE__, __LINE__, test_suite)

    call test_case_create('REP%GET_NUM_CHILDREN', test_suite)

    call assert_equal(rep%get_num_children(),2,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total'),2,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_01'),1,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_01.subcalc_01'),1,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_01.subcalc_01.subsub_01'),1,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_01.subcalc_01.subsub_01.subsubsub_01'),0,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_02'),0,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_num_children('total.calc_02.sub_01'),0,__FILE__, __LINE__, test_suite)

    call test_case_create('REP%GET_SIZE_VALUES', test_suite)

    call assert_equal(rep%get_size_values(),1,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_01'),3,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_01.subcalc_01'),3,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_01.subcalc_01.subsub_01'),3,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_01.subcalc_01.subsub_01.subsubsub_01'),3,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_02'),3,__FILE__, __LINE__, test_suite)
    call assert_equal(rep%get_size_values('total.calc_02.sub_01'),0,__FILE__, __LINE__, test_suite)

    call test_case_create('REP%GET_VALUES', test_suite)
    call rep%get_values(values)
    call assert_approximate(real(values(1),4),real(1.2,4),__FILE__, __LINE__,0.05, test_suite)
    call rep%get_values('total',values)
    call assert_approximate(real(values(1),4),real(1.2,4),__FILE__, __LINE__,0.05, test_suite)
    
    call rep%get_child_values(values)

    call assert_approximate(real(values(1),4),real(3.4,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(2),4),real(5.6,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(3),4),real(7.8,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(4),4),real(9.1,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(5),4),real(11.12,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(6),4),real(13.14,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(7),4),real(1,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(8),4),real(2,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(9),4),real(3,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(10),4),real(4,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(11),4),real(5,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(12),4),real(6,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(13),4),real(11.2,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(14),4),real(23.4,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(15),4),real(35.6,4),__FILE__, __LINE__,0.05, test_suite)

    call test_case_create('REP%GET_VALUES APPEND', test_suite)

    call rep%get_values(values)
    call assert_equal(size(values),1,__FILE__, __LINE__,test_suite)
    call assert_approximate(real(values(1),4),real(1.2,4),__FILE__, __LINE__,0.05, test_suite)
    
    call rep%get_values(values,append=.true.)
    call assert_equal(size(values),2,__FILE__, __LINE__,test_suite)
    call assert_approximate(real(values(1),4),real(1.2,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(2),4),real(1.2,4),__FILE__, __LINE__,0.05, test_suite)

    call rep%get_values('total.calc_01.subcalc_01.subsub_01.subsubsub_01',values,append=.true.)
    call assert_equal(size(values),5,__FILE__, __LINE__,test_suite)
    call assert_approximate(real(values(1),4),real(1.2,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(2),4),real(1.2,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(3),4),real(4,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(4),4),real(5,4),__FILE__, __LINE__,0.05, test_suite)
    call assert_approximate(real(values(5),4),real(6,4),__FILE__, __LINE__,0.05, test_suite)

    call rep%get_child_values('total',values)
    call assert_approximate(real(values(1),4),real(3.4,4),__FILE__, __LINE__,0.05, test_suite)

    call test_case_create('REP%GET_PARENT', test_suite)

    call assert_false(rep%get_parent('total',child),__FILE__, __LINE__, test_suite)
    
    found = rep%get_parent('total.calc_01',child)
    call assert_true(found,__FILE__, __LINE__, test_suite)
    if (found) call assert_equal(child%label,'total',__FILE__, __LINE__, test_suite)
    
    found = rep%get_parent('total.calc_01.subcalc_01.subsub_01.subsubsub_01',child)
    call assert_true(found,__FILE__, __LINE__, test_suite)
    if (found) call assert_equal(child%label,'subsub_01',__FILE__, __LINE__, test_suite)

    call test_case_create('REPORT EXAMPLE', test_suite)
    call rep%print_report('REPORT EXAMPLE 01')
    call assert_true(.True.,__FILE__, __LINE__, test_suite)
    
    ! (4) Finish your test with the two following lines
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
