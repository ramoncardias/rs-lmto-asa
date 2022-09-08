program Main
    use lattice_mod
    use unit_test
    use symbolic_atom_mod
    type(symbolic_atom), dimension(:), allocatable :: ar

    type(test_suite_type) :: test_suite
    type(lattice) :: lattice_obj
    
    ! (1) Create a test
    call test_suite_init('Symbolic Atom Test Reading From File', test_suite)

    ! (2) Create as many cases as you need
    call test_case_create('CASE TITLE', test_suite)

    ! (3) Test your variables:

    lattice_obj%ntype = 4
    ar = array_of_symbolic_atoms('input.nml',lattice_obj)

    call assert_equal(ar(1)%element%symbol, 'Fe', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%atomic_number, 26, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%valence, 8, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(1)%potential%center_band_s_up),  11.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_p_up), -21.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_d_up),  31.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_s_dw), -41.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_p_dw),  51.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_d_dw), -61.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_s_up),   71.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_p_up),  -81.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_d_up),   91.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_s_dw), -101.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_p_dw),  111.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_d_dw), -121.5, __FILE__, __LINE__, test_suite)

    
    call assert_equal(ar(2)%element%symbol, 'Fe', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%atomic_number, 26, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%valence, 8, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(2)%potential%center_band_s_up),  12.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_p_up), -22.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_d_up),  32.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_s_dw), -42.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_p_dw),  52.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_d_dw), -62.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_s_up),   72.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_p_up),  -82.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_d_up),   92.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_s_dw), -102.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_p_dw),  112.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_d_dw), -122.5, __FILE__, __LINE__, test_suite)


    call assert_equal(ar(3)%element%symbol, 'Fe', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%atomic_number, 26, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%valence, 8, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(3)%potential%center_band_s_up),  13.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_p_up), -23.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_d_up),  33.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_s_dw), -43.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_p_dw),  53.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_d_dw), -63.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_s_up),   73.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_p_up),  -83.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_d_up),   93.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_s_dw), -103.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_p_dw),  113.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_d_dw), -123.5, __FILE__, __LINE__, test_suite)

    

    call assert_equal(ar(4)%element%symbol, 'Fe', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%atomic_number, 26, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%valence, 8, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%num_quant_d, 5, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(4)%potential%center_band_s_up),  12.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_p_up), -22.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_d_up),  32.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_s_dw), -42.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_p_dw),  52.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_d_dw), -62.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_s_up),   72.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_p_up),  -82.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_d_up),   92.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_s_dw), -102.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_p_dw),  112.5, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_d_dw), -122.5, __FILE__, __LINE__, test_suite)


    ! (4) Finish your test with the two following lines
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
