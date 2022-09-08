program Main
    use symbolic_atom_mod
    use unit_test

    type(test_suite_type) :: test_suite
    
    type(symbolic_atom) :: sym_atm
    type(symbolic_atom), dimension(4) :: ar

    ! (1) Create a test
    call test_suite_init('Symbolic atom module initial test', test_suite)

    ! (2) Create as many cases as you need
    call test_case_create('Reading from current directory', test_suite)

    ! (3) Test your variables:

    sym_atm = symbolic_atom('Co','./')

    call assert_equal(sym_atm%element%symbol, 'Co', __FILE__, __LINE__, test_suite)
    call assert_equal(sym_atm%element%atomic_number, 27, __FILE__, __LINE__, test_suite)
    call assert_equal(sym_atm%element%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(sym_atm%element%valence, 9, __FILE__, __LINE__, test_suite)
    call assert_equal(sym_atm%element%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(sym_atm%element%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(sym_atm%element%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(sym_atm%element%num_quant_d, 3, __FILE__, __LINE__, test_suite)
    
    call assert_equal(real(sym_atm%potential%center_band_p_up),  0.448747850, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%center_band_d_up), -0.169678507, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%center_band_s_dw), -0.253518500, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%center_band_p_dw),  0.492954098, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%center_band_d_dw),  -0.056204107, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%width_band_s_up), 0.421162577, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%width_band_p_up), 0.275865043, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%width_band_d_up), 0.122886151, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%width_band_s_dw), 0.425030588, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%width_band_p_dw), 0.283340379, __FILE__, __LINE__, test_suite)
    call assert_equal(real(sym_atm%potential%width_band_d_dw), 0.137919686, __FILE__, __LINE__, test_suite)

    call test_case_create('Loading arrays of symbolic_atoms', test_suite)

    ar = array_of_symbolic_atoms((/'Fe','Co','Cu','Ir'/),'.')

    call assert_equal(ar(1)%element%symbol, 'Fe', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%atomic_number, 26, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%valence, 8, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(1)%element%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(1)%potential%center_band_s_up), -0.304811719, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_p_up),  0.342417756, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_d_up), -0.213141873, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_s_dw), -0.277714747, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_p_dw),  0.397445417, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%center_band_d_dw),  0.012723659, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_s_up), 0.399729997, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_p_up), 0.260669836, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_d_up), 0.117641990, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_s_dw), 0.403198147, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_p_dw), 0.268128801, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%potential%width_band_d_dw), 0.137180112, __FILE__, __LINE__, test_suite)

    call assert_equal(ar(2)%element%symbol, 'Co', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%atomic_number, 27, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%valence, 9, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(2)%element%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(2)%potential%center_band_s_up), -0.266523012, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_p_up),  0.448747850, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_d_up), -0.169678507, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_s_dw), -0.253518500, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_p_dw),  0.492954098, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%center_band_d_dw),  -0.056204107, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_s_up), 0.421162577, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_p_up), 0.275865043, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_d_up), 0.122886151, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_s_dw), 0.425030588, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_p_dw), 0.283340379, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%potential%width_band_d_dw), 0.137919686, __FILE__, __LINE__, test_suite)

    call assert_equal(ar(3)%element%symbol, 'Cu', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%atomic_number, 29, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%core, 18, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%valence, 11, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%f_core, 0, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%num_quant_s, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%num_quant_p, 4, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(3)%element%num_quant_d, 3, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(3)%potential%center_band_s_up), -0.430414987, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_p_up),  0.249868163, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_d_up), -0.299825466, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_s_dw), -0.430414983, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_p_dw),  0.249868170, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%center_band_d_dw),  -0.299825467, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_s_up), 0.391394572, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_p_up), 0.260459319, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_d_up), 0.095901360, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_s_dw), 0.391394573, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_p_dw), 0.260459319, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%potential%width_band_d_dw), 0.095901360, __FILE__, __LINE__, test_suite)

    call assert_equal(ar(4)%element%symbol, 'Ir', __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%atomic_number, 77, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%core, 54, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%valence, 9, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%f_core, 14, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%num_quant_s, 6, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%num_quant_p, 6, __FILE__, __LINE__, test_suite)
    call assert_equal(ar(4)%element%num_quant_d, 5, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(4)%potential%center_band_s_up), -0.463475546, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_p_up),  0.300925645, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_d_up), -0.209816425, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_s_dw), -0.463475546, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_p_dw),  0.300925645, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%center_band_d_dw),  -0.209816425, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_s_up), 0.387432631, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_p_up), 0.241564895, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_d_up), 0.168137666, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_s_dw), 0.387432631, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_p_dw), 0.241564895, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%potential%width_band_d_dw), 0.168137666, __FILE__, __LINE__, test_suite)
    
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
