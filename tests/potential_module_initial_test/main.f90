program Main
    use potential_mod
    use unit_test

    type(test_suite_type) :: test_suite

    type(potential) :: pot
    type(potential), dimension(4) :: ar
    
    ! (1) Create a test
    call test_suite_init('Potential module initial test', test_suite)

    ! (2) Create as many cases as you need
    call test_case_create('Reading from current directory', test_suite)
    
    ! (3) Test your variables:
    
    pot = potential('Co','./')
    
    call assert_equal(real(pot%center_band_p_up),  0.448747850, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%center_band_d_up), -0.169678507, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%center_band_s_dw), -0.253518500, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%center_band_p_dw),  0.492954098, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%center_band_d_dw),  -0.056204107, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%width_band_s_up), 0.421162577, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%width_band_p_up), 0.275865043, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%width_band_d_up), 0.122886151, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%width_band_s_dw), 0.425030588, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%width_band_p_dw), 0.283340379, __FILE__, __LINE__, test_suite)
    call assert_equal(real(pot%width_band_d_dw), 0.137919686, __FILE__, __LINE__, test_suite)
    
    call test_case_create('Loading arrays of potential', test_suite)

    ar = array_of_potentials((/'Fe','Co','Cu','Ir'/),'.')

    call assert_equal(real(ar(1)%center_band_s_up), -0.304811719, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%center_band_p_up),  0.342417756, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%center_band_d_up), -0.213141873, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%center_band_s_dw), -0.277714747, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%center_band_p_dw),  0.397445417, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%center_band_d_dw),  0.012723659, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%width_band_s_up), 0.399729997, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%width_band_p_up), 0.260669836, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%width_band_d_up), 0.117641990, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%width_band_s_dw), 0.403198147, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%width_band_p_dw), 0.268128801, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(1)%width_band_d_dw), 0.137180112, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(2)%center_band_s_up), -0.266523012, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%center_band_p_up),  0.448747850, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%center_band_d_up), -0.169678507, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%center_band_s_dw), -0.253518500, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%center_band_p_dw),  0.492954098, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%center_band_d_dw),  -0.056204107, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%width_band_s_up), 0.421162577, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%width_band_p_up), 0.275865043, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%width_band_d_up), 0.122886151, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%width_band_s_dw), 0.425030588, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%width_band_p_dw), 0.283340379, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(2)%width_band_d_dw), 0.137919686, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(3)%center_band_s_up), -0.430414987, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%center_band_p_up),  0.249868163, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%center_band_d_up), -0.299825466, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%center_band_s_dw), -0.430414983, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%center_band_p_dw),  0.249868170, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%center_band_d_dw),  -0.299825467, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%width_band_s_up), 0.391394572, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%width_band_p_up), 0.260459319, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%width_band_d_up), 0.095901360, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%width_band_s_dw), 0.391394573, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%width_band_p_dw), 0.260459319, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(3)%width_band_d_dw), 0.095901360, __FILE__, __LINE__, test_suite)

    call assert_equal(real(ar(4)%center_band_s_up), -0.463475546, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%center_band_p_up),  0.300925645, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%center_band_d_up), -0.209816425, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%center_band_s_dw), -0.463475546, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%center_band_p_dw),  0.300925645, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%center_band_d_dw),  -0.209816425, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%width_band_s_up), 0.387432631, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%width_band_p_up), 0.241564895, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%width_band_d_up), 0.168137666, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%width_band_s_dw), 0.387432631, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%width_band_p_dw), 0.241564895, __FILE__, __LINE__, test_suite)
    call assert_equal(real(ar(4)%width_band_d_dw), 0.168137666, __FILE__, __LINE__, test_suite)
    
    ! (4) Finish your test with the two following lines
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
