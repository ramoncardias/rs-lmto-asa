program Main
    use self_mod
    use control_mod
    use unit_test

    type(test_suite_type) :: test_suite_control, test_suite_self_default, test_suite_self_from_file

    call test_suite_init('Control input test', test_suite_control)


    control_obj = control('input.nml')

    call test_case_create('Control Module from input.nml')

    call assert_equal(control_obj%nrec     ,2,   __FILE__, __LINE__, suite=test_suite_control)
    call assert_equal(control_obj%calctype ,'B', __FILE__, __LINE__, suite=test_suite_control)
    call assert_equal(control_obj%ntot     ,10,  __FILE__, __LINE__, suite=test_suite_control)
  
    call test_suite_report(test_suite_control)
    call test_suite_final(test_suite_control)

    
    
    call test_suite_init('SELF MODULE TEST [Default]', test_suite_self_default)

    call self_obj%restore_to_default()
    
    call test_case_create('Control variables')

    call assert_true(self_obj%all_inequivalent, __FILE__, __LINE__, suite=test_suite_self_default)

    call test_case_create('Variables according to calculation type')

    call assert_equal(self_obj%channels_ldos,    6000, __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_equal(self_obj%nbas,             control_obj%NTOT, __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_approximate(self_obj%energy_min, real(-5.5,8), __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_approximate(self_obj%energy_max, real( 5.5,8), __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_approximate(self_obj%fermi,      real(-0.05,8), __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_false(self_obj%fix_fermi, __FILE__, __LINE__, suite=test_suite_self_default)

    call test_case_create('Wigner Seitz Radius')
    
    call assert_true(self_obj%ws_all, __FILE__, __LINE__, suite=test_suite_self_default)

    call assert_equal(size(self_obj%ws), control_obj%nrec, __FILE__, __LINE__, suite=test_suite_self_default)
    do i=1, size(self_obj%ws)
        val = real(control_obj%wav * 1.8897259886,4)
        call assert_approximate(real(self_obj%ws(i),4),  val, __FILE__, __LINE__, eps=0.01, suite=test_suite_self_default)
    enddo

    
    call test_case_create('Mixing parameters')
    
    call assert_true(self_obj%mix_all, __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_equal(size(self_obj%mix),control_obj%nrec, __FILE__, __LINE__, suite=test_suite_self_default)
    do i=1, size(self_obj%mix)
        call assert_approximate(real(self_obj%mix(i),4),  0.01, __FILE__, __LINE__, suite=test_suite_self_default)
    enddo

    
    call test_case_create('Magnetic mixing parameters')
    
    call assert_false(self_obj%magnetic_mixing, __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_true(self_obj%mixmag_all, __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_equal(size(self_obj%mixmag),control_obj%nrec, __FILE__, __LINE__, suite=test_suite_self_default)
    do i=1, size(self_obj%mixmag)
        call assert_approximate(real(self_obj%mixmag(i),4),  0.05, __FILE__, __LINE__, suite=test_suite_self_default)
    enddo


    call test_case_create('Convergence parameters')
    
    call assert_approximate(self_obj%conv_thr, 0.5d-10,__FILE__, __LINE__, suite=test_suite_self_default)
    call assert_equal(self_obj%nstep, 1,__FILE__, __LINE__, suite=test_suite_self_default)

    call test_case_create('Constrains variables')
    
    call assert_false(self_obj%freeze, __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_false(self_obj%rigid_band, __FILE__, __LINE__, suite=test_suite_self_default)
    call assert_equal(size(self_obj%rb),control_obj%nrec,__FILE__, __LINE__, suite=test_suite_self_default)
    do i=1, size(self_obj%rb)
        call assert_equal(self_obj%rb(i),  2, __FILE__, __LINE__, suite=test_suite_self_default)
    enddo

    call assert_false(self_obj%orbital_polarization, __FILE__, __LINE__, suite=test_suite_self_default)

    
    call test_suite_report(test_suite_self_default)
    call test_suite_final(test_suite_self_default)
    
    
    call test_suite_init('SELF MODULE TEST [Default]', test_suite_self_from_file)
    
    self_obj = self('input.nml')

    call test_case_create('Control variables')
    call assert_false(self_obj%all_inequivalent, __FILE__, __LINE__, suite=test_suite_self_from_file)
  
    call test_case_create('Variables according to calculation type')
    call assert_equal(self_obj%channels_ldos, 1234, __FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_equal(self_obj%nbas, 1, __FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%energy_min, real(-3.0,8), __FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%energy_max, real( 3.0,8), __FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%fermi, real(-0.123,8), __FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_true(self_obj%fix_fermi, __FILE__, __LINE__,suite=test_suite_self_from_file)
  
    call test_case_create('Wigner Seitz Radius')
    call assert_false(self_obj%ws_all, __FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_equal(size(self_obj%ws), control_obj%nrec, __FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%ws(1),real(1.23,8), __FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%ws(2),real(4.56,8), __FILE__, __LINE__,suite=test_suite_self_from_file)
    ! self_obj%ws = 1.23, 4.56
  
    call test_case_create('Mixing parameters')
    call assert_false(self_obj%mix_all,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_equal(size(self_obj%mix),control_obj%nrec,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%mix(1),real(4.560,8),__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%mix(2),real(7.89,8),__FILE__, __LINE__,suite=test_suite_self_from_file)
  
    call test_case_create('Magnetic mixing parameters')
    call assert_true(self_obj%magnetic_mixing,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_false(self_obj%mixmag_all,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_equal(size(self_obj%mixmag), control_obj%nrec,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%mixmag(1), real(2.34,8),__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%mixmag(2), real(5.67,8),__FILE__, __LINE__,suite=test_suite_self_from_file)
  
    call test_case_create('Convergence parameters')
    call assert_equal(self_obj%nstep, 1000,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%conv_thr, 0.5d-10,__FILE__, __LINE__,suite=test_suite_self_from_file)
  
    call test_case_create('Constrains variables')
    call assert_true(self_obj%freeze,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_true(self_obj%rigid_band,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_equal(size(self_obj%rb), control_obj%nrec,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_equal(self_obj%rb(1), 123,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_equal(self_obj%rb(2), 456,__FILE__, __LINE__,suite=test_suite_self_from_file)
  
    call test_case_create('Other Variables')
    
    call assert_equal(self_obj%init, 123,__FILE__, __LINE__,suite=test_suite_self_from_file)
    call assert_approximate(self_obj%valence, real(8.123,8),__FILE__, __LINE__,suite=test_suite_self_from_file)

    call test_suite_report(test_suite_self_from_file)
    call test_suite_final(test_suite_self_from_file)
    

end program Main
