program Main
    use string_mod
    use unit_test

    type(test_suite_type) :: test_suite
    character(len=255) :: a
    
    call test_suite_init('TEST STRING MOD', test_suite)

    call test_case_create('STR_CONTAINS', test_suite)
    call assert_true(str_contains('ABCD','BC'), __FILE__, __LINE__, test_suite)
    call assert_false(str_contains('ABCD','ABD'), __FILE__, __LINE__, test_suite)
    call assert_true(startswith('ABCD','ABC'), __FILE__, __LINE__, test_suite)
    
    call test_case_create('STARTSWITH', test_suite)
    call assert_false(startswith('ABCD','BC'), __FILE__, __LINE__, test_suite)
    call assert_false(startswith('ABCD','ABD'), __FILE__, __LINE__, test_suite)
    call assert_true(startswith('ABCD','ABC'), __FILE__, __LINE__, test_suite)
    call assert_true(startswith('TOTAL.Calc_0.subcalc1','TOTAL'), __FILE__, __LINE__, test_suite)

    call test_case_create('ENDSWITH', test_suite)
    call assert_false(endswith('ABCD','BC'), __FILE__, __LINE__, test_suite)
    call assert_false(endswith('ABCD','ABD'), __FILE__, __LINE__, test_suite)
    call assert_false(endswith('ABCD','ABC'), __FILE__, __LINE__, test_suite)
    call assert_true(endswith('ABCD','CD'), __FILE__, __LINE__, test_suite)
    call assert_true(endswith('ABCD','BCD'), __FILE__, __LINE__, test_suite)

    call test_case_create('JOIN', test_suite)
    call assert_equal(fjoin(['A','B','C','D']), 'ABCD', __FILE__, __LINE__, test_suite)
    call assert_equal(fjoin(['A','B','C','D'],'.'), 'A.B.C.D', __FILE__, __LINE__, test_suite)
    call assert_equal(fjoin(['TOTAL','ABCD ','CALC1'],'-'), 'TOTAL-ABCD-CALC1', __FILE__, __LINE__, test_suite)

    call test_case_create('REPLACE', test_suite)
    call assert_equal(freplace('ABCD','A','1'), '1BCD', __FILE__, __LINE__, test_suite)
    call assert_equal(freplace('ABCDABCD','A','1'), '1BCD1BCD', __FILE__, __LINE__, test_suite)
    call assert_equal(freplace('AAAAABBB','B','2'), 'AAAAA222', __FILE__, __LINE__, test_suite)
    call assert_equal(freplace('abcABCabc','a','123'), '123bcABC123bc', __FILE__, __LINE__, test_suite)
    call assert_equal(freplace('abcABCabc','abc','3'), '3ABC3', __FILE__, __LINE__, test_suite)
    
    call test_case_create('PATH_JOIN', test_suite)
    call assert_equal(path_join(['DIR1','DIR2','DIR3']), 'DIR1/DIR2/DIR3', __FILE__, __LINE__, test_suite)
    call assert_equal(path_join(['DIR1/','DIR2 ','DIR3 ']), 'DIR1/DIR2/DIR3', __FILE__, __LINE__, test_suite)
    call assert_equal(path_join(['DIR1/','DIR2/','DIR3 ']), 'DIR1/DIR2/DIR3', __FILE__, __LINE__, test_suite)
    call assert_equal(path_join(['DIR1/','DIR2/','DIR3/']), 'DIR1/DIR2/DIR3', __FILE__, __LINE__, test_suite)
    
    call test_case_create('SPLIT', test_suite)
    call assert_equal(fsplit('ABC','.'), ['ABC'], __FILE__, __LINE__, test_suite)
    call assert_equal(fsplit('ABC.','.'), ['ABC','   '], __FILE__, __LINE__, test_suite)

    call test_case_create('JOIN SPLIT', test_suite)
    call assert_equal(fjoin(fsplit('ABC.DEF.GHI','.'),'-'), 'ABC-DEF-GHI', __FILE__, __LINE__, test_suite)
    
    call test_case_create('FCOUNT_STR', test_suite)
    call assert_equal(fcount_str('ABC.DEF.GHI','.'), 2, __FILE__, __LINE__, test_suite)
    call assert_equal(fcount_str('....','.'), 4, __FILE__, __LINE__, test_suite)
    call assert_equal(fcount_str('ABABDABB.','AB'), 3, __FILE__, __LINE__, test_suite)
    call assert_equal(fcount_str('ABABDABB.','D'), 1, __FILE__, __LINE__, test_suite)
    call assert_equal(fcount_str('ABABDABB.','B'), 4, __FILE__, __LINE__, test_suite)

    call test_case_create('FINDENT_STR', test_suite)
    call assert_equal(findent_str('ASD',0),'ASD', __FILE__, __LINE__, test_suite)
    call assert_equal(findent_str('123',1),' 123', __FILE__, __LINE__, test_suite)
    call assert_equal(findent_str('   FILE   ',2),'     FILE   ', __FILE__, __LINE__, test_suite)
    call assert_equal(findent_str('ABC DEF',3),'   ABC DEF', __FILE__, __LINE__, test_suite)

    
    call test_suite_report(test_suite)
    call test_suite_final(test_suite)

end program Main
