program Main
    ! use globals_mod, only : g_timer
    use timer_mod, only: time_reporter

    type(time_reporter) :: tr
    integer :: i,j,sum

    tr = time_reporter()

    ! (3) Test your variables:
    do j=0,2
        call tr%start('Selfcon.calculation_1')
        sum=0
        do i=0,1000000000
            if (mod(i,2) == 0) then
                sum = sum/2
            else
                sum = 3*sum+1
            endif
        enddo
        call tr%finish('Selfcon.calculation_1')

        call tr%start('Selfcon.calculation_2')
        sum=0
        do i=0,1000000000
            sum = sum + i
        enddo
        call tr%finish('Selfcon.calculation_2')
    enddo
    call tr%report()

end program Main
