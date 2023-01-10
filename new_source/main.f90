program main

  use calculation_mod
  use os_mod
  use, intrinsic :: iso_fortran_env, only: output_unit
  use precision_mod, only: rp
  use logger_mod, only: g_logger
  use timer_mod, only: g_timer, timer

  implicit none
  type(calculation) :: calculation_obj
  type(argument_parser) :: args
  integer :: nomp
  real(rp) :: start, finish
#ifdef OMP
  ! External functions
  integer, external :: omp_get_num_threads
#endif


#ifdef OMP
   !$omp parallel 
   !$omp master
   nomp=omp_get_num_threads()
   !$omp end master
   !$omp end parallel
#endif
  g_timer = timer()
  !call g_logger%debug('Initializing with DEBUG=ON',__FILE__,__LINE__)

  ! Input
  args = argument_parser()
  calculation_obj = calculation(args%input)


  ! Run
  call g_timer%start('Calculation')
  call cpu_time(start)
  call calculation_obj%process
  call cpu_time(finish)
  call g_timer%stop('Calculation')
  print '("Total calculatio time = ",f10.3," seconds.")',(finish-start)!/32

  call g_timer%print_report()
end program main
