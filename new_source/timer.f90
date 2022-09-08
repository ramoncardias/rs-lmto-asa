!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Timer
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to mesure timer between code
!------------------------------------------------------------------------------


module timer_mod
  use, intrinsic :: iso_fortran_env, only: output_unit
  use string_mod, only: sl, clean_str
  use precision_mod, only: rp
  use logger_mod, only: g_logger
  implicit none

  ! private

  type, public :: timer
    character(len=sl) :: label
    real, private :: t_start, t_finish
  contains
    procedure :: start => timer_start
    procedure :: finish => timer_finish
    procedure :: has_finished => timer_has_finished
    procedure :: delta => timer_delta
    procedure :: report => timer_report
  endtype timer

  type, public :: time_reporter
    type(timer), allocatable, dimension(:), private :: timer_list
    type(timer), private :: t_global
    contains
    procedure :: start => time_reporter_start
    procedure :: finish => time_reporter_finish
    procedure :: report => time_reporter_report
    final :: time_reporter_destructor
  endtype time_reporter

  interface time_reporter
    procedure :: time_reporter_constructor
  end interface time_reporter
  
  ! Global
  type(time_reporter) :: g_timer

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> time_reporter_constructor
  !
  !> @return type(time_reporter)
  !---------------------------------------------------------------------------
  function time_reporter_constructor() result(obj)
    type(time_reporter) :: obj
    call obj%t_global%start('global')
    allocate(obj%timer_list(0))
  end function time_reporter_constructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine time_reporter_destructor(this)
    type(time_reporter) :: this
    if (allocated(this%timer_list)) then
      print *, 'dealocando'
      deallocate(this%timer_list)
    endif
  end subroutine time_reporter_destructor


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Measures the start time
  !---------------------------------------------------------------------------
  subroutine timer_start(this,label)
    class(timer) :: this
    character(len=*), intent(in) :: label
    this%label = label
    call cpu_time(this%t_start)
    this%t_finish = this%t_start
  end subroutine timer_start

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Measures the final time
  !---------------------------------------------------------------------------
  subroutine timer_finish(this)
    class(timer) :: this
    call cpu_time(this%t_finish)
  end subroutine timer_finish

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Reports the time spent with label
  !---------------------------------------------------------------------------
  subroutine timer_report(this)
    class(timer) :: this
    call g_logger%log(this%label,this%delta(),__FILE__,__LINE__)
  end subroutine timer_report

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Test if t_finish .gt. t_start
  !---------------------------------------------------------------------------
  function timer_has_finished(this) result(has_finished)
    class(timer) :: this
    logical :: has_finished
    has_finished = this%t_finish .gt. this%t_start
  end function timer_has_finished

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Returns the difference between t_finish and t_start
  !---------------------------------------------------------------------------
  function timer_delta(this) result(delta)
    class(timer) :: this
    real(rp) :: delta
    delta = this%t_finish-this%t_start
  end function timer_delta

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Measures the start time
  !---------------------------------------------------------------------------
  subroutine time_reporter_start(this,label)
    class(time_reporter) :: this
    character(len=*), intent(in) :: label
    type(timer), allocatable, dimension(:) :: new_timer_list
    type(timer) :: new_timer
    integer :: new_size
    new_size = size(this%timer_list) + 1
    allocate(new_timer_list(new_size))
    new_timer_list(:new_size-1) = this%timer_list
    call new_timer%start(label)
    new_timer_list(new_size) = new_timer
    call move_alloc(new_timer_list,this%timer_list)
  end subroutine time_reporter_start

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Measures the final time
  !---------------------------------------------------------------------------
  subroutine time_reporter_finish(this,label)
    class(time_reporter) :: this
    character(len=*), intent(in) :: label
    logical :: has_found
    integer :: i
    has_found = .False.
    do i = 1, size(this%timer_list)
      if (trim(this%timer_list(i)%label) == trim(label) .and. .not. this%timer_list(i)%has_finished()) then
        call this%timer_list(i)%finish
        has_found = .True.
        continue
      endif
    enddo
    if (.not. has_found) then
      call g_logger%error('Finishing '//trim(label)//' without starting first',__FILE__,__LINE__)
    endif
  end subroutine time_reporter_finish

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Reports the time spent with label
  !---------------------------------------------------------------------------
  subroutine time_reporter_report(this)
    class(time_reporter) :: this
    integer :: i, j, labels_list_size
    character(len=sl), dimension(:), allocatable :: labels_list, aux_labels_list
    integer :: n_unique_labels
    real(rp) :: t_mean, t_max, t_min, t_delta, total_time
    integer :: t_total
    call this%t_global%finish
    total_time = this%t_global%delta()
    allocate(labels_list(0))
    labels_list_size = 0
    print *,'--------------------------------- TIME CONSUPTION REPORT --------------------------------'
    print *, 'CATEGORY                          Counts   Max (sec)   Min (sec)  Mean (sec)    Mean (%)'
    do i = 1, size(this%timer_list)
      if (all(labels_list .ne. this%timer_list(i)%label)) then
        call move_alloc(labels_list, aux_labels_list)
        labels_list_size = labels_list_size + 1
        allocate(labels_list(labels_list_size))
        labels_list(:labels_list_size-1) = aux_labels_list 
        labels_list(labels_list_size) = this%timer_list(i)%label
        deallocate(aux_labels_list)
      endif
    enddo
    do j = 1, labels_list_size
      t_total = 0
      t_max = 0
      t_min = 1e23
      do i = 1, size(this%timer_list)
        if (labels_list(j) == this%timer_list(i)%label) then
          t_delta = this%timer_list(i)%delta()
          t_mean = t_mean + t_delta
          t_total = t_total + 1
          t_max = max(t_max,t_delta)
          t_min = min(t_min,t_delta)
        endif
      enddo
      t_mean = t_mean / t_total
      print '("  ",A28,"  ",I8,"  ",E10.3,"  ",E10.3,"  ",E10.3,"  ",F10.2)', clean_str(labels_list(j)), t_total, t_max, t_min, t_mean, 100*t_mean/total_time
    enddo
#ifdef DEBUG
    print *,'-------------------------- FULL DESCRIPTION OF TIME CONSUPTION --------------------------'
    do i = 1, size(this%timer_list)
      call this%timer_list(i)%report()
    enddo
#endif
    print *,'-----------------------------------------------------------------------------------------'
    print '(A,E10.3)','Total CPU TIME   ',total_time
  end subroutine time_reporter_report

end module timer_mod


