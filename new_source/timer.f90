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
!> Module to mesure time between code
!------------------------------------------------------------------------------


module timer_mod
  use string_mod, only: sl, clean_str, split
  use precision_mod, only: rp
  use logger_mod, only: g_logger
  use report_mod
  implicit none

  ! Private by default
  private

  type :: simple_timer
    character(len=sl) :: label
    real, private :: t_start, t_finish
  contains
    procedure :: start => simple_timer_start
    procedure :: finish => simple_timer_finish
    procedure :: has_finished => simple_timer_has_finished
    procedure :: delta => simpletimer_delta
    procedure :: report => simple_timer_report
  endtype simple_timer

  type :: timer
    type(simple_timer), allocatable, dimension(:), private :: timer_list
    type(simple_timer), private :: t_global
    contains
    procedure :: start => timer_start
    procedure :: finish => timer_finish
    procedure :: print_report => timer_report
    final :: timer_destructor
  endtype timer

  interface timer
    procedure :: timer_constructor
  end interface timer
  
  
  ! Global
  type(timer) :: g_timer
  
  ! Public
  public :: timer, g_timer

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> timer_constructor
  !
  !> @return type(timer)
  !---------------------------------------------------------------------------
  function timer_constructor() result(obj)
    type(timer) :: obj
    call obj%t_global%start('global')
    allocate(obj%timer_list(0))
  end function timer_constructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine timer_destructor(this)
    type(timer) :: this
    if (allocated(this%timer_list)) then
      print *, 'dealocando'
      deallocate(this%timer_list)
    endif
  end subroutine timer_destructor


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Measures the start time
  !---------------------------------------------------------------------------
  subroutine simple_timer_start(this,label)
    class(simple_timer) :: this
    character(len=*), intent(in) :: label
    this%label = label
    call cpu_time(this%t_start)
    this%t_finish = this%t_start
  end subroutine simple_timer_start

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Measures the final time
  !---------------------------------------------------------------------------
  subroutine simple_timer_finish(this)
    class(simple_timer) :: this
    call cpu_time(this%t_finish)
  end subroutine simple_timer_finish

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Reports the time spent with label
  !---------------------------------------------------------------------------
  subroutine simple_timer_report(this)
    class(simple_timer) :: this
    call g_logger%log(this%label,this%delta(),__FILE__,__LINE__)
  end subroutine simple_timer_report

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Test if t_finish .gt. t_start
  !---------------------------------------------------------------------------
  function simple_timer_has_finished(this) result(has_finished)
    class(simple_timer) :: this
    logical :: has_finished
    has_finished = this%t_finish .gt. this%t_start
  end function simple_timer_has_finished

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Returns the difference between t_finish and t_start
  !---------------------------------------------------------------------------
  function simpletimer_delta(this) result(delta)
    class(simple_timer) :: this
    real(rp) :: delta
    delta = this%t_finish-this%t_start
  end function simpletimer_delta

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Measures the start time
  !---------------------------------------------------------------------------
  subroutine timer_start(this,label)
    class(timer) :: this
    character(len=*), intent(in) :: label
    type(simple_timer), allocatable, dimension(:) :: new_timer_list
    type(simple_timer) :: new_timer
    integer :: new_size
    new_size = size(this%timer_list) + 1
    allocate(new_timer_list(new_size))
    new_timer_list(:new_size-1) = this%timer_list
    call new_timer%start(label)
    new_timer_list(new_size) = new_timer
    call move_alloc(new_timer_list,this%timer_list)
  end subroutine timer_start

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Measures the final time
  !---------------------------------------------------------------------------
  subroutine timer_finish(this,label)
    class(timer) :: this
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
  end subroutine timer_finish

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief Reports the time spent with label
  !---------------------------------------------------------------------------
  subroutine timer_report(this)
    class(timer) :: this
    type(report) :: rep
    integer :: i
    rep = report('total')
    do i=1,size(this%timer_list)
      call rep%add_value('total.'//trim(this%timer_list(i)%label),this%timer_list(i)%delta())
    enddo
    call rep%print_report('TIMER REPORT')
  end subroutine timer_report

end module timer_mod


