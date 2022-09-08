!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Namelist Generator
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
!> lorem ipsum
!------------------------------------------------------------------------------


module namelist_generator_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use precision_mod, only: rp
  use string_mod, only: sl, int2str, real2str
  implicit none

  private

  type :: integer_variable
    character(len=sl) :: name
    integer :: value
    contains
      procedure :: to_string => integer_variable_to_string
  end type integer_variable
  
  type :: real_variable
    character(len=sl) :: name
    real(rp) :: value
    contains
      procedure :: to_string => real_variable_to_string
  end type real_variable
  
  type :: string_variable
    character(len=sl) :: name
    character(len=sl) :: value
    contains
      procedure :: to_string => string_variable_to_string
  end type string_variable

  type :: logical_variable
    character(len=sl) :: name
    logical :: value
    contains
      procedure :: to_string => logical_variable_to_string
  end type logical_variable

  type :: array_real_variable
    character(len=sl) :: name
    real(rp), dimension(:), allocatable :: value
    contains
      procedure :: to_string => array_real_variable_to_string
  end type array_real_variable


  !> Module's main structure
  type, public :: namelist_generator
    !> Name of namelist
    !> 
    !> Name of namelist
    character(len=sl) :: name

    !> Hold all variables
    !> 
    !> Hold all variables
    type(integer_variable), dimension(:), allocatable :: list_of_integer_variables
    type(real_variable), dimension(:), allocatable :: list_of_real_variables
    type(string_variable), dimension(:), allocatable :: list_of_string_variables
    type(logical_variable), dimension(:), allocatable :: list_of_logical_variables
    type(array_real_variable), dimension(:), allocatable :: list_of_array_real_variables
  contains
    procedure :: add_integer_variable, add_real4_variable, add_real8_variable, add_string_variable, add_logical_variable, add_array_real4_variable, add_array_real8_variable
    generic :: add_variable => add_integer_variable, add_real4_variable, add_real8_variable, add_string_variable, add_logical_variable, add_array_real4_variable, add_array_real8_variable
    procedure :: generate_namelist
    final :: destructor
  end type namelist_generator

  interface namelist_generator
    procedure :: constructor
  end interface namelist_generator

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] fname Input file with 'namelist_generator' namelist
  !> @return type(namelist_generator)
  !---------------------------------------------------------------------------
  function constructor(name) result(obj)
    type(namelist_generator) :: obj
    character(len=*), intent(in) :: name
    obj%name = name
    allocate(obj%list_of_integer_variables(0))
    allocate(obj%list_of_real_variables(0))
    allocate(obj%list_of_string_variables(0))
    allocate(obj%list_of_logical_variables(0))
    allocate(obj%list_of_array_real_variables(0))
  end function constructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine destructor(this)
    type(namelist_generator) :: this
  end subroutine destructor

  ! Member functions

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine add_integer_variable(this,name,value)
    implicit none
    class(namelist_generator),intent(inout) :: this
    type(integer_variable) :: variable
    character(len=*) :: name
    integer,intent(in) :: value
    type(integer_variable), dimension(:), allocatable :: list_of_integer_variables
    integer :: new_size

    variable%name = name
    variable%value = value

    if(allocated(this%list_of_integer_variables)) then
      call move_alloc(this%list_of_integer_variables,list_of_integer_variables)
    endif
    new_size = size(this%list_of_integer_variables) + 1
    allocate(this%list_of_integer_variables(new_size))

    if(new_size > 1) then
      this%list_of_integer_variables(:new_size-1) = list_of_integer_variables
    endif
    this%list_of_integer_variables(new_size) = variable
  end subroutine add_integer_variable
  
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine add_real4_variable(this,name,value)
    implicit none
    class(namelist_generator),intent(inout) :: this
    type(real_variable) :: variable
    character(len=*) :: name
    real,intent(in) :: value
    type(real_variable), dimension(:), allocatable :: list_of_real_variables
    integer :: new_size

    variable%name = name
    variable%value = value

    if(allocated(this%list_of_real_variables)) then
      call move_alloc(this%list_of_real_variables,list_of_real_variables)
    endif
    new_size = size(this%list_of_real_variables) + 1
    allocate(this%list_of_real_variables(new_size))

    if(new_size > 1) then
      this%list_of_real_variables(:new_size-1) = list_of_real_variables
    endif
    this%list_of_real_variables(new_size) = variable
  end subroutine add_real4_variable

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine add_real8_variable(this,name,value)
    implicit none
    class(namelist_generator),intent(inout) :: this
    type(real_variable) :: variable
    character(len=*) :: name
    real(rp),intent(in) :: value
    type(real_variable), dimension(:), allocatable :: list_of_real_variables
    integer :: new_size

    variable%name = name
    variable%value = value

    if(allocated(this%list_of_real_variables)) then
      call move_alloc(this%list_of_real_variables,list_of_real_variables)
    endif
    new_size = size(this%list_of_real_variables) + 1
    allocate(this%list_of_real_variables(new_size))

    if(new_size > 1) then
      this%list_of_real_variables(:new_size-1) = list_of_real_variables
    endif
    this%list_of_real_variables(new_size) = variable
  end subroutine add_real8_variable

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine add_string_variable(this,name,value)
    implicit none
    class(namelist_generator),intent(inout) :: this
    type(string_variable) :: variable
    character(len=*) :: name
    character(len=*),intent(in) :: value
    type(string_variable), dimension(:), allocatable :: list_of_string_variables
    integer :: new_size

    variable%name = name
    variable%value = value

    if(allocated(this%list_of_string_variables)) then
      call move_alloc(this%list_of_string_variables,list_of_string_variables)
    endif
    new_size = size(this%list_of_string_variables) + 1
    allocate(this%list_of_string_variables(new_size))

    if(new_size > 1) then
      this%list_of_string_variables(:new_size-1) = list_of_string_variables
    endif
    this%list_of_string_variables(new_size) = variable
  end subroutine add_string_variable

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine add_logical_variable(this,name,value)
    implicit none
    class(namelist_generator),intent(inout) :: this
    type(logical_variable) :: variable
    character(len=*) :: name
    logical,intent(in) :: value
    type(logical_variable), dimension(:), allocatable :: list_of_logical_variables
    integer :: new_size

    variable%name = name
    variable%value = value

    if(allocated(this%list_of_logical_variables)) then
      call move_alloc(this%list_of_logical_variables,list_of_logical_variables)
    endif
    new_size = size(this%list_of_logical_variables) + 1
    allocate(this%list_of_logical_variables(new_size))

    if(new_size > 1) then
      this%list_of_logical_variables(:new_size-1) = list_of_logical_variables
    endif
    this%list_of_logical_variables(new_size) = variable
  end subroutine add_logical_variable

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine add_array_real4_variable(this,name,value)
    implicit none
    class(namelist_generator),intent(inout) :: this
    type(array_real_variable) :: variable
    character(len=*) :: name
    real,dimension(:),intent(in) :: value
    type(array_real_variable), dimension(:), allocatable :: list_of_array_real_variables
    integer :: new_size

    variable%name = name
    allocate(variable%value(size(value)))
    variable%value = value

    if(allocated(this%list_of_array_real_variables)) then
      call move_alloc(this%list_of_array_real_variables,list_of_array_real_variables)
    endif
    new_size = size(this%list_of_array_real_variables) + 1
    allocate(this%list_of_array_real_variables(new_size))

    if(new_size > 1) then
      this%list_of_array_real_variables(:new_size-1) = list_of_array_real_variables
    endif
    this%list_of_array_real_variables(new_size) = variable
  end subroutine add_array_real4_variable

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine add_array_real8_variable(this,name,value)
    implicit none
    class(namelist_generator),intent(inout) :: this
    type(array_real_variable) :: variable
    character(len=*) :: name
    real(rp),dimension(:),intent(in) :: value
    type(array_real_variable), dimension(:), allocatable :: list_of_array_real_variables
    integer :: new_size

    variable%name = name
    allocate(variable%value(size(value)))
    variable%value = value

    if(allocated(this%list_of_array_real_variables)) then
      call move_alloc(this%list_of_array_real_variables,list_of_array_real_variables)
    endif
    new_size = size(this%list_of_array_real_variables) + 1
    allocate(this%list_of_array_real_variables(new_size))

    if(new_size > 1) then
      this%list_of_array_real_variables(:new_size-1) = list_of_array_real_variables
    endif
    this%list_of_array_real_variables(new_size) = variable
  end subroutine add_array_real8_variable

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine generate_namelist(this,fname)
    implicit none
    class(namelist_generator),intent(in) :: this
    character(len=*) :: fname
    integer :: funit, i

    open(unit=funit,file=trim(fname)//'.nml',action='write')

    write(funit,'("&",A)') this%name

    do i=1, size(this%list_of_integer_variables)
      write(funit,'(" ",A)') trim(this%list_of_integer_variables(i)%to_string())
    enddo
    do i=1, size(this%list_of_string_variables)
      write(funit,'(" ",A)') trim(this%list_of_string_variables(i)%to_string())
    enddo
    do i=1, size(this%list_of_real_variables)
      write(funit,'(" ",A)') trim(this%list_of_real_variables(i)%to_string())
    enddo
    do i=1, size(this%list_of_logical_variables)
      write(funit,'(" ",A)') trim(this%list_of_logical_variables(i)%to_string())
    enddo
    do i=1, size(this%list_of_array_real_variables)
      write(funit,'(" ",A)') trim(this%list_of_array_real_variables(i)%to_string())
    enddo

    write(funit,'("/")')

    close(funit)
  end subroutine generate_namelist

  function integer_variable_to_string(this) result(output)
    implicit none
    class(integer_variable),intent(in) :: this
    character(len=sl) :: output
    output = trim(this%name)//' = '//trim(int2str(this%value))
  end function integer_variable_to_string

  function string_variable_to_string(this) result(output)
    implicit none
    class(string_variable),intent(in) :: this
    character(len=sl) :: output
    output = trim(this%name)//" = '"//trim(this%value)//"'"
  end function string_variable_to_string

  function real_variable_to_string(this) result(output)
    implicit none
    class(real_variable),intent(in) :: this
    character(len=sl) :: output
    output = trim(this%name)//' = '//trim(real2str(this%value))
  end function real_variable_to_string

  function logical_variable_to_string(this) result(output)
    implicit none
    class(logical_variable),intent(in) :: this
    character(len=sl) :: output
    if(this%value) then
      output = trim(this%name)//' = .True.'
    else
      output = trim(this%name)//' = .False.'
    endif
  end function logical_variable_to_string

  function array_real_variable_to_string(this) result(output)
    implicit none
    class(array_real_variable),intent(in) :: this
    character(len=sl) :: output
    integer i
    output = ''
    do i=1,size(this%value)
      output = trim(output)//', '//trim(real2str(this%value(i)))
    enddo
    output = trim(this%name)//' = '//output(3:)
  end function array_real_variable_to_string

end module namelist_generator_mod
  
  
  