!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Mix
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


module mix_mod
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use precision_mod, only: rp
    implicit none
  
    private
  
    !> Module's main structure
    type, public :: mix
      !> Description
      !> 
      !> Description
      integer :: var
    contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state
      final :: destructor
    end type mix
  
    interface mix
      procedure :: constructor
    end interface mix
  
  contains
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Constructor
    !
    !> @param[in] fname Input file with 'mix' namelist
    !> @return type(mix)
    !---------------------------------------------------------------------------
    function constructor(fname) result(obj)
      type(mix) :: obj
      character(len=*), intent(in) :: fname
  
      call obj%restore_to_default()
      call obj%build_from_file(fname)
    end function constructor
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Destructor
    !---------------------------------------------------------------------------
    subroutine destructor(this)
      type(mix) :: this
    end subroutine destructor
  
    ! Member functions
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Read parameters from input file
    !
    !> @param[in] fname Namelist file
    !---------------------------------------------------------------------------
    subroutine build_from_file(this, fname)
      class(mix),intent(inout) :: this
      character(len=*), intent(in) :: fname
  
      integer :: var

      ! variables associated with the reading processes
      integer :: iostatus, funit



      namelist /mix/ var
  
  
      var = this%var
  
      open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
      if(iostatus /= 0) then
        write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
        error stop
      endif
  
      read(funit,nml=mix,iostat=iostatus)
      if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
        write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
        write(error_unit,'("iostatus = ",I0)') iostatus
      endif
      close(funit)
  
      this%var = 0
  
    end subroutine build_from_file
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
      implicit none
      class(mix), intent(out) :: this
  
      this%var = 0
      
    end subroutine restore_to_default
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Print class members values in namelist format 
    !>
    !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
    !> @param[in] unit File unit used to write namelist
    !> @param[in] file File name used to write namelist
    !---------------------------------------------------------------------------
    subroutine print_state(this,unit,file)
      implicit none
      class(mix), intent(in) :: this
  
      integer,intent(in),optional :: unit
      character(len=*),intent(in),optional :: file
      integer :: newunit
  
      integer :: var
  
      namelist /mix/ var
  
      var = this%var
  
      if(present(unit) .and. present(file)) then
        write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
        error stop
      else if(present(unit)) then
          write(unit,nml=mix)
      else if(present(file)) then
          open(unit=newunit,file=file)
          write(newunit,nml=mix)
          close(newunit)
      else
          write(*,nml=mix)
      endif
      close(newunit)
    end subroutine print_state
  
  end module mix_mod