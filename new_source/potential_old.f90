!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Potential
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
!> Module to handle data relative to potential parameters of the
!Hamiltonian
!------------------------------------------------------------------------------

module potential_mod
  use lattice_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use precision_mod, only: rp
  use string_mod, only: sl
  implicit none
  
  private

  type, public :: potential

    !> Potential parameters filename
    character(len=sl),dimension(:),allocatable :: potfile

    !> Center of the band
    real(rp), dimension(:), allocatable :: center_band_s_up, center_band_p_up, center_band_d_up
    real(rp), dimension(:), allocatable :: center_band_s_dw, center_band_p_dw, center_band_d_dw
    !> Width of the band
    real(rp), dimension(:), allocatable :: width_band_s_up, width_band_p_up, width_band_d_up
    real(rp), dimension(:), allocatable :: width_band_s_dw, width_band_p_dw, width_band_d_dw


    !> Lattice class
    type(lattice), pointer :: lattice
  contains
    procedure :: build_from_file
    procedure :: restore_to_default
    procedure :: read_potential
    final :: destructor

  end type potential

  interface potential
    procedure :: constructor
  end interface potential

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] fname Namelist file
  !> @return type(pootential)
  !---------------------------------------------------------------------------
  function constructor(fname,lattice_obj) result(obj)
    type(potential) :: obj
    type(lattice), target, intent(in) :: lattice_obj
    character(len=*), intent(in) :: fname
    
    obj%lattice => lattice_obj
    
    call obj%restore_to_default()
    call obj%build_from_file(fname)
  end function constructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine destructor(this)
    type(potential) :: this
  end subroutine


  ! Member functions
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Read parameters from input file
  !
  !> @param[in] fname Namelist file
  !---------------------------------------------------------------------------
  subroutine build_from_file(this,fname)
    class(potential), intent(inout) :: this
    character(len=*), intent(in) :: fname
    ! Namelist variables
    character(len=sl),dimension(:),allocatable :: potfile    
    ! variables associated with the reading processes
    integer :: iostatus, funit

    ! Namelist
    namelist /potential/ potfile 

    open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
    if(iostatus /= 0) then
      write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
      error stop
    endif

    allocate(potfile(this%lattice%ntype))
    read(funit,nml=potential,iostat=iostatus)

    call move_alloc(potfile,this%potfile)

    call this%read_potential()
  end subroutine  build_from_file


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Reset all members to default
  !---------------------------------------------------------------------------
  subroutine restore_to_default(this)
    class(potential) :: this

    allocate(this%center_band_s_up(this%lattice%ntype),this%center_band_p_up(this%lattice%ntype),this%center_band_d_up(this%lattice%ntype))
    allocate(this%center_band_s_dw(this%lattice%ntype),this%center_band_p_dw(this%lattice%ntype),this%center_band_d_dw(this%lattice%ntype))
    allocate(this%width_band_s_up(this%lattice%ntype),this%width_band_p_up(this%lattice%ntype),this%width_band_d_up(this%lattice%ntype))
    allocate(this%width_band_s_dw(this%lattice%ntype),this%width_band_p_dw(this%lattice%ntype),this%width_band_d_dw(this%lattice%ntype))

  end subroutine restore_to_default


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Read the potential parameter from a given file
  !---------------------------------------------------------------------------
  subroutine read_potential(this)
    class(potential),intent(inout) :: this
    integer :: iostatus, funit, i, j
    !> Namelist variables
    real(rp) :: center_band_s_up, center_band_p_up, center_band_d_up
    real(rp) :: center_band_s_dw, center_band_p_dw, center_band_d_dw
    real(rp) :: width_band_s_up, width_band_p_up, width_band_d_up
    real(rp) :: width_band_s_dw, width_band_p_dw, width_band_d_dw

    namelist /par/ center_band_s_up, center_band_p_up, center_band_d_up,&
                   center_band_s_dw, center_band_p_dw, center_band_d_dw,&
                   width_band_s_up, width_band_p_up, width_band_d_up,&
                   width_band_s_dw, width_band_p_dw, width_band_d_dw


    do i=1,this%lattice%ntype
      open(newunit=funit,file=this%potfile(i),action='read',iostat=iostatus,status='old')
      if(iostatus /= 0) then
        write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
        error stop
      end if
      read(funit,nml=par,iostat=iostatus)

      this%center_band_s_up(i) = center_band_s_up
      this%center_band_p_up(i) = center_band_p_up
      this%center_band_d_up(i) = center_band_d_up

      this%center_band_s_dw(i) = center_band_s_dw
      this%center_band_p_dw(i) = center_band_p_dw
      this%center_band_d_dw(i) = center_band_d_dw

      this%width_band_s_up(i) = width_band_s_up
      this%width_band_p_up(i) = width_band_p_up
      this%width_band_d_up(i) = width_band_d_up

      this%width_band_s_dw(i) = width_band_s_dw
      this%width_band_p_dw(i) = width_band_p_dw
      this%width_band_d_dw(i) = width_band_d_dw
    
      close(funit)
    end do

  end subroutine read_potential   
end module potential_mod
