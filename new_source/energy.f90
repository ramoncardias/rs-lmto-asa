!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Energy
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas f. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle self-consistent process
!------------------------------------------------------------------------------

module energy_mod

  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use lattice_mod
  use control_mod
  use logger_mod, only: g_logger
  use math_mod
  use precision_mod
  use string_mod
  implicit none

  private

  !> Module's main structure
  type, public :: energy
    !> Lattice
    class(lattice), pointer :: lattice

    !> Specify the number of Local Density of States (ldos) energy points used in the calculation.
    integer ::  channels_ldos 
    !> Fermi energy (Ry).
    real(rp) :: fermi, chebfermi
    !> Lower energy limit
    real(rp) :: energy_min
    !> Upper energy limit
    real(rp) :: energy_max 
    !> Energy mesh
    real(rp) :: edel
    real(rp), dimension(:), allocatable :: ene
    integer :: ik1, nv1
    integer :: enpt
    !> Logical fix fermi level
    logical :: fix_fermi
  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: build_from_file
    procedure :: restore_to_default
    procedure :: e_mesh  
  end type energy 

  ! Constructor
  interface energy
    procedure :: constructor
  end interface energy 

contains
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] lattice_obj Pointer to system's lattice
  !> @return type(energy)
  !---------------------------------------------------------------------------
  function constructor(lattice_obj) result(obj)
    type(energy) :: obj
    type(lattice), target, intent(in) :: lattice_obj
  
    obj%lattice => lattice_obj
   
    call obj%restore_to_default()
    call obj%build_from_file()
  end function


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine destructor(this)
    type(energy) :: this
    if(allocated(this%ene)) deallocate(this%ene)
  end subroutine destructor


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Read parameters from input file
  !---------------------------------------------------------------------------
  subroutine build_from_file(this)
    class(energy), intent(inout) :: this

    ! Namelist variables
    logical :: fix_fermi
    integer :: channels_ldos
    real(rp) :: fermi, energy_min, energy_max
  
    ! Reading process variables
    integer :: iostatus, funit, i
  
    namelist /energy/ fix_fermi, channels_ldos, fermi, energy_min, energy_max 
        

    ! Save previous values
    fix_fermi = this%fix_fermi
    channels_ldos = this%channels_ldos
    fermi = this%fermi
    energy_min = this%energy_min
    energy_max = this%energy_max

    ! Reading
    open(newunit=funit,file=this%lattice%control%fname,action='read',iostat=iostatus,status='old')
    if(iostatus /= 0) then
      call g_logger%fatal('file '//trim(this%lattice%control%fname)//' not found',__FILE__,__LINE__)
    endif
   
    read(funit,nml=energy,iostat=iostatus)
    if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
      call g_logger%error('Error while reading namelist',__FILE__,__LINE__)
      call g_logger%error('iostatus = '//int2str(iostatus),__FILE__,__LINE__)
    endif
    close(funit)

    this%fix_fermi = fix_fermi
    this%channels_ldos = channels_ldos
    this%fermi = fermi
    this%energy_min = energy_min
    this%energy_max = energy_max
  end subroutine build_from_file

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Reset all members to default
  !---------------------------------------------------------------------------
  subroutine restore_to_default(this)
    class(energy), intent(inout) :: this
 
    select case(this%lattice%control%calctype)
      case( 'B' ) 
          this%channels_ldos = 6000
          this%energy_min = -5.5
          this%energy_max = 5.5
          this%fermi = -0.05
          this%fix_fermi = .false.
      case ( 'I' )
          this%channels_ldos = 3000
          this%energy_min = -1.5
          this%energy_max = 0.5
          ! this%fermi = readed from bulk calculation
          this%fix_fermi = .true.
      case ( 'S' )
          this%energy_min = -1.5
          this%energy_max = 0.5
          ! this%fermi = readed from bulk calculation
          this%fix_fermi = .true.
          this%channels_ldos = 6000
    end select
  end subroutine restore_to_default

  subroutine e_mesh(this)
    class(energy), intent(inout) :: this
    integer :: i
    real(rp) :: prevfermi

    prevfermi = this%fermi
 
    call this%build_from_file()

    this%fermi = prevfermi
    !Make sure nv1 is odd
    if(mod(this%channels_ldos,2)==0)then
      this%nv1 = this%channels_ldos + 1
    else
      this%nv1 = this%channels_ldos
      this%channels_ldos = this%channels_ldos - 1
    end if
    this%ik1 = this%nv1

    if(allocated(this%ene)) deallocate(this%ene)
    allocate(this%ene(this%channels_ldos+10))

    this%edel = (this%energy_max-this%energy_min)/this%channels_ldos
    this%enpt = nint((this%fermi-this%energy_min)/this%edel)
    this%edel = (this%fermi-this%energy_min)/nint((this%fermi-this%energy_min)/this%edel)

    do i=0,this%channels_ldos+9
      this%ene(i+1) = this%energy_min + this%edel*i
    end do
  end subroutine e_mesh
end module energy_mod
