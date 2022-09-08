!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Self
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
!> Module to handle self-consistent process
!------------------------------------------------------------------------------


module self_mod

  use control_mod
  use lattice_mod
  use math_mod, only : ang2au
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use precision_mod, only: rp
  implicit none

  private

  !> Module's main structure
  type, public :: self
    !TODO: check description
    !> If true treats all atoms as inequivalents. Default: true.
    !>
    !> If true treats all atoms as inequivalents.
    !>
    !> Default: true for all kind of calculations. If false, the user may provide the final lines of old format self file in a different file.
    logical :: all_inequivalent
    
    ! Calculation type dependents

    !TODO: check description
    !> Specify if calculation is performed with fixed fermi level.
    !>
    !> Specify if calculation is performed with fixed fermi level.
    !>
    !> Default: false for 'bulk' calculation and
    !> true for 'surface' and 'impurity' calculations.
    logical :: fix_fermi

    !TODO: lack description
    !> Default: 8
    !>
    !> Default: 8
    integer :: init

    !TODO: check description
    !> Specify the number of ldos channels in calculation. Default: 6000.
    !>
    !> Specify the number of Local Density of States (ldos) channels used in calculation.
    !>
    !> Default: 6000.
    !>
    !> Suggestitive values: \n
    !> - 6000, if calctype = 'B'; \n
    !> - 3000, if calctype = 'I'.
    integer :: channels_ldos


    !TODO: check description
    !> Fermi energy (Ry). Default: 0.05 for bulk calculation.
    !>
    !> Fermi energy (Ry)
    !> 
    !> Default: 0.05 for bulk calculation. For 'S' and 'I' calculations read from previous 'bulk' calculation.
    real(rp) :: fermi, chebfermi

    !> Minimum value of energy in calculation (Ry). Default: -5.5.
    !>
    !> Minimum value of energy in calculation (Ry).
    !>
    !> Default: -5.5.
    !> 
    !> Suggestitive values: \n
    !> - -1.5, if calctype = 'B' or calctype = 'I'.
    real(rp) :: energy_min

    !> Maximum value of energy in calculation (Ry). Default: 5.5.
    !>
    !> Maximum value of energy in calculation (Ry).
    !>
    !> Default: 5.5.
    !> 
    !> Suggestitive values: \n
    !> - 0.5, if calctype = 'B' or calctype = 'I'.
    real(rp) :: energy_max

    !> Energy mesh
    real(rp) :: edel              
    real(rp), dimension(:), allocatable :: ene
    integer :: enpt
    ! Wigner Seitz Radius

    !> Use same value of @ref ws for all atoms. Default: true.
    !> 
    !> Use same value of @ref ws for all atoms. If true, @ref ws's size is one slot of memory, else it is @ref lattice.nrec
    !> 
    !> Default: true.
    logical :: ws_all

    !> Wigner Seitz Radius
    !> 
    !> Wigner Seitz Radius. If @ref ws_all is true, @ref ws's size is one slot of memory, else it is @ref lattice.nrec
    real(rp), dimension(:), allocatable :: ws
    
    ! Mixing parameters

    !> Use same value of @ref mix for all atoms. Default: true.
    !> 
    !> Use same value of @ref mix for all atoms. If true, @ref mix's size is one slot of memory, else it is @ref lattice.nrec
    !> 
    !> Default: true.
    logical :: mix_all

    !> Mixture occupation in self-consistent calculation. Default: 0.01.
    !> 
    !> Mixture occupation in self-consistent calculation.
    !> 
    !> Default: 0.01.
    real(rp), dimension(:), allocatable :: mix

    ! Magnetic mixing parameters

    !> If true enable spin-(up/down) mixing. Default: false.
    !>
    !> If true enable spin-(up/down) mixing, else disable.
    !>
    !> Default: false.
    logical :: magnetic_mixing

    !> Use same value of @ref mixmag for all atoms. Default: true.
    !> 
    !> Use same value of @ref mixmag for all atoms. If true, @ref mixmag's size is one slot of memory, else it is @ref lattice.nrec
    !> 
    !> Default: true.
    logical :: mixmag_all

    !> Spin-(up/down) occupation in self-consistent calculation. Default: 0.05.
    !> 
    !> Spin-(up/down) occupation in self-consistent calculation.
    !> 
    !> Default: 0.05.
    real(rp), dimension(:), allocatable :: mixmag
    
    ! Convergence parameters

    !> Number of steps in calculation. Default: 1.
    !>
    !> Number of steps in self-consistent calculation.
    !>
    !> Default: 1.
    integer :: nstep

    !> Convergency precision threshold. Default: 0.5d-10.
    !>
    !> Specify the precision required in charge density to stop the self-consistent process.
    !>
    !> Default: 0.5d-10.
    real(rp) :: conv_thr

    ! Constrains variables

    !> Freezes magnetic moment direction during nonlinear calculation. Default: false.
    !>
    !> Freezes magnetic moment direction during nonlinear calculation. Do not set true in collinear calculation.
    !>
    !> Default: false.
    logical :: freeze

    !> If true rigid band calculations will performed for each atom. Default: false.
    !>
    !> If true rigid band calculations will performed for each atom. Using this option you should specify the variable @ref rb.
    !>
    !> Default: false.
    logical :: rigid_band

    !> Specify the number of rigid band calculation used for each atom.
    !>
    !> Specify the number of rigid band calculation used for each atom. The number of values provided should be the same as @ref lattice.nrec.
    integer, dimension(:), allocatable :: rb

    !> If true, performs a calculation with moment orbital polarization. Default: false.
    !>
    !> If true, performs a calculation with moment orbital polarization.
    !>
    !> Default: false.
    logical :: orbital_polarization

    !> Pointer to system's lattice.
    !>
    !> Pointer to system's lattice.
    type(lattice), pointer :: lattice

  contains
    procedure :: build_from_file
    procedure :: restore_to_default
    procedure :: print_state_formatted
    procedure :: print_state
    procedure :: print_state_full
    final :: destructor
  end type self

  interface self
    procedure :: constructor
  end interface self

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] fname Namelist file
  !> @param[in] lattice_obj Pointer to system's lattice
  !> @return type(self)
  !---------------------------------------------------------------------------
  function constructor(fname,lattice_obj) result(obj)
    type(self) :: obj
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
    type(self) :: this
    if(allocated(this%ws)) deallocate(this%ws)
    if(allocated(this%mix)) deallocate(this%mix)
    if(allocated(this%mixmag)) deallocate(this%mixmag)
    if(allocated(this%rb)) deallocate(this%rb)
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
    implicit none
    class(self), intent(inout) :: this
    character(len=*), intent(in) :: fname

    ! variables associated with the object
    logical :: all_inequivalent, fix_fermi, ws_all, mix_all, &
    magnetic_mixing, mixmag_all, freeze, rigid_band, orbital_polarization
    integer :: init, channels_ldos, nstep 
    integer, dimension(:), allocatable :: rb
    real(rp) :: fermi, energy_min, energy_max, conv_thr
    real(rp), dimension(:), allocatable :: mix, ws, mixmag

    ! variables associated with the reading processes
    integer :: iostatus, funit, i

    namelist /self/ ws_all, fix_fermi, all_inequivalent, &
    mix_all, magnetic_mixing, mixmag_all, freeze, orbital_polarization, &
    rigid_band, rb, channels_ldos, nstep, init, energy_min, & 
    energy_max, fermi, conv_thr, ws, mix, mixmag

    ! Save previous values
    all_inequivalent = this%all_inequivalent
    channels_ldos = this%channels_ldos
    energy_min = this%energy_min
    energy_max = this%energy_max
    fermi = this%fermi
    fix_fermi = this%fix_fermi
    ws_all = this%ws_all
    mix_all = this%mix_all
    magnetic_mixing = this%magnetic_mixing
    mixmag_all = this%mixmag_all
    freeze = this%freeze
    rigid_band = this%rigid_band
    orbital_polarization = this%orbital_polarization
    init = this%init
    nstep = this%nstep
    conv_thr = this%conv_thr

    call move_alloc(this%mix, mix)
    call move_alloc(this%ws, ws)
    call move_alloc(this%mixmag, mixmag)
    call move_alloc(this%rb, rb)

    ! Check if size is right
    if(size(mix) .ne. this%lattice%nrec) then
      write(*,*) '[',__FILE__,':',__LINE__,']: resizing array "mix"' 
      write(*,*) 'from ',size(mix),' to ',this%lattice%nrec
      write(*,*) 'If this variable is not in the namelist, it may cause errors'
      deallocate(mix)
      allocate(mix(this%lattice%nrec))
    endif
    if(size(mixmag) .ne. this%lattice%nrec) then
      write(*,*) '[',__FILE__,':',__LINE__,']: resizing array "mixmag"' 
      write(*,*) 'from ',size(mixmag),' to ',this%lattice%nrec
      write(*,*) 'If this variable is not in the namelist, it may cause errors'
      deallocate(mixmag)
      allocate(mixmag(this%lattice%nrec))
    endif
    if(size(ws) .ne. this%lattice%nrec) then
      write(*,*) '[',__FILE__,':',__LINE__,']: resizing array "ws"' 
      write(*,*) 'from ',size(ws),' to ',this%lattice%nrec
      write(*,*) 'If this variable is not in the namelist, it may cause errors'
      deallocate(ws)
      allocate(ws(this%lattice%nrec))
    endif
    if(size(rb) .ne. this%lattice%nrec) then
      write(*,*) '[',__FILE__,':',__LINE__,']: resizing array "rb"' 
      write(*,*) 'from ',size(rb),' to ',this%lattice%nrec
      write(*,*) 'If this variable is not in the namelist, it may cause errors'
      deallocate(rb)
      allocate(rb(this%lattice%nrec))
    endif

    ! Reading
    open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
    if(iostatus /= 0) then
      write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
      error stop
    endif
    
    read(funit,nml=self,iostat=iostatus)
    if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
      write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
      write(error_unit,'("iostatus = ",I0)') iostatus
    endif
    close(funit)

    ! Setting user values

    ! Control variables
    this%all_inequivalent = all_inequivalent
    this%channels_ldos = channels_ldos
    this%energy_min = energy_min
    this%energy_max = energy_max
    this%fermi = fermi
    this%fix_fermi = fix_fermi

    ! Wigner Seitz Radius
    this%ws_all = ws_all
    if( ws_all ) then
      allocate(this%ws(1))
      this%ws(1) = ws(1)
    else
      call move_alloc(ws, this%ws)
    endif
    
    ! Mixing parameters
    this%mix_all = mix_all
    if( mix_all ) then
      allocate(this%mix(1))
      this%mix(1) = mix(1)
    else
      call move_alloc(mix, this%mix)
    endif
    
    ! Magnetic mixing parameters
    this%magnetic_mixing = magnetic_mixing
    this%mixmag_all = mixmag_all
    if( magnetic_mixing ) then
      if( mixmag_all ) then
        allocate(this%mixmag(1))
        this%mixmag(1) = mixmag(1)
      else
        call move_alloc(mixmag, this%mixmag)
      endif
    else
      allocate(this%mixmag(0))
    endif

    ! Convergence parameters
    this%conv_thr = conv_thr
    this%nstep = nstep

    ! Constrains variables
    ! static magnetic momentum (non-linar calculation)
    this%freeze = freeze
    this%rigid_band = rigid_band

    ! number of loops to use rigid bands per atom
    call move_alloc(rb,this%rb)

    ! other variables
    this%orbital_polarization = orbital_polarization
    this%init = init

    ! Solve the energy mesh
    allocate(this%ene(this%channels_ldos+10))
    this%edel = (this%energy_max-this%energy_min)/this%channels_ldos
    this%enpt = nint((this%fermi-this%energy_min)/this%edel)
    this%edel = (this%fermi-this%energy_min)/nint((this%fermi-this%energy_min)/this%edel)

    do i=0,this%channels_ldos+9
      this%ene(i+1) = this%energy_min + this%edel*i
    end do
 end subroutine build_from_file

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Reset all members to default
  !---------------------------------------------------------------------------
 subroutine restore_to_default(this)
    implicit none
    class(self), intent(inout):: this

    ! Control variables
    ! if false force to read the original self file
    this%all_inequivalent = .true.
    ! variables according to calculation type
    select case (this%lattice%control%calctype)
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
          this%init = 8
    end select

    ! Wigner Seitz Radius
    this%ws_all = .true.
    allocate(this%ws(this%lattice%nrec))
    this%ws(:) = this%lattice%wav * ang2au
    
    ! Mixing parameters
    this%mix_all = .true.
    allocate(this%mix(this%lattice%nrec))
    this%mix(:) = 0.01
    
    ! Magnetic mixing parameters
    this%magnetic_mixing = .false.
    this%mixmag_all = .true.
    allocate(this%mixmag(this%lattice%nrec))
    this%mixmag(:) = 0.05

    ! Convergence parameters
    this%conv_thr = 0.5d-10
    this%nstep = 1

    ! Constrains variables
    ! static magnetic momentum (non-linar calculation)
    this%freeze = .false.
    this%rigid_band = .false.
    ! number of loops to use rigid bands per atom
    allocate(this%rb(this%lattice%nrec))
    this%rb(:) = 2

    this%orbital_polarization = .false.
  end subroutine restore_to_default

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Print class members values formatted 
  !>
  !> Print class members values formatted 
  !---------------------------------------------------------------------------
  subroutine print_state_formatted(this)
    implicit none
    class(self), intent(in) :: this

    print*, 'Printing SELF object'
    print*, ''
    print*, '[Control Variables]'
    print*, 'all_inequivalent ', this%all_inequivalent
    print*, ''
    print*, '[Calculation Type Dependent Variables]'
    print*, 'channels_ldos ', this%channels_ldos
    print*, 'energy_min    ', this%energy_min
    print*, 'energy_max    ', this%energy_max
    print*, 'fermi         ', this%fermi
    print*, 'fix_fermi     ', this%fix_fermi
    print*, ''
    print*, '[Wigner Seitz Radius]'
    print*, 'ws_all ', this%ws_all
    print*, 'ws     ', this%ws
    print*, ''
    print*, '[Mixing Parameters]'
    print*, 'mix_all ', this%mix_all
    print*, 'mix     ', this%mix
    print*, ''
    print*, '[Magnetic Mixing Parameters]'
    print*, 'magnetic_mixing ', this%magnetic_mixing
    print*, 'mixmag_all      ', this%mixmag_all
    print*, 'mixmag          ', this%mixmag
    print*, ''
    print*, '[Convergence Parameters]'
    print*, 'conv_thr ', this%conv_thr
    print*, 'nstep    ', this%nstep
    print*, ''
    print*, '[Constrain Variables]'
    print*, 'freeze     ', this%freeze
    print*, 'rigid_band ', this%rigid_band
    print*, 'rb         ', this%rb
    print*, ''
    print*, '[Other Variables]'
    print*, 'orbital_polarization ', this%orbital_polarization
    print*, 'init                 ', this%init

  end subroutine print_state_formatted


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Print class members values in namelist format 
  !>
  !> Print class members values in namelist format 
  !---------------------------------------------------------------------------
  subroutine print_state_full(this)
    implicit none
    class(self), intent(in) :: this

    real(rp) :: fermi, energy_min, energy_max, conv_thr
    logical :: ws_all, rigid_band, orbital_polarization, mixmag_all, mix_all, magnetic_mixing, freeze, fix_fermi, all_inequivalent
    integer :: nstep, init, channels_ldos
    real(rp), dimension(:), allocatable :: ws, mixmag, mix
    integer, dimension(:), allocatable :: rb

    namelist /self/ fermi, energy_min, energy_max, &
    conv_thr, ws_all, rigid_band, orbital_polarization,&
    mixmag_all, mix_all, magnetic_mixing, freeze, fix_fermi, &
    all_inequivalent, nstep, init, channels_ldos, ws, mixmag, &
    mix, rb

    ! scalar

    fermi = this%fermi
    energy_min = this%energy_min
    energy_max = this%energy_max
    conv_thr = this%conv_thr
    ws_all = this%ws_all
    rigid_band = this%rigid_band
    orbital_polarization = this%orbital_polarization
    mixmag_all = this%mixmag_all
    mix_all = this%mix_all
    magnetic_mixing = this%magnetic_mixing
    freeze = this%freeze
    fix_fermi = this%fix_fermi
    all_inequivalent = this%all_inequivalent
    nstep = this%nstep
    init = this%init
    channels_ldos = this%channels_ldos

    ! 1d allocatable
    
    if(allocated(this%ws)) then
      allocate(ws,mold=this%ws)
      ws = this%ws
    else
      allocate(ws(0))
    endif
    if(allocated(this%mixmag)) then
      allocate(mixmag,mold=this%mixmag)
      mixmag = this%mixmag
    else
      allocate(mixmag(0))
    endif
    if(allocated(this%mix)) then
      allocate(mix,mold=this%mix)
      mix = this%mix
    else
      allocate(mix(0))
    endif
    if(allocated(this%rb)) then
      allocate(rb,mold=this%rb)
      rb = this%rb
    else
      allocate(rb(0))
    endif

    write(*,nml=self)

  end subroutine print_state_full

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Print input class members values in namelist format 
  !>
  !> Print input class members values in namelist format 
  !---------------------------------------------------------------------------
  subroutine print_state(this)
    implicit none
    class(self), intent(in) :: this

    real(rp) :: fermi, energy_min, energy_max, conv_thr
    logical :: ws_all, rigid_band, orbital_polarization, mixmag_all, mix_all, magnetic_mixing, freeze, fix_fermi, all_inequivalent
    integer :: nstep, init, channels_ldos
    real(rp), dimension(:), allocatable :: ws, mixmag, mix
    integer, dimension(:), allocatable :: rb

    namelist /self/ ws_all, fix_fermi, all_inequivalent, &
    mix_all, magnetic_mixing, mixmag_all, freeze, orbital_polarization, &
    rigid_band, rb, channels_ldos, nstep, init, energy_min, & 
    energy_max, fermi, conv_thr, ws, mix, mixmag


    ! scalar

    fermi = this%fermi
    energy_min = this%energy_min
    energy_max = this%energy_max
    conv_thr = this%conv_thr
    ws_all = this%ws_all
    rigid_band = this%rigid_band
    orbital_polarization = this%orbital_polarization
    mixmag_all = this%mixmag_all
    mix_all = this%mix_all
    magnetic_mixing = this%magnetic_mixing
    freeze = this%freeze
    fix_fermi = this%fix_fermi
    all_inequivalent = this%all_inequivalent
    nstep = this%nstep
    init = this%init
    channels_ldos = this%channels_ldos

    ! 1d allocatable
    
    if(allocated(this%ws)) then
      allocate(ws,mold=this%ws)
      ws = this%ws
    else
      allocate(ws(0))
    endif
    if(allocated(this%mixmag)) then
      allocate(mixmag,mold=this%mixmag)
      mixmag = this%mixmag
    else
      allocate(mixmag(0))
    endif
    if(allocated(this%mix)) then
      allocate(mix,mold=this%mix)
      mix = this%mix
    else
      allocate(mix(0))
    endif
    if(allocated(this%rb)) then
      allocate(rb,mold=this%rb)
      rb = this%rb
    else
      allocate(rb(0))
    endif

    write(*,nml=self)

  end subroutine print_state
end module self_mod


