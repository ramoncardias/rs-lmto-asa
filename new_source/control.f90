!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Control
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
!> Module to handle data relative to basic control operations
!------------------------------------------------------------------------------


module control_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use precision_mod, only: rp
  use string_mod, only: sl
  implicit none

  private

  !> Module's main structure
  type, public :: control
    !> Recursion cutoff LL for d electrons
    !> 
    !> Recursion cutoff LL for d electrons
    integer :: lld

    !> Recursion cutoff LL for s-p electrons
    !> 
    !> Recursion cutoff LL for s-p electrons
    integer :: llsp

    !> Number of atoms in clust.
    !>
    !> Number of atoms in clust, where all first neighbors of the atoms under consideration will bee included.
    !>
    !> Use \ref nlim \f$= 0\f$ (zero) forbulk and surface.
    integer :: nlim
    integer :: npold

    !> Set calculation collinearity and relativistic type
    !>
    !> Set calculation collinearity and relativistic type
    !>
    !> Values:
    !>
    !> - \ref nsp \f$= 1\f$: Collinear scalar relativistic
    !> - \ref nsp \f$= 2\f$: Collinear fully relativistic (l.s) OP and noOP
    !> - \ref nsp \f$= 3\f$: Non-collinear scalar relativistic
    !> - \ref nsp \f$= 4\f$: Non-collinear fully relativistic (l.s)-only no OP
    integer :: nsp

    !> Type of LDOS s, p, d output
    !>
    !> Type of LDOS s, p, d output
    !> 
    !> Values:
    !>
    !> - \ref idos \f$= 0\f$: no LDOS output
    !> - \ref idos \f$= 1\f$: LDOS s, p, d only for the first type
    !> - \ref idos \f$= 2\f$: LDOS s, p, d for each type of atom
    integer :: idos
    
    !integer, dimension(:),allocatable :: ifc ! common_caro

    !> Set spin rotation
    !>
    !> If true, rotates in spin and real space so that local spin axis is along z-axis (beta version, use always false).
    !>
    !> **WARNIG:** Set true only for non-collinear calculations.
    logical :: lrot

    !> Includes orbital moment when calculating local spin axis
    !>
    !> Includes orbital moment when calculating local spin axis.
    !>
    !> **WARNIG:** Set true only for non-collinear calculations
    logical :: incorb

    !> Acceleration of rotation of spins
    !>
    !> Performs acceleration of rotation of spins
    !>
    !> Values:
    !>
    !> - \ref mext \f$= 0\f$: no
    !> - \ref mext \f$= 1\f$: linear extrapolation
    !> - \ref mext \f$= 2\f$: Broyden
    integer :: mext

    ! TODO
    logical :: do_asd

    ! TODO
    integer :: asd_atom

    logical :: do_comom ! common_cnstr

    !> Shift d-band levels.
    !>
    !> Shift d-band levels far up in energy for empty spheres (ES) (Emulates only sp-basis for ES) hoh
    !>
    !> Values:
    !>
    !> - \ref svac \f$= false\f$: Hamiltonian without hoh term
    !> - \ref svac \f$= true\f$: Hamiltonian with hoh term
    logical :: svac ! common_blb6

    !> Set type of calculation.
    !>
    !> Set type of calculation.
    !>
    !> Values:
    !>
    !> - \ref calctype \f$= B\f$: bulk
    !> - \ref calctype \f$= I\f$: impurity
    !> - \ref calctype \f$= S\f$: surface
    character :: calctype

    !> Description.
    !>
    !> Description.
    !> 
    !> Allowable values: 'lanczos', 'chebyshev'
    character(len=9) :: recur

    integer :: txc ! xcdata
    logical :: blockrec ! common_defs
    ! 0=over orbitals, 1=over atoms
    integer :: partype ! mpi_global_data
    logical :: do_cochg ! common_cnstr
    logical :: asd_jij ! common_asd 
    integer :: terminator ! common_asd 
    real(rp) :: conca,concb,ruban ! common_asr
    !real(rp), dimension(:,:), allocatable :: mom ! common_cnc
    integer :: njij ! common_c0
    integer :: nmdir
    character(len=sl) :: fname
  contains
    procedure :: build_from_file
    procedure :: restore_to_default
    procedure :: print_state
    procedure, private :: check_all
    final :: destructor
  end type control

  interface control
    procedure :: constructor
  end interface control

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] fname Input file with 'control' namelist
  !> @return type(control)
  !---------------------------------------------------------------------------
  function constructor(fname) result(obj)
    type(control) :: obj
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
    type(control) :: this
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
    class(control),intent(inout) :: this
    character(len=*), intent(in), optional :: fname

    ! variables associated with the object
    integer :: lld,llsp,nlim,npold,nsp,idos,mext,txc,partype,terminator,njij
    !integer, dimension(:), allocatable :: ifc                 
    real(rp) :: conca,concb,ruban
    logical :: lrot,incorb,do_asd,svac,blockrec,do_cochg,asd_jij,do_comom
    character(len=sl) :: calctype, recur

    ! variables associated with the reading processes
    integer :: iostatus, funit
    character(len=sl) :: fname_

    namelist /control/ recur, lld, llsp, nlim, npold, nsp,  &
    idos, lrot, incorb, do_asd, mext ,&
    svac, calctype, txc, blockrec, partype, do_cochg, asd_jij, terminator, conca,&
    concb, ruban, njij,do_comom !ifc, mom

    if (present(fname)) then
      fname_ = fname
      this%fname = fname
    else
      fname_ = this%fname
    endif

    ! Save previous values
    nlim = this%nlim
    npold = this%npold
    nsp = this%nsp
    calctype = this%calctype
    llsp = this%llsp
    lld = this%lld
    idos = this%idos
    lrot = this%lrot
    incorb = this%incorb
    mext = this%mext
    svac = this%svac
    txc = this%txc
    blockrec = this%blockrec
    partype = this%partype
    do_asd = this%do_asd
    do_cochg = this%do_cochg
    asd_jij = this%asd_jij
    terminator = this%terminator
    conca = this%conca
    concb = this%concb
    ruban = this%ruban
    njij = this%njij
    do_comom = this%do_comom
    recur = this%recur

    open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
    if(iostatus /= 0) then
      write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
      error stop
    endif

    read(funit,nml=control,iostat=iostatus)
    if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
      write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
      write(error_unit,'("iostatus = ",I0)') iostatus
    endif
    close(funit)

    ! Setting user values
    this%npold = npold
    this%llsp = llsp
    this%lld = lld
    this%idos = idos
    this%lrot = lrot
    this%incorb = incorb
    this%mext = mext
    this%svac = svac
    this%txc = txc
    this%blockrec = blockrec
    this%partype = partype
    this%do_asd = do_asd
    this%do_cochg = do_cochg
    this%asd_jij = asd_jij
    this%terminator = terminator
    this%conca = conca
    this%concb = concb
    this%ruban = ruban
    this%njij = njij
    this%do_comom = do_comom
    this%recur = recur
    ! end default

    ! Mandatory statements
    this%calctype = calctype
    this%nsp = nsp
    this%nlim = nlim
    if(this%nsp<=2)then
      this%nmdir=1
    else
      this%nmdir=3
    end if
    ! end reading the control file

    ! check input
    call this%check_all()
  end subroutine build_from_file

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Reset all members to default
  !---------------------------------------------------------------------------
  subroutine restore_to_default(this)
    implicit none
    class(control), intent(out) :: this

    this%calctype = ''
    this%npold = 9
    this%llsp = 16  
    this%lld = 16 
    this%idos = 0   
    this%lrot = .false.
    this%incorb = .false.
    this%mext = 0   
    this%svac = .false.
    this%txc = 1  
    this%blockrec = .false. 
    this%partype = 1         
    this%do_asd = .false.
    this%do_cochg = .false. 
    this%asd_jij = .false.
    this%terminator = 5         
    this%conca = 0.0d0
    this%concb = 0.0d0
    this%ruban = 0.0d0
    this%njij = 0   
    this%do_comom = .false.
    this%recur = 'lanczos'
    this%fname = ''
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
    class(control), intent(in) :: this

    integer,intent(in),optional :: unit
    character(len=*),intent(in),optional :: file
    integer :: newunit

    integer :: lld,llsp,nlim,npold,nsp,idos,mext,txc,partype,terminator,njij
    real(rp) :: conca,concb,ruban
    logical :: lrot,incorb,do_asd,svac,blockrec,do_cochg,asd_jij,do_comom
    character :: calctype

    namelist /control/ lld,llsp, nlim, npold, nsp,  &
    idos, lrot, incorb, do_asd, mext ,&
    svac, calctype, txc, blockrec, partype, do_cochg, asd_jij, terminator, conca,&
    concb, ruban, njij,do_comom

    lld = this%lld
    llsp = this%llsp
    nlim = this%nlim
    npold = this%npold
    nsp = this%nsp
    idos = this%idos
    lrot = this%lrot
    incorb = this%incorb
    do_asd = this%do_asd
    mext = this%mext 
    svac = this%svac
    calctype = this%calctype
    txc = this%txc
    blockrec = this%blockrec
    partype = this%partype
    do_cochg = this%do_cochg
    asd_jij = this%asd_jij
    terminator = this%terminator
    conca = this%conca
    concb = this%concb
    ruban = this%ruban
    njij = this%njij
    do_comom = this%do_comom

    if(present(unit) .and. present(file)) then
      write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
      error stop
    else if(present(unit)) then
        write(unit,nml=control)
    else if(present(file)) then
        open(unit=newunit,file=file)
        write(newunit,nml=control)
        close(newunit)
    else
        write(*,nml=control)
    endif
    close(newunit)
  end subroutine print_state

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Check if object is well fulfilled
  !
  !> @param[in] fname
  !> @return type(calculation)
  !---------------------------------------------------------------------------
  subroutine check_all(this)
    implicit none
    class(control) :: this
    if(this%calctype /= 'B' &
      .and. this%calctype /= 'S' &
      .and. this%calctype /= 'I' ) then 
      write(error_unit,'("[",A,":",I0,"]: lattice%calctype must be one of: ''B'', ''S'', ''I''")') __FILE__,__LINE__
      error stop
    end if
    if(this%recur /= 'lanczos' &
      .and. this%recur /= 'chebyshev' ) then
      write(error_unit,'("[",A,":",I0,"]: lattice%recur must be one of: ''lanczos'' or ''chebyshev''")') __FILE__,__LINE__
      error stop
    end if
  end subroutine check_all

end module control_mod


