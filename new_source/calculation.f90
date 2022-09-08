!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Calculation
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
!> Module to handle pre-processing, processing and post-processing calculations
!------------------------------------------------------------------------------


module calculation_mod

  use control_mod
  use self_mod
  use lattice_mod
  use charge_mod
  use symbolic_atom_mod
  use hamiltonian_mod
  use recursion_mod
  use green_mod
  use density_of_states_mod
  use bands_mod
  use precision_mod, only: rp
  use string_mod, only: sl
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  implicit none

  private

  type, public :: calculation
    !> Pre-processing. Options are:
    !> 'none' (default)
    !> 'bravais' : Builds the bulk clust
    !> 'buildsurf' : Builds the surface clust
    !> 'newclubulk' : Builds the imputiry clust from the bluk clust
    !> 'newclusurf' : Builds the impurity clust from the surface clust
    character(len=10) :: pre_processing

    !> Processing. Options are
    !> 'none' (default)
    character(len=10) :: processing

    !> Post-processing. Options are
    !> 'none' (default)
    character(len=10) :: post_processing
    
    !> Controller for preprocessing verbosity.
    !>
    !> Controller for preprocessing verbosity. If true, call the subroutines:
    !> @ref control.print_state, @ref lattice.print_state,
    !> @ref self.print_state and @ref charge.print_state after pre-processing calls.
    logical :: verbose

    !> name list input file
    character(len=sl) :: fname
  contains
    procedure :: build_from_file
    procedure :: restore_to_default
    procedure, private :: pre_processing_bravais
    procedure, private :: pre_processing_buildsurf
    procedure, private :: pre_processing_newclubulk
    procedure, private :: pre_processing_newclusurf
    procedure :: process
    final :: destructor
  end type calculation

  interface calculation
    procedure :: constructor
  end interface calculation

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] fname Namelist file
  !> @return type(calculation)
  !---------------------------------------------------------------------------
  function constructor(fname) result(obj)
    type(calculation) :: obj
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
    type(calculation) :: this
  end subroutine destructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Read parameters from input file
  !
  !> @param[in] fname Namelist file with '&calculation' table
  !---------------------------------------------------------------------------
  subroutine build_from_file(this,fname)
    class(calculation),intent(inout) :: this
    character(len=*),intent(in) :: fname
    ! Namelist variables
    character(len=10) :: pre_processing
    character(len=10) :: processing
    character(len=10) :: post_processing
    logical :: verbose
      ! variables associated with the reading processes
    integer :: iostatus, funit

    ! Namelist
    namelist /calculation/ verbose, pre_processing, processing, post_processing

    verbose = this%verbose
    pre_processing = this%pre_processing
    processing = this%processing
    post_processing = this%post_processing

    open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
    if(iostatus /= 0) then
      write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
      error stop
    endif

    read(funit,nml=calculation,iostat=iostatus)
    if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
      write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
      write(error_unit,'("iostatus = ",I0)') iostatus
    endif

    ! Pre-processing
    call check_pre_processing(trim(pre_processing))
    ! Processing
    call check_processing(trim(processing))
    ! Post-processing
    call check_post_processing(trim(post_processing))

    this%verbose = verbose
    this%fname = fname
    this%pre_processing = pre_processing
    this%processing = processing
    this%post_processing = post_processing

    close(funit)
  end subroutine build_from_file

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Handle for general process
  !---------------------------------------------------------------------------
  subroutine process(this)
    class(calculation), intent(in) :: this

    ! Pre-processing
    select case(this%pre_processing)
    case('bravais')
      call this%pre_processing_bravais()
    case('buildsurf')
      call this%pre_processing_buildsurf()
    case('newclubulk')
      call this%pre_processing_newclubulk()
    case('newclusurf')
      call this%pre_processing_newclusurf()
    end select

  end subroutine

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Pre-process for new clust bulk calculation
  !---------------------------------------------------------------------------
  subroutine pre_processing_newclubulk(this)
    class(calculation),intent(in) :: this

    type(control) :: control_obj
    type(lattice) :: lattice_obj
    type(self) :: self_obj
    type(charge) :: charge_obj
    type(symbolic_atom), dimension(:), allocatable :: symbolic_atoms_obj
    type(hamiltonian) :: hamiltonian_obj
    type(recursion) :: recursion_obj
    real(rp) :: start, finish
    integer :: i

    call control_obj%restore_to_default
    call control_obj%build_from_file(this%fname)
    !control_obj = control_constructor(this%fname)
    lattice_obj = lattice(this%fname,control_obj)
    call lattice_obj%build_data()
    call lattice_obj%bravais()
    call lattice_obj%newclu()
    call lattice_obj%structb()
    call lattice_obj%atomlist()
    self_obj = self(this%fname,lattice_obj)
    charge_obj = charge(this%fname,lattice_obj)
    call charge_obj%impmad
    symbolic_atoms_obj = array_of_symbolic_atoms(this%fname,lattice_obj)
    hamiltonian_obj = hamiltonian(symbolic_atoms_obj,charge_obj)
    call hamiltonian_obj%build_bulkham
    call hamiltonian_obj%build_locham
    recursion_obj = recursion(hamiltonian_obj,self_obj)
    call cpu_time(start)
    call recursion_obj%recur
    call cpu_time(finish)
    print '("Recursion time = ",f10.3," seconds.")',(finish-start)/32

    if(this%verbose) then
      call control_obj%print_state()
      call lattice_obj%print_state()
      call self_obj%print_state()
      call charge_obj%print_state()
      call save_state(symbolic_atoms_obj)
    endif

end subroutine pre_processing_newclubulk

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Pre-process for new clust surface calculation
  !---------------------------------------------------------------------------
  subroutine pre_processing_newclusurf(this)
    class(calculation),intent(in) :: this

    type(control) :: control_obj
    type(lattice) :: lattice_obj
    type(self) :: self_obj
    type(charge) :: charge_obj
    type(symbolic_atom), dimension(:), allocatable :: symbolic_atoms_obj
    type(hamiltonian) :: hamiltonian_obj
    type(recursion) :: recursion_obj
    real(rp) :: start, finish
    integer :: i

    call control_obj%restore_to_default()
    call control_obj%build_from_file(this%fname)
    !control_obj = control(this%fname)
    lattice_obj = lattice(this%fname,control_obj)
    call lattice_obj%build_data()
    call lattice_obj%bravais()
    call lattice_obj%build_clusup()
    call lattice_obj%build_surf()
    call lattice_obj%newclu()
    call lattice_obj%structb()
    call lattice_obj%atomlist
    self_obj = self(this%fname,lattice_obj)
    charge_obj = charge(this%fname,lattice_obj)
    call charge_obj%impmad
    symbolic_atoms_obj = array_of_symbolic_atoms(this%fname,lattice_obj)
    hamiltonian_obj = hamiltonian(symbolic_atoms_obj,charge_obj)
    call hamiltonian_obj%build_bulkham
    call hamiltonian_obj%build_locham
    recursion_obj = recursion(hamiltonian_obj,self_obj)
    call cpu_time(start)
    call recursion_obj%recur
    call cpu_time(finish)
    print '("Recursion time = ",f10.3," seconds.")',(finish-start)/32

    if(this%verbose) then
      call control_obj%print_state()
      call lattice_obj%print_state()
      call self_obj%print_state()
      call charge_obj%print_state()
      call save_state(symbolic_atoms_obj)
    endif

end subroutine pre_processing_newclusurf

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Pre-process for build surface calculation
  !---------------------------------------------------------------------------
  subroutine pre_processing_buildsurf(this)
    class(calculation),intent(in) :: this

    type(control) :: control_obj
    type(lattice) :: lattice_obj
    type(self) :: self_obj
    type(charge) :: charge_obj
    type(symbolic_atom), dimension(:), allocatable :: symbolic_atoms_obj
    type(hamiltonian) :: hamiltonian_obj
    type(recursion) :: recursion_obj
    real(rp) :: start, finish
    integer :: i

    call control_obj%restore_to_default
    call control_obj%build_from_file(this%fname)
    !control_obj = control(this%fname)
    lattice_obj = lattice(this%fname,control_obj)
    call lattice_obj%build_data()
    call lattice_obj%bravais()
    call lattice_obj%build_clusup()
    call lattice_obj%build_surf()
    call lattice_obj%structb()
    call lattice_obj%atomlist()
    self_obj = self(this%fname,lattice_obj)
    charge_obj = charge(this%fname,lattice_obj)
    call charge_obj%build_alelay
    call charge_obj%surfmat
    symbolic_atoms_obj = array_of_symbolic_atoms(this%fname,lattice_obj)
    hamiltonian_obj = hamiltonian(symbolic_atoms_obj,charge_obj)
    call hamiltonian_obj%build_bulkham
    recursion_obj = recursion(hamiltonian_obj,self_obj)
    call cpu_time(start)
    call recursion_obj%recur
    call cpu_time(finish)
    print '("Recursion time = ",f10.3," seconds.")',(finish-start)/32

    if(this%verbose) then
      call control_obj%print_state()
      call lattice_obj%print_state()
      call self_obj%print_state()
      call charge_obj%print_state()
      call save_state(symbolic_atoms_obj)
    endif

end subroutine pre_processing_buildsurf

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Pre-process for bravais calculation
  !---------------------------------------------------------------------------
  subroutine pre_processing_bravais(this)
    class(calculation),intent(in) :: this

    type(control) :: control_obj
    type(lattice) :: lattice_obj
    type(self) :: self_obj
    type(charge) :: charge_obj
    type(symbolic_atom), dimension(:), allocatable :: symbolic_atoms_obj
    type(hamiltonian) :: hamiltonian_obj
    type(recursion) :: recursion_obj
    type(green) :: green_obj
    type(dos) :: dos_obj
    type(bands) :: bands_obj
    real(rp) :: start, finish
    integer :: i

 
    call control_obj%restore_to_default
    call control_obj%build_from_file(this%fname)
    !control_obj = control(this%fname)
    lattice_obj = lattice(this%fname,control_obj)
    call cpu_time(start)
    call lattice_obj%build_data()
    call lattice_obj%bravais()
    call lattice_obj%structb()
    call cpu_time(finish)
    print '("Pre-processing time = ",f10.3," seconds.")',(finish-start)!/32

    call lattice_obj%atomlist()
    self_obj = self(this%fname,lattice_obj)
    charge_obj = charge(this%fname,lattice_obj)
    call charge_obj%bulkmat
    symbolic_atoms_obj = array_of_symbolic_atoms(this%fname,lattice_obj)
    hamiltonian_obj = hamiltonian(symbolic_atoms_obj,charge_obj)
    call hamiltonian_obj%build_bulkham
    !call hamiltonian_obj%build_lsham
    recursion_obj = recursion(hamiltonian_obj,self_obj)
    !call cpu_time(start)
    call recursion_obj%recur()
    !call cpu_time(finish)
    !print '("Recursion time = ",f10.3," seconds.")',(finish-start)!/32

    dos_obj = dos(recursion_obj,self_obj)
    !Chebyshev test
    call cpu_time(start)
!    call recursion_obj%chebyshev_recur()
!    call recursion_obj%chebyshev_recur_full()
    call cpu_time(finish)
    print '("Recursion Chebyshev time = ",f10.3," seconds.")',(finish-start)!/32

!    call dos_obj%chebyshev_dos()
!    call dos_obj%chebyshev_dos_full()

    green_obj = green(dos_obj)
    call green_obj%sgreen()
    bands_obj = bands(green_obj)
    call bands_obj%calculate_fermi()
    call bands_obj%calculate_moments()
    call bands_obj%calculate_moments_chebgauss()

    if(this%verbose) then
      call control_obj%print_state()
      call lattice_obj%print_state()
      call self_obj%print_state()
      call charge_obj%print_state()
      call save_state(symbolic_atoms_obj)
      call print_state(symbolic_atoms_obj)
    endif

end subroutine pre_processing_bravais

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Reset all members to default ('none') value
  !---------------------------------------------------------------------------
  subroutine restore_to_default(this)
    class(calculation) :: this
    this%pre_processing      = 'none'
    this%processing          = 'none'
    this%post_processing     = 'none'
  end subroutine restore_to_default


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Check availability for post-processing
  !
  !> @param post_processing Type of processing. Allowed values: 'none'
  !---------------------------------------------------------------------------
  subroutine check_post_processing(post_processing)
    character(len=*),intent(in) :: post_processing
    if(post_processing /= 'none') then
      write(error_unit,*) '[calculation.check_post_processing]: &
      calculation%post_processing must be one of: ''none'''
      error stop
    end if
  end subroutine check_post_processing


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Check availability for pre-processing
  !
  !> @param[in] pre_processing Type of pre-processing. Allowed values:
  !> 'bravais', 'buildsurf', 'newclubulk', 'newclusurf', 'none'
  !---------------------------------------------------------------------------
  subroutine check_pre_processing(pre_processing)
    character(len=*),intent(in) :: pre_processing
    if(pre_processing /= 'none' &
     .and. pre_processing /= 'bravais' &
     .and. pre_processing /= 'buildsurf' &
     .and. pre_processing /= 'newclubulk'&
     .and. pre_processing /= 'newclusurf' ) then
      write(error_unit,*) '[calculation.check_pre_processing]: &
       calculation%pre_processing must be one of: ''none'', ''bravais'', &
       ''buildsurf'', ''newclusurf'', ''newcluimp'''
      error stop
    end if
  end subroutine check_pre_processing

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Check availability for processing
  !
  !> @param[in] processing Type of processing. Allowed values: 'none'
  !---------------------------------------------------------------------------
  subroutine check_processing(processing)
    character(len=*),intent(in) :: processing
    if(processing /= 'none') then
      write(error_unit,*) '[calculation.check_processing]: &
       calculation%processing must be one of: ''none'''
      error stop
    end if
  end subroutine check_processing
end module calculation_mod
