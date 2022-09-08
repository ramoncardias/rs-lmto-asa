!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Lattice
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
!> Module to handle structural properties
!------------------------------------------------------------------------------


module lattice_mod

  use control_mod
  use string_mod, only : clean_str
  use math_mod, only: cross_product
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use precision_mod, only: rp
  implicit none

  private

  !> Module's main structure
  type, public :: lattice
    ! General variables

    !> Lattice parameter (in \f$Å\f$)
    !>
    !> Lattice parameter (in \f$Å\f$)
    real(rp) :: alat
    
    !> Sphere radius to cut cluster (in \f$Å\f$)
    !>
    !> Use \ref r2 = \ref ct\f$^2\f$. This radius (\ref ct) refers to the distance which includes all first neighbors,
    !> or all second nearest neighbors, etc. (use \ref ct and \ref r2 including 5 th neighs. to run newclu.x)
    !> Example: Pd fcc: let \f$R_1\f$ the distance to first neighbors, and \f$R_2\f$ the distance to second neighbors,
    !> 
    !> \f$ R_1 = a \frac{\sqrt{2}}{2} = 2.7506 \, Å \f$,
    !>
    !> then,
    !>
    !> \f$ R_1^2 = 7.566 \, Å^2 \f$;
    !>
    !> and
    !> 
    !> \f$R_2 = a = 3.89 \, Å\f$
    !> 
    !> then,
    !>
    !> \f$R_2^2 = 15.13 \, Å^2\f$.
    !>
    !> Since \ref r2 \f$= 13\f$ is between \f$R_1^2\f$ and \f$R_2^2\f$ it will include all first neighs., but not second
    !> neighs.
    real(rp) :: r2

    !> TODO
    real(rp) :: celldm
    
    !> Wigner Seitz Radius (in \f$Å\f$)
    !>
    !> Wigner Seitz Radius (in \f$Å\f$)
    real(rp) :: wav

    !> Cell volume
    real(rp) :: vol

    !> Near neighbors cut radius (in \f$Å\f$)
    !>
    !> Near neighbors cut radius (in \f$Å\f$). See \ref r2 description for more details.
    real(rp), dimension(:), allocatable :: ct                  

    !> Auxiliar variables to set the correct ntype, nbulk and ntot
    integer :: nbulk_bulk

    !> Number of atoms of type bulk.
    !>
    !> Number of atoms of type bulk.
    !>
    !> Values:
    !> 
    !> - \ref nbulk \f$ = 0\f$: for bulk
    !> - \ref nbulk \f$ = 1\f$: for impurities in bulk
    !> - \ref nbulk \f$ = 1\f$: for surface withou defects
    !> - \ref nbulk \f$ = 6\f$: for for a system with defects on surface, where this surface has been calculated with 5 layers plus bulk.
    integer :: nbulk

    !> \ref ntot \f$= 1\f$ for fcc and bcc, without relaxation; \ref ntot \f$= 2\f$ for hcp.
    !>
    !> - \ref ntot \f$= 1\f$ for fcc and bcc, without relaxation;
    !> - \ref ntot \f$= 2\f$ for hcp.
    !>
    !> TODO: poor description
    integer :: ntot

    !> Number of inequivalent atoms.
    !>
    !> Number of inequivalent atoms.
    !>
    !> Always equal to \ref nbulk + \ref nrec
    !>
    !> Usage:
    !>
    !> - \ref ntype \f$= 1\f$: bulk material
    !> - \ref ntype \f$= 2\f$: impurity embedded in bulk (single site calculation)
    !> - \ref ntype \f$= 3\f$: impurity (at1) ambedded in bulk (plus nearest neighs. at2)
    !> - \ref ntype \f$= 6\f$: free surface system with 5 layers (at1, at2, at3, at4, at5) (since \ref nbulk \f$= 1\f$)
    !> - \ref ntype \f$= 7\f$: adatom on surface (single site calculation), where the free surface were converged with 5 layers + bulk (since \ref nbulk \f$= 6\f$)
    integer :: ntype
    
    !> Number of atoms to be considered in the recursion.
    !>
    !> Number of atoms to be considered in the recursion (for bulk \ref nrec\f$= 1\f$, for surface with 5 layers being calculated self-con. \ref nrec\f$= 5\f$, etc.)
    integer :: nrec

    !> TODO
    !> Clust size
    !>
    !> Clust size
    integer :: kk

    !> TODO
    !> Clust coordinates
    !>
    !> Clust coordinates
    real(rp), dimension(:,:), allocatable :: cr
    
    !> TODO
    !> Clust coordinates
    !>
    !> Clust coordinates
    real(rp), dimension(:,:), allocatable :: crd
    !> Clust atom number
    integer, dimension(:), allocatable :: ham_i
    !> TODO
    !> Atoms coordinates of the primitive cell
    !>
    !> Atoms coordinates of the primitive cell
    real(rp), dimension(:,:), allocatable :: primcell

    !> Number of shells/atoms to distribute charge
    !>
    !> Number of shells/atoms to distribute charge
    integer :: nbas

    !> Reduced \ref nbas
    !>
    !> Reduced \ref nbas
    integer :: reduced_nbas

    !>TODO
    !> Nmax
    integer :: nmax

    !> TODO
    !> Atom number for calculation.
    !>
    !> Vector of size \ref ntot.
    !>
    !> - \f$iu = 1\f$: for bulk material
    !> - Surface: this number IU is chosen in the clust file, looking for a site which can characterize the bulk layers, i.e. far from surface sites. The output on screen from buildsurf.x program gives this number.
    !> - Defects on surface: number given in newclu.x output (on screen).
    integer, dimension(:), allocatable :: iu

    !> TODO
    !> Atom number for calculation.
    !>
    !> Vector of size \ref ntot.
    !>
    !> - \f$ib = 1\f$: for bulk material
    !> - Surface: this number IB is chosen in the clust file, looking for a site which can characterize the bulk layers, i.e. far from surface sites. The output on screen from buildsurf.x program gives this number.
    !> - Defects on surface: number given in newclu.x output (on screen).
    integer, dimension(:), allocatable :: ib

    !> Refers to the atoms that will be calculated self-consistently, or to represent all equivalent atoms in the same neighboring shell.
    !>
    !> Refers to the atoms (in the clust file) that will be calculated self-consistently, or to represent all equivalent atoms in the same neighboring shell.
    !>
    !> - Those numbers are given by newclu.x output (on screen), for defects on surface and by buildsurf.x output (atomch.d file) for surface systems.
    integer, dimension(:), allocatable :: irec

    !> List containing the number of each atom inside the clust file
    integer, dimension(:), allocatable :: atlist

    !> Neighbouring map for each atom type
    integer, dimension(:,:), allocatable :: nn

    !> Structure constant
    complex(rp), dimension(:,:,:,:), allocatable :: sbar 
    !> Vectors in the structure constant
    real(rp), dimension(:,:), allocatable :: sbarvec
    ! Variables to build clust for bulk calculation

    !> TODO
    !> Cut radius
    !>
    !> Cut radius
    real(rp) :: rc

    !> Primitive vectors in units of lattice parameter \ref alat
    !>
    !> Primitive vectors in units of lattice parameter \ref alat
    real(rp), dimension(3,3) :: a
    
    !> TODO
    !> TBD
    integer, dimension(:), allocatable :: izp
    
    !> TODO
    integer, dimension(:), allocatable :: iz

    !> TODO
    integer, dimension(:), allocatable :: num
    
    !> TODO
    integer, dimension(:), allocatable :: no 
    
    !> TODO
    !> Paramemter that determines the clust size before cut. Should be large
    !enough, rarely changed. Default (ndim = 9900000, npe = 49)
    integer :: ndim, npe

    !> TODO
    !> Crystal symmetry. Options are 'bcc', 'fcc', 'hcp' and 'nsy'
    !>
    !> Crystal symmetry. Options are 'bcc', 'fcc', 'hcp' and 'nsy'
    character(len=4) :: crystal_sym

    ! Variables to build clust for surface calculation
    !> TODO
    !> Surface parameters
    real(rp) :: zmin, zmax, zstep              
    
    !> TODO
    !> Layer coordinate
    real(rp), dimension(:), allocatable :: z

    !> TODO
    !> Surface symmetry. Options are '111', '110' and '001'.
    !>
    !> Surface symmetry. Options are '111', '110' and '001'.
    character(len=10) :: surftype      

    !> TODO
    !> Number of layers
    !>
    !> Number of layers
    integer :: nlay

    !> TODO
    !> Surface indexes 
    integer :: dx, dy, dz, dw

    !> TODO
    !> TDB
    integer, dimension(:), allocatable :: izpsurf, izsurf, nosurf
    
    !> TODO
    !> Clust variables from surface   
    real(rp), dimension(:,:), allocatable :: crsurf     

    ! Variables to build clust for impurity calculation
    !> New coordinates
    real(rp), dimension(:,:), allocatable :: inclu
    
    !> TDB
    integer, dimension(:), allocatable :: izpo                  
    
    !> Impurity number
    integer :: nclu
    
    !> Clust variables from impurity
    real(rp), dimension(:,:), allocatable :: acr

    !> TODO
    integer, dimension(:), allocatable :: reduced_acr

    type(control), pointer :: control
  contains
    procedure :: build_from_file
    procedure :: restore_to_default
    procedure :: bravais
    procedure :: build_data
    procedure :: build_clusup
    procedure :: build_surf
    procedure :: newclu
    procedure :: structb
    procedure :: atomlist
    procedure, private :: dbar1
    procedure :: clusba
    procedure :: calculate_nbas
    procedure :: print_state
    procedure :: print_state_full
    procedure, private :: check_all
    final :: destructor
  end type lattice 

  interface lattice 
    procedure :: constructor
  end interface lattice

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] fname Namelist file
  !> @return type(lattice)
  !---------------------------------------------------------------------------
  function constructor(fname,control_obj) result(obj)
    type(lattice) :: obj
    type(control), target, intent(in) :: control_obj
    character(len=*), intent(in) :: fname

    obj%control => control_obj
    
    call obj%restore_to_default()
    call obj%build_from_file(fname)
  end function constructor
  
  
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine destructor(this)
    type(lattice) :: this
    if(allocated(this%ct)) deallocate(this%ct)
    if(allocated(this%cr)) deallocate(this%cr)
    if(allocated(this%crd)) deallocate(this%crd)
    if(allocated(this%iu)) deallocate(this%iu)
    if(allocated(this%ib)) deallocate(this%ib)
    if(allocated(this%irec)) deallocate(this%irec)
    if(allocated(this%izp)) deallocate(this%izp)
    if(allocated(this%iz)) deallocate(this%iz)
    if(allocated(this%num)) deallocate(this%num)
    if(allocated(this%no)) deallocate(this%no)
    if(allocated(this%z)) deallocate(this%z)
    if(allocated(this%izpsurf)) deallocate(this%izpsurf)
    if(allocated(this%izsurf)) deallocate(this%izsurf)
    if(allocated(this%nosurf)) deallocate(this%nosurf)
    if(allocated(this%crsurf)) deallocate(this%crsurf)
    if(allocated(this%inclu)) deallocate(this%inclu)
    if(allocated(this%izpo)) deallocate(this%izpo)
    if(allocated(this%acr)) deallocate(this%acr)
    if(allocated(this%atlist)) deallocate(this%atlist)
    if(allocated(this%nn)) deallocate(this%nn)
    if(allocated(this%sbar)) deallocate(this%sbar)
    if(allocated(this%sbarvec)) deallocate(this%sbarvec)
  end subroutine destructor


  ! Member functions
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Read parameters from input file
  !
  !> @param[in] fname Namelist file
  !---------------------------------------------------------------------------
  subroutine build_from_file(this,fname) 
    class(lattice),intent(inout) :: this
    character(len=*),intent(in) :: fname
    ! Namelist for general parameters
    integer, dimension(:), allocatable :: iu, ib, irec
    real(rp), dimension(:), allocatable :: ct                 
    ! Namelist variables for bulk
    real(rp) :: rc, r2, alat, celldm, wav
    real(rp), dimension(3,3) :: a
    integer, dimension(:), allocatable :: izp, no
    real(rp), dimension(:,:), allocatable :: crd
    integer :: ndim, npe
    character(len=4) :: crystal_sym
    ! Namelist variable for surface
    character(len=10) :: surftype
    integer :: nlay
    ! Namelist variables for impurity
    real(rp), dimension(:,:), allocatable :: inclu
    integer :: nclu
    ! variables associated with the reading processes
    integer :: iostatus, funit

    ! Namelist
    namelist /lattice/ crystal_sym, rc, a, izp, no, crd, ndim, npe,& ! Bulk parameters
                       surftype, nlay,&
                       inclu, nclu, & ! Impurity varables 
                       alat, celldm, wav        ! General parameters 

    ! Save previous values
    ! Bulk initialization
    ndim = this%ndim
    npe = this%npe
    rc = this%rc
    crystal_sym = this%crystal_sym
    a = this%a
    wav = this%wav
    celldm = this%celldm

    if(size(this%izp) .ne. this%ndim) then
      deallocate(this%izp)
      allocate(this%izp(this%ndim))
    end if
    if(size(this%no) .ne. this%ndim) then
      deallocate(this%no)
      allocate(this%no(this%ndim))
    end if
    if(size(this%crd) .ne. 3*this%ndim) then
      deallocate(this%crd)
      allocate(this%crd(3,this%ndim))
    end if

    call move_alloc(this%izp,izp)
    call move_alloc(this%no,no)
    call move_alloc(this%crd,crd)

    ! Impurity initialization
    nclu = this%nclu
    if(size(this%inclu) .ne. 3*this%nclu) then
      deallocate(this%inclu)
      allocate(this%inclu(this%nclu,3))
    end if
    call move_alloc(this%inclu,inclu)    

    ! Surface initialization
    surftype = this%surftype
    nlay = this%nlay

    open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
    if(iostatus /= 0) then
      write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
      error stop
    endif
                   
    read(funit,nml=lattice,iostat=iostatus) 

    if(size(izp) .ne. ndim) then
      deallocate(izp)
      allocate(izp(ndim))
    end if
    if(size(no) .ne. ndim) then
      deallocate(no)
      allocate(no(ndim))
    end if
    if(size(crd) .ne. 3*ndim) then
      deallocate(crd)
      allocate(crd(3,ndim))
    end if
    if(size(inclu) .ne. 3*nclu) then
      deallocate(inclu)
      allocate(inclu(nclu,3))
    end if
    
    rewind(funit)
    read(funit,nml=lattice,iostat=iostatus)
    if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
      write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
      write(error_unit,'("iostatus = ",I0)') iostatus
    endif
    close(funit)

    ! General intialization
    r2 = ceiling(alat)**2

    this%r2 = r2
    this%alat = alat
    this%celldm = celldm
  
    ! Reads Wigner-Seitz radius if available
    this%wav = wav
    ! Bulk initialization
    this%ndim = ndim
    this%npe = npe
    this%rc = rc
    this%crystal_sym = crystal_sym
    this%a = a
    call move_alloc(izp,this%izp)
    call move_alloc(no,this%no)
    call move_alloc(crd,this%crd)

    ! Impurity initialization
    this%nclu = nclu
    call move_alloc(inclu,this%inclu)    

    ! Surface initialization
    this%surftype = surftype
    this%nlay = nlay

    ! checks
    call this%check_all()

  end subroutine build_from_file     

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Build structure according to crystal symetry provided
  !---------------------------------------------------------------------------
  subroutine build_data(this)
    class(lattice), intent(inout) :: this
    !> Local variables
    real(rp), dimension(3,3) :: a
    integer :: i, j

    select case(this%crystal_sym)
    case('bcc')
      a(:,1) = [-0.50000000,  0.50000000,  0.50000000]
      a(:,2) = [ 0.50000000, -0.50000000,  0.50000000]
      a(:,3) = [ 0.50000000,  0.50000000, -0.50000000]
      this%crd(:,1) = [0.0, 0.0, 0.0]
      this%izp(1) = 1
      this%no(1) = 1
      this%nbulk_bulk = 1
      this%ntot = 1
      if(this%control%calctype=='B')then
        this%nrec = this%nbulk_bulk
        this%nbulk = 0
        this%ntype = 1
        this%nbas = this%ntot
        allocate(this%ib(this%ntot),this%iu(this%ntot),this%irec(this%nrec),this%ct(this%ntype))
        this%iu(1) = 1
        this%ib(1) = 1
        this%irec(1) = 1
        this%ct(:) = this%alat+0.1d0
      end if
      this%a = a
    case('fcc')
      a(:,1) =  [0.00000000,  0.50000000,  0.50000000]
      a(:,2) =  [0.50000000,  0.00000000,  0.50000000]
      a(:,3) =  [0.50000000,  0.50000000,  0.00000000]
      this%crd(:,1) = [0.0, 0.0, 0.0]
      this%izp(1) = 1
      this%no(1) = 1
      this%nbulk_bulk = 1
      this%ntot = 1
      if(this%control%calctype=='B')then
        this%nrec = this%nbulk_bulk
        this%nbulk = 0
        this%ntype = 1
        this%nbas = this%ntot
        allocate(this%ib(this%ntot),this%iu(this%ntot),this%irec(this%nrec),this%ct(this%ntype))
        this%iu(1) = 1
        this%ib(1) = 1
        this%irec(1) = 1
        this%ct(:) = this%alat+0.1d0
      end if
      this%a = a
    case('fcc2')
      a(:,1) =  [0.00000000,  0.50000000,  0.50000000]
      a(:,2) =  [0.50000000,  0.00000000,  0.50000000]
      a(:,3) =  [0.50000000,  0.50000000,  0.00000000]
      this%crd(:,1) = [0.00,0.00,0.00]
      this%crd(:,2) = [0.50,0.50,0.50]
      this%izp(:) = [1,2]
      this%no(:) = [1,2]
      this%nbulk_bulk = 2
      this%ntot = 2
      if(this%control%calctype=='B')then
        this%nrec = this%nbulk_bulk
        this%nbulk = 0
        this%ntype = 2
        this%nbas = this%ntot
        allocate(this%ib(this%ntot),this%iu(this%ntot),this%irec(this%nrec),this%ct(this%ntype))
        this%iu(:) = [1,2]
        this%ib(:) = [1,2]
        this%irec(:) = [1,2]
        this%ct(:) = this%alat-1.5d0
        this%r2 = this%ct(1)**2
      end if
      this%a = a
    case('fcc3')
      a(:,1) =  [0.00000000,  0.50000000,  0.50000000]
      a(:,2) =  [0.50000000,  0.00000000,  0.50000000]
      a(:,3) =  [0.50000000,  0.50000000,  0.00000000]
      this%crd(:,1) = [0.00,0.00,0.00]
      this%crd(:,2) = [0.50,0.50,0.50]
      this%crd(:,3) = [0.25,0.25,0.25]
      this%crd(:,4) = [-.25,-.25,-.25]
      this%izp(:) = [1,2,3,4]
      this%no(:) = [1,2,3,4]
      this%nbulk_bulk = 4
      this%ntot = 4
      if(this%control%calctype=='B')then
        this%nrec = this%nbulk_bulk
        this%nbulk = 0
        this%ntype = 4
        this%nbas = this%ntot
        allocate(this%ib(this%ntot),this%iu(this%ntot),this%irec(this%nrec),this%ct(this%ntype))
        this%iu(:) = [1,2,3,4]
        this%ib(:) = [1,2,3,4]
        this%irec(:) = [1,2,3,4]
        this%ct(:) = this%alat-1.5d0
        this%r2 = this%ct(1)**2
      end if
      this%a = a
    case('hcp')
      if (this%celldm==0.0d0) then
         print *, 'WARNING: hcp structure using as default c/a = 1.633'
         this%celldm = 1.633d0
      end if
      a(:,1) =  [1.00000000,  0.00000000,  0.00000000]
      a(:,2) =  [-.50000000,  0.86603000,  0.00000000]
      a(:,3) =  [0.00000000,  0.00000000,  0.00000000]
      a(3,3) = this%celldm
      this%crd(:,1) = [0.0, 0.0, 0.0]
      this%crd(:,2) = [-.5, -.28868000, 0.00000000]
      this%crd(3,2) = (0.5d0)*this%celldm
      this%nbulk_bulk = 2
      this%ntot = 2
      this%izp(:) = [1,2]
      this%no(:) = [1,2]
      if(this%control%calctype=='B')then
        this%nrec = this%nbulk_bulk
        this%nbulk = 0
        this%ntype = 2
        this%nbas = this%ntot
        allocate(this%ib(this%ntot),this%iu(this%ntot),this%irec(this%nrec),this%ct(this%ntype))
        this%iu(:) = [1,2]
        this%ib(:) = [1,2]
        this%irec(:) = [1,2]
        this%ct(:) = this%alat+0.1d0
      end if
      this%a = a
    case('nsy')
      ! reads from input
    end select
    ! Volume in cubic Angstroms
    this%vol = abs(dot_product(a(:,3),cross_product(a(:,1),a(:,2))))*(this%alat**3)
    if(this%wav.eq.0) then
       this%wav = (this%vol/((16.0d0/3.0d0)*atan(1.0d0)*this%ntot))**(1.0d0/3.0d0)
    end if
    if(this%control%calctype=='B'.or.this%control%calctype=='S') this%nmax = 0
  end subroutine build_data

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Reset all members to default
  !---------------------------------------------------------------------------
  subroutine restore_to_default(this)
    class(lattice) :: this

    this%ndim = 9900000
    this%npe = 49
    this%nclu = 0
    this%surftype = 'none'
    this%nlay = 0
    this%wav = 0
    this%celldm = 0.0d0

    allocate(this%izp(this%ndim),this%no(this%ndim))
    allocate(this%crd(3,this%ndim))
    allocate(this%inclu(this%nclu,3))

  end subroutine restore_to_default

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !---------------------------------------------------------------------------
  subroutine bravais(this)
    class(lattice),intent(inout) :: this
    ! Local variables
    real(rp) :: rc, rs, lc
    integer, dimension(:), allocatable :: iz, num
    real(rp), dimension(:,:), allocatable :: cr, crbravais 
    integer :: npe, ndim, nx, ny, nz, npr, l, n, i, nl, k, kk
    logical :: isopen
    integer :: iostatus

    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'lattice%bravais, file clust: Unit 10 is already open'
      error stop
    else
      open(unit=10,file='clust')
    endif

    n = this%ntot
    ! Defining a size-like parameter for the clust before the cut
    npe = this%npe
    ndim = this%ndim
    rc = this%rc 
    allocate(iz(ndim),num(ndim))
    allocate(cr(3,ndim),crbravais(3,ndim))
   ! Defining the cut radius
    rs = (0.8d0*int((npe/(1.0d0))/2.0d0))**2
    rs = dmin1(rs,rc)        
    
    npr = int((ndim/(n*1.0d0))**(1.0d0/3.0d0))
  
    crbravais(:,:) = this%crd(:,:)
    lc = (npr+1)/2
    l=n
    do i=1,l                    
      do nx=1,npr
        do ny=1,npr
          do nz=1,npr
            if(nx.eq.lc.and.ny.eq.lc.and.nz.eq.lc)go to 13
            n=n+1
!  ...........Verifies dimension NDIM.........................
            if(n.gt.ndim)then
              write(6,600) n
              stop
            end if
! ...............................................................
     crbravais(1,n)=this%crd(1,i)+(nx-lc)*this%a(1,1)+(ny-lc)*this%a(1,2)+(nz-lc)*this%a(1,3)
     crbravais(2,n)=this%crd(2,i)+(nx-lc)*this%a(2,1)+(ny-lc)*this%a(2,2)+(nz-lc)*this%a(2,3)
     crbravais(3,n)=this%crd(3,i)+(nx-lc)*this%a(3,1)+(ny-lc)*this%a(3,2)+(nz-lc)*this%a(3,3)
            this%no(n)=this%no(i)
            this%izp(n)=this%izp(i)
 13         continue
          end do
        end do
      end do
    end do
    nl = l*(npr**3)

    kk=0
    if(rc==0.0d0)rs=npr**3
    
    do i=1,nl
      call cut(i,l,ndim,crbravais,cr,this%izp,iz,num,this%no,rs,kk)
    end do

    if(int(kk/2)/=kk/2.d0) kk = kk-1

    write(10,300)kk
    do k=1,kk,2
      write(10,200)(cr(1,k+i-1),cr(2,k+i-1),cr(3,k+i-1),iz(k+i-1),num(k+i-1),i=1,2)
    end do
!    write(900,*) kk
!    do k=1,kk
!      write(900,*) 'Fe', cr(:,k)
!    end do

200 format(3F14.8,2I4,3F14.8,2I4)
300 format(3X,"II =",I7)
600 format(1X,'K =',I7,5X,'redefine the dimension ndim')

    this%kk = kk
    call move_alloc(cr,this%cr)
    call move_alloc(iz,this%iz)                    
    call move_alloc(num,this%num)                 
    close(10)
  end subroutine bravais

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !---------------------------------------------------------------------------
  subroutine build_clusup(this)
    class(lattice),intent(inout) :: this
    character(len=10) :: surftype
    ! Local variables
    real(rp), dimension(:), allocatable :: z
    integer, dimension(:), allocatable :: izpsurf            
    real(rp) :: ds, ds2, new, one
    integer :: i, j, n
    
    select case(this%crystal_sym)
      case('hcp')
        read(this%surftype, *) this%dx, this%dy, this%dz, this%dw
        if ((this%dz).ne.(-1*(this%dx+this%dy))) then
          this%dz = -1*(this%dx + this%dy)
          write(error_unit,*) 'Hexagonal Miller indices not right. Does it should be [', &
          this%dx, this%dy, this%dz, this%dw,']?'
          error stop
        end if
        this%dx = (2*this%dx+this%dy)
        this%dy = (this%dx+2*this%dy)
        this%dz = this%dw
      case default
        read(this%surftype, *) this%dx, this%dy, this%dz
    end select

!  ............Find zstep, zmin, zmax....................
    ds=1000.0d0 ; ds2=1000.0d0
    select case(this%crystal_sym)
      case('hcp')
        do i=1,this%kk
           do j=1,this%kk
              new = this%dx*this%cr(1,i)+this%dy*this%cr(2,i)+this%dz*this%cr(3,i)
              one = this%dx*this%cr(2,j)+this%dy*this%cr(2,j)+this%dz*this%cr(3,j)
              if ((abs(new-one).gt.1.d-6).and.(abs(new-one).le.ds)) then
                 ds = abs(new-one)
              end if
              if ((abs(new).le.ds2)) then
                 ds2 = abs(new)
              end if
           end do
        end do
      case default
        do i=1,this%kk
           new = this%dx*this%cr(1,i)+this%dy*this%cr(2,i)+this%dz*this%cr(3,i)
           if ((abs(new).gt.1.d-6).and.(abs(new).le.ds)) then
              ds = abs(new)
           end if
           if ((abs(new).le.ds2)) then
              ds2 = abs(new)
           end if
        end do
      end select
    
    this%zstep = ds
    this%zmin = ds2 - this%zstep
    this%zmax = ds2 + 15*this%zstep
    
!  ......................................................
    
    n = int((this%zmax-this%zmin)/this%zstep)+1
    
    allocate(z(n),izpsurf(n))
    do i=1,n
       z(i)=this%zmin+((i-1)*this%zstep)
    end do
    call move_alloc(z,this%z)
  
    do i=1,n
       izpsurf(i) = mod(i,this%ntot)+1
    end do
  
    j = this%ntot
    do i=1,this%nlay
       j = j+1
       izpsurf(i) = j
    end do
    call move_alloc(izpsurf,this%izpsurf)
    
    if(this%control%calctype=='S')then 
      this%ntype = this%nbulk_bulk+this%nlay
      this%nbulk = this%nbulk_bulk
      this%nrec = this%ntype-this%nbulk
      this%nbas = 17        
      allocate(this%ib(this%nbulk),this%irec(this%nrec),this%iu(this%ntot),this%ct(this%ntype))
      this%ct(:) = this%alat+0.1d0            
    end if
    
    this%surftype = clean_str(this%surftype)
  end subroutine build_clusup

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !---------------------------------------------------------------------------
  subroutine build_surf(this)
    class(lattice),intent(inout) :: this
    ! Local variables
    integer :: i, j, n, k, kk, ichoice, icont, ichoicetype          
    real(rp), dimension(this%kk) :: crh, crhd
    real(rp) :: disf, disi_min, disi
    real(rp), dimension(:,:), allocatable :: crsurf
    integer, dimension(:), allocatable :: izp, no, ichoicen, ichoicetypen
    logical :: isopen

    n = int((this%zmax-this%zmin)/this%zstep)+1
    kk = this%kk
    allocate(crsurf(3,kk),izp(kk),no(kk),ichoicen(kk),ichoicetypen(kk))
    icont = 0
    crh = 0.0d0
    crhd = 0.0d0
    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'lattice%build_surf, file clust: Unit 10 is already open'
      error stop
    else
      open(unit=10,file='clust')
    endif

    do i=1,n
      disf = 1.0d0
      disi_min = sqrt(this%z(i)**2)+0.5d0
      do k=1,kk
        crh(k) = 0.0d0
        crh(k) = this%dx*this%cr(1,k)+this%dy*this%cr(2,k)+this%dz*this%cr(3,k)
        if(abs(crh(k)-this%z(i))<1.0d-6)then
          icont = icont+1
          do j=1,3
            crsurf(j,icont) = this%cr(j,k)
          end do
          crhd(icont) = this%dx*crsurf(1,icont)+this%dy*crsurf(2,icont)+this%dz*crsurf(3,icont) 
          if(abs(crhd(icont)-this%z(i))<1.0d-6)then
            izp(icont) = this%izpsurf(i) 
            disi = norm2(crsurf(:,icont))
            if(disi<disi_min)then
              disi_min = disi
              ichoice = icont
              ichoicetype = izp(icont)   
            end if
            no(icont) = this%num(k) 
          else
          end if
        end if
      end do
      ichoicetypen(i) = ichoicetype
      ichoicen(i) = ichoice
    end do


    if(this%control%calctype=='S')then
      do i=1,this%nrec
        this%irec(i) = ichoicen(i)
      end do
      do i=1,this%ntot
        this%iu(ichoicetypen(i+this%nrec)) = ichoicen(i+this%nrec)
      end do
      do i=1,this%nbulk
        this%ib(ichoicetypen(i+this%nrec)) = ichoicen(i+this%nrec)
      end do
    end if

    if(int(icont/2)/=icont/2.d0) icont = icont-1

    write (10,10004) icont
    do k= 1,icont-1,2
      write (10,10002) &
          crsurf(1,k), crsurf(2,k), crsurf(3,k), izp(k), no(k), crsurf(1,k+1), crsurf(2,k+1), &
          crsurf(3,k+1), izp(k+1), no(k+1)
    end do

    deallocate(this%cr,this%iz,this%num)

    call move_alloc(crsurf,this%cr)
    call move_alloc(izp,this%iz)
    call move_alloc(no,this%num)

    this%kk = icont

    10002 format (3f14.6,2i4,3f14.6,2i4)
    10004 format (3x,"II =",i7)
    close(10)
  end subroutine build_surf

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Adds the cluster/impurity atoms in 'inclu' to 'clu0' and outputs to 'clust'
  !---------------------------------------------------------------------------
  subroutine newclu(this)
    class(lattice),intent(inout) :: this
    ! Local variables 
    integer :: ndi, nnmx, ncnt, nmax                    
    logical :: isopen
    integer :: i,j, kk, k, ntypecount, ireccount, inclucheck
    integer, dimension(:), allocatable ::  ibulk                        
    integer, dimension(:), allocatable :: izpo, izp, no, nnmax, izimp, noimp           
    integer, dimension(:,:), allocatable :: nn, nn2            
    real(rp) :: nnscale
    real(rp), dimension(:,:), allocatable :: acr, crd, crimp             
    real(rp), dimension(:), allocatable :: ctnew   


    ! Set clust variables
    this%nbulk = this%nbulk_bulk+this%nlay
    this%ntype = this%nbulk_bulk+this%nlay+this%nclu
    this%nrec = this%ntype-this%nbulk
    ! Set the clust dimension
    kk = this%kk  
    ndi = 150000
    nnmx = 200
    ! Allocating clust dimension
    allocate(this%ib(this%nbulk),this%iu(this%ntot),this%irec(this%nrec),this%ct(this%ntype))
    allocate(ctnew(this%ntype),ibulk(this%nbulk))
    allocate(nn(ndi,nnmx),nn2(ndi,nnmx))
    allocate(izpo(kk),izp(kk),no(kk),nnmax(kk),izimp(kk),noimp(kk))
    allocate(acr(kk,7),crd(3,kk),crimp(3,kk))
    ! Setting ct values for impurity
    this%ct(:) = this%alat+0.1d0
    ! Identify impurity atoms from 'inclu'
    do i=1,kk
      izpo(i)=this%iz(i)
    end do
    ntypecount = this%nbulk
    inclucheck = 0
    do j=1,this%nclu
      do i=1,kk
        if(abs(this%cr(1,i)-this%inclu(j,1))<1.0d-6.and.&
           abs(this%cr(2,i)-this%inclu(j,2))<1.0d-6.and.&
           abs(this%cr(3,i)-this%inclu(j,3))<1.0d-6)then
          inclucheck = inclucheck + 1
          ntypecount = ntypecount + 1
          this%iz(i) = ntypecount 
        end if
      end do
    end do

    if(inclucheck/=this%nclu)then
      write(error_unit,'(a)') 'Atoms chosen for impurity were not found inside the clust. Please, check the inclu input'
      error stop
    end if

    acr=0.0d0
    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'lattice%newclu, file clust: Unit 10 is already open'
      error stop
    else
      open(unit=10,file='clust')
    endif

    open(unit=11,file='outnewclu')

    do k=1,kk
      do i=1,3
        acr(k,i) = this%cr(i,k)
        acr(k,6) = acr(k,6)+(this%cr(i,k)-this%inclu(1,i))**2
      end do
      acr(k,4) = (this%iz(k))
      acr(k,5) = (this%num(k))
      acr(k,7) = (izpo(k)) 
    end do
    call bubble(this%nclu,this%nclu,acr(1:this%nclu,1:7),7,4)
    call bubble(kk-this%nclu,kk-this%nclu,acr(this%nclu+1:kk,1:7),7,6)
    
     write (10,10004) kk
    do k = 1,kk,2
      write (10,10002) &
        ( &
        acr(k+i-1,1),acr(k+i-1,2),acr(k+i-1,3),int(acr(k+i-1,4)) &
        ,int(acr(k+i-1,5)), i = 1,2)
    end do

    rewind(10)

    call leia(this%alat,kk,crd,izp,no,10)

    close(10)

    nnscale = 0.95d0
 
    ctnew(:) = this%ct(:)*nnscale   

    call nncal(ctnew,crd,3,kk,izp,nn,ndi,nnmx,mapa,this%ntot) 
    nnmx = 200
    call nncal(this%ct,crd,3,kk,izp,nn2,ndi,nnmx,mapa,this%ntot)
    nnmx = 200

    do i=1,kk
      if(acr(i,4)>this%nbulk.and.acr(i,4)<=this%ntype) then
        do j=2,nn2(i,1)
          if(acr(nn2(i,j),4)<=this%nbulk)  acr(nn2(i,j),4)=2000+acr(i,7)
        end do
      end if
    end do
    do i=1,kk
      if(acr(i,4)>this%nbulk.and.acr(I,4)<=this%ntype) then
        do j=2,nn(i,1)
          if(acr(nn(i,j),4)<=this%nbulk.or.acr(nn(i,j),4)>2000) acr(nn(i,j),4)=1000+acr(i,7)
        end do
      end if
    end do
    do i=1,kk
      if(acr(i,4)==1) acr(i,4)=4000+acr(I,7)
    end do
    do i=1,kk
      if(acr(i,4)<=this%nbulk) acr(i,4)=3000+acr(I,7)
    end do
    call bubble(kk,kk,acr(1:kk,1:7),7,4)
    ncnt=0
    do I=1,kk
      if(acr(i,4)<2000) ncnt=ncnt+1
    end do
    do i=1,kk
      if(acr(i,4)>this%ntype) acr(i,4)=acr(i,7)
    end do
    call bubble(kk-ncnt,kk-ncnt,acr(ncnt+1:kk,1:7),7,6)
    open(unit=10,file='clust')
    write (10,10004) kk
    do k = 1,kk,2
      write (10,10002) &
           ( &
           acr(k+i-1,1),acr(k+i-1,2),acr(k+i-1,3),int(acr(k+i-1,4)) &
           ,int(acr(k+i-1,5)), i = 1,2)
    end do
    rewind(10)
    call leia(this%alat,kk,crd,izp,no,10)
    call nncal(this%ct,crd,3,kk,izp,nn,ndi,nnmx,mapa,this%ntot)
    call outmap(11,izp,nn,no,ndi,nnmx,this%nrec)
    close(10)
    nmax=0
    do i=1,this%nrec
      do j=2,nn(i,1)
        if(nmax<nn(i,j)) nmax=nn(i,j)
      end do
    end do
    write(11,*) "--Info-for-bulcri------------------"
    do k=1,ncnt
      write(11,'(i5)') int(acr(k,7))   !,int(ACR(K,7))
    end do
    this%nbas = ncnt
    this%nmax = nmax
    write(11,*) "--Info-for-control-----------------"
    write(11,'(a6,i4,a6,i6,a6,i6,a6,i6)') 'NTYPE=',this%ntype, 'NMAX=',this%nmax,'NBAS=',this%nbas,'NREC=',this%nrec

    nnmax=0
    do i=nmax+1,kk
      if(nnmax(izp(i))<nn(i,1)) then
        nnmax(izp(i))=nn(i,1)
        ibulk(izp(i))=i
      end if
    end do
    write(11,'(50i7)') (nnmax(k),k=1,this%nbulk)
    write(11,'(50i7)') (ibulk(k),k=1,this%nbulk)
    do i=1,this%nbulk
      do j=2,nn(ibulk(i),1)
        if(nn(ibulk(i),j)==0) write(*,*) "Warning! Atom",ibulk(i)," is probably wrong."
      end do
    end do

    deallocate(this%cr,this%iz,this%num)

    do k=1,kk
     crimp(1:3,k) = acr(k,1:3)
    end do

    izimp(:) = int(acr(:,4))
    noimp(:) = int(acr(:,5))

    call move_alloc(izimp,this%iz)
    call move_alloc(noimp,this%num)
    call move_alloc(crimp,this%cr)
 
    do i=1,this%nbulk
      this%ib(i) = ibulk(i)
    end do

    do i=1,this%ntot
      this%iu(i) = ibulk(i)
    end do

    ireccount = 0
    do j=1,this%nclu
      do i=1,kk
        if(abs(this%cr(1,i)-this%inclu(j,1))<1.0d-6.and.&
           abs(this%cr(2,i)-this%inclu(j,2))<1.0d-6.and.&
           abs(this%cr(3,i)-this%inclu(j,3))<1.0d-6)then
           ireccount = ireccount + 1
           this%irec(ireccount) = i
        end if
      end do
    end do

    close(11)
    10002 format (3f14.6,2i4,3f14.6,2i4)
    10004 format (3x,"II =",i7)
  end subroutine newclu  

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Calculates the tight-binding structure constant matrix elements from
  !> information in 'control' and 'clust.
  !---------------------------------------------------------------------------
  subroutine structb(this)
    class(lattice),intent(inout) :: this
    ! Local variables
    integer :: i,ia,nr,ii,j,nm,np,nlim,nomx,ncut,kk,nnmx
    integer, dimension(:,:), allocatable :: nn
    integer, dimension(:), allocatable :: idnn
    real(rp), dimension(:,:,:), allocatable :: set
    real(rp), dimension(3) :: ret

    ! Open files
    open (12,file = 'map', form = 'unformatted')
    open (13,file = "sbar",FORM = "unformatted")
    open (16,file = "view.sbar")
    open (17,file = 'str.out')
    
    ! Clust parameters
    ncut = 9
    nnmx = 5250
    nomx = this%ntot
    kk = this%kk
    allocate(set(3,nomx,nnmx)); set = 0.0d0
    allocate(nn(kk,nnmx))
    allocate(idnn(nnmx))
    nm = nnmx
    write (17,*) 'irec',this%nrec, this%irec
    write (17,*) 'ndi=',kk
    write (17,10000) kk
    write (17,10001)
    write (17,10002) (i, (this%cr(j,i)*this%alat, j = 1,3), i = 1,max(this%nmax,this%ntype))
    call nncal(this%ct,this%cr*this%alat,3,kk,this%iz,nn,kk,nm,mapa,this%ntype)

    allocate(this%nn(this%kk,nm+1))

    do ii = 1,nm+1
      this%nn(:,ii) = nn(:,ii)     
    end do

    allocate(this%sbar(9,9,nm,this%ntot))

    write (17,*)'ndi=',kk
    write (17,*)'remd'
    call remd(this%cr*this%alat,this%num,this%iu,this%nn,kk,this%ntot,nomx,kk,nnmx,set,idnn,ret)
    write (17,*)'outmap',this%nmax,maxval(this%irec)
    call outmap(17,this%iz,this%nn,this%num,kk,nnmx,max(this%nmax,maxval(this%irec)))
    write (17,10003) kk, nm
    do ii=1,this%ntot
      ia = this%iu(ii)
      nr = this%nn(ia,1)
      write (17,'(1x,a,i5,a,i5)') 'Sbar atom no:',ii,' Ntot:',this%ntot
      call this%dbar1(ia,ncut*this%r2,this%wav,this%cr*this%alat,kk,kk,this%control%npold,nr,ii)
    end do

    10000 format (i5)
    10001 format (" LATTICE COORDINATES")
    10002 format (2(i5,3f8.4))
    10003 format (3i5)
    10004 format (3i5)
    10005 format (7x,i7)
  end subroutine structb


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Makes a list of atoms with their respective number inside the clust file 
  !---------------------------------------------------------------------------
  subroutine atomlist(this)
    class(lattice), intent(inout) :: this
    ! Local variables
    integer :: i, j

    allocate(this%atlist(this%ntype))
    allocate(this%ham_i(this%kk))

    do i=1,this%kk
      this%ham_i(i)=i
    end do

    j = 0
    do i = 1,this%nbulk
      j = j+1
      this%atlist(j) = this%ib(i)
    end do
    do i = 1,this%nrec
      j = j+1
      this%atlist(j) = this%irec(i)
    end do
  end subroutine atomlist

  ! Local library
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> TODO
  !---------------------------------------------------------------------------
  subroutine dbar1(this,ia,r2,wav,crd,nat,ndi,np,nr,ii)
    implicit none
    class(lattice), intent(inout) :: this
    ! Inputs
    integer, intent(in) :: ia,nat,ndi,np
    integer, intent(in) :: nr,ii
    real(rp), intent(in) :: r2,wav
    real(rp), dimension(3,ndi), intent(in) :: crd
    ! Local Scalars    
    integer :: i,j,k,m,na,nrl,nt
    real(rp), dimension(:),allocatable :: bet,wk
    real(rp), dimension(:),allocatable :: a
    real(rp), dimension(:,:),allocatable :: cr
    real(rp), dimension(:,:),allocatable :: s
    real(rp), dimension(:,:,:),allocatable :: sbar
    !
    ! External Calls    
    !external CLUSBA, MICHA

    nt=5250 !350
    allocate(cr(3,nt))
    allocate(sbar(np,np,nt))
    call this%clusba(r2,crd,ia,nat,ndi,nt)
    write (17,10000) nt
    write (17,10001) ((this%sbarvec(j,i), j = 1,3), i = 1,nt)
    nrl = np * nt
    na = (nrl*(nrl+1)) / 2
    allocate(a(na))
    allocate(bet(nrl))
    allocate(wk(nrl))
    allocate(s(nrl,nrl))
    call micha(wav,this%sbarvec,nt,np,nrl,na,sbar,a,wk,bet,s,ia,r2)

    ! Saving parameters to be used in the Hamiltonian build
    call this%clusba((r2/9.0d0),crd,ia,nat,ndi,nt)

    do m = 1,nt
      do i = 1,9
        do j = 1,9
          this%sbar(i,j,m,ii) = sbar(i,j,m)
        end do
      end do
    end do


    !do m=1,nt
    !  write(*,*) ia, m 
    !  do i=1,9
    !    write(*,'(9f10.4)')((real(this%sbar(i,j,m,ia))),j=1,9)
    !  end do
    !end do
    deallocate(a,bet,cr,wk,s,sbar)
    return
     
    10000 format (i5)
    10001 format (3f8.4)
  end subroutine dbar1

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> TODO
  !
  !> @param[in] r2
  !> @param[in] crd
  !> @param[in] ia
  !> @param[in] nat
  !> @param[in] ndi
  !> @param[inout] n
  !> @param[inout] cr
  !> @return type(calculation)
  !---------------------------------------------------------------------------
  subroutine clusba(this,r2,crd,ia,nat,ndi,n)
    implicit none
    class(lattice), intent(inout) :: this
    ! Inputs               
    integer, intent(in) :: ia,nat,ndi
    real(rp), intent(in) :: r2
    real(rp), dimension(3,ndi), intent(in) :: crd
    ! Output
    integer, intent(inout) :: n
!    real(rp), dimension(3,n), intent(inout) :: cr
    ! Local variables   
    integer :: i,ii,k,nn
    real(rp) :: s1
    real(rp), dimension(3) :: dum
 
    if(allocated(this%sbarvec)) deallocate(this%sbarvec)
    allocate(this%sbarvec(3,this%kk))
    this%sbarvec(:,:) = 0.0d0

    ii = 1
    do k = 1,3
      this%sbarvec(k,1) = 0.0d0
    end do
    do nn = 1,nat
      s1 = 0.0
      do i = 1,3
        dum(i) = (crd(i,nn)-crd(i,ia))**2
        s1 = s1 + dum(i)
      end do
      if (s1<r2 .and. s1>0.0001) then
        ii = ii + 1
        this%sbarvec(1,ii) = crd(1,nn) - crd(1,ia)
        this%sbarvec(2,ii) = crd(2,nn) - crd(2,ia)
        this%sbarvec(3,ii) = crd(3,nn) - crd(3,ia)
      end if
!!    if(ii>n) stop "Too large sbar cutoff, decrease NCUT in MAIN or increase NA
!!    in DBAR1."
    end do
    n = ii
  end subroutine clusba

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Program to make screened structure constants
  !> S(ME)=-0.5*S(OLE)   J(ME)=2*J(OLE)  QBAR(ME)=2*QBAR(OLE)
  !
  !> @param[in] rws
  !> @param[in] r
  !> @param[in] nr
  !> @param[in] nlm
  !> @param[in] nrl
  !> @param[in] na
  !> @param[inout] sbar
  !> @param[inout] a
  !> @param[inout] wk
  !> @param[inout] bet
  !> @param[inout] s
  !> @param[in] iclus
  !> @param[in] r2
  !---------------------------------------------------------------------------
  subroutine micha(rws,r,nr,nlm,nrl,na,sbar,a,wk,bet,s,iclus,r2)
    implicit none
    ! Inputs              .
    integer, intent(in) :: iclus,na,nlm,nr,nrl
    real(rp), intent(in) :: rws,r2
    real(rp), dimension(3,nr), intent(in) :: r
    ! Outputs
    real(rp), dimension(na), intent(inout) :: a
    real(rp), dimension(nrl), intent(inout) :: bet,wk
    real(rp), dimension(nrl,nrl), intent(inout) :: s
    real(rp), dimension(nlm,nlm,nr), intent(inout) :: sbar
    ! Local variables  
    real(rp) :: fak,pi
    real(rp), dimension(3) :: q
    ! External Calls   
    !external SHLDCH, STREZE
    ! Intrinsic Functions   
    intrinsic ATAN

    pi = 4.d0 * ATAN(1.d0)
    ! ----------------------------------
    call STREZE(rws,r,nr,s,nrl,nlm)
    ! -------------------------------------------
    fak = 2.d0
    !Original faktors
    q(1) = 0.3485d0 * fak
    q(2) = 0.05303d0 * fak
    q(3) = 0.010714d0 * fak
    ! Factors from LMTO47
    !q(1) = 0.33727d0 * fak
    !q(2) = 0.05115d0 * fak
    !q(3) = 0.01397d0 * fak
  !! ! New from LMTO47 (fcc Pt/Es)    data QM/0.3500000d0,0.0571667d0,0.0168070d0/
    !q(1) = 0.3500000d0* fak
    !q(2) = 0.0571667d0* fak
    !q(3) = 0.0168070d0* fak
    !     NA=(NRL*(NRL+1))/2
    call SHLDCH(r,nr,nlm,nrl,s,a,na,q,bet,wk,sbar,iclus,r2)
  end subroutine micha

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !>  Makes big matrix of all structure constants connecting sites R
  !
  !> @param[in] w
  !> @param[in] r
  !> @param[in] nr
  !> @param[inout] s
  !> @param[in] nrl
  !> @param[in] nlm
  !---------------------------------------------------------------------------
  subroutine STREZE(w,r,nr,s,nrl,nlm)
    implicit none
    ! Input                       
    integer, intent(in) :: nlm,nr,nrl
    real(rp), intent(in) :: w
    real(rp), dimension(3,nr), intent(in) :: r
    ! Output 
    real(rp), dimension(nrl,nrl), intent(inout) :: s
    ! Local variables    
    integer :: ilm,ir,irl0,jlm,jr,jrl0
    real(rp) :: rr,w1
    real(rp), dimension(3) :: dr
    real(rp), dimension(81) :: s0
    ! External calls      
    !external CANSO
    ! Intrinsic Functions   
    intrinsic SQRT

    if (nlm > 9) stop "**** CHANGE DIMS IN STRMAT"
    do ir = 1,nr
      irl0 = (ir-1) * nlm
      do jr = 1,nr
        jrl0 = (jr-1) * nlm
        dr(1) = (r(1,jr)-r(1,ir)) / w
        dr(2) = (r(2,jr)-r(2,ir)) / w
        dr(3) = (r(3,jr)-r(3,ir)) / w
        rr = SQRT(dr(1)**2+dr(2)**2+dr(3)**2)
        w1 = 1.d0
        call CANSO(w1,dr,s0)
        do jlm = 1,nlm
          do ilm = 1,nlm
            s(ilm+irl0,jlm+jrl0) = s0(ilm+(jlm-1)*nlm)
          end do
        end do
      end do
    end do
  end subroutine streze

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !>  Solves for the screened structure constants
  !
  !> @param[in] r
  !> @param[in] nr
  !> @param[in] nlm
  !> @param[in] nrl
  !> @param[inout] s
  !> @param[inout] a
  !> @param[in] na
  !> @param[in] q
  !> @param[inout] bet
  !> @param[inout] wk
  !> @param[inout] sbar
  !> @param[in] iclus
  !> @param[in] r2
  !---------------------------------------------------------------------------
  subroutine shldch(r,nr,nlm,nrl,s,a,na,q,bet,wk,sbar,iclus,r2)
    implicit none
    !parameter for cutoff of sbar construction
    integer, parameter :: ncut = 9
    ! Input                   
    integer, intent(in) :: iclus,na,nlm,nr,nrl
    real(rp), intent(in) ::r2
    real(rp), dimension(3), intent(in) :: q
    real(rp), dimension(3,nr), intent(in) :: r
    ! Output
    real(rp), dimension(na), intent(inout) :: a
    real(rp), dimension(nrl), intent(inout) :: bet,wk
    real(rp), dimension(nrl,nrl), intent(inout) :: s
    real(rp), dimension(nlm,nlm,nr), intent(inout) :: sbar
    ! Local variables    
    integer :: i,ia,ilm,ir,irl,irl0,isb,j,jlm,jsb,l,l2,lmax,m,ndef,hitc, info
    real(rp), dimension(:,:), allocatable :: s_temp
    ! External Calls   
    !external chlr2f, chlr2s
    ! External Functions  .
    !integer, external :: LL
    ndef = 0 
    lmax = LL(nlm)
    allocate(s_temp(nrl,nrl))
    write (17,10000) lmax, q
    irl = 0
    do ir = 1,nr
      do l = 0,lmax
        do m = 1,2*l+1
          irl = irl + 1
          bet(irl) = 1.d0 / q(l+1)
        end do
      end do
    end do
    s_temp=s
    do i=1,nrl
       s_temp(i,i)=s_temp(i,i)+bet(i)
    end do
!    ia = 0
!    do i = 1,nrl
!      do j = 1,i
!        ia = ia + 1
!        a(ia) = s(i,j)
!        if (i == j) then
!          a(ia) = a(ia) + bet(i)
!        end if
!      end do
!    end do
!    call chlr2f(a,na,wk,nrl,ndef)
    call DPOTRF('U',nrl,s_temp,nrl,info)
    write (17,10001) ndef
!    call chlr2s(a,na,s,nrl,nlm)
    call DPOTRS( 'U', nrl, nlm, s_temp, nrl, s, nrl, INFO )
    deallocate(s_temp)
    do ilm = 1,nlm
      do irl = 1,nrl
        s(irl,ilm) = -bet(irl)*s(irl,ilm)
      end do
    end do
    ! --------------------------------
    ir = 0
    hitc=0
    !print *,' nr = ', nr
    do ir = 1,nr
       if(abs(R(1,IR)**2+R(2,IR)**2+R(3,IR)**2 - &
            R(1,1)**2+R(2,1)**2+R(3,1)**2)<=(R2/ncut) ) then
          write (16,10002) r(1,ir), r(2,ir), r(3,ir), iclus
          hitc=hitc+1
          irl0 = (ir-1) * nlm
          do ilm = 1,nlm
             do jlm = 1,nlm
                sbar(ilm,jlm,hitc) = 2. * s(ilm+irl0,jlm)
             end do
          end do
          !!! Changed for 'symmetrizing' Sbar
          !! NOT NEEDED SINCE WITH LARGER CUTOFF
          !! SBAR IS CALCULATED PROPERLY
          ! do l = 2,9
          !    l2 = l - 1
          !    do j = 1,l2
          !       sbar(l,j,ir) = sbar(j,l,ir)
          !    end do
          ! end do
          ! do l = 1,3
          !    sbar(l+1,1,ir) = -sbar(l+1,1,ir)
          ! end do
          ! do l = 5,9
          !    do j = 2,4
          !       sbar(l,j,ir) = -sbar(l,j,ir)
          !    end do
          ! end do
          !
          do isb = 1,nlm
             write (13) (sbar(isb,jsb,hitc), jsb = 1,nlm)
          end do
          ! 208 FORMAT(5F12.6:/(12X,4F12.6))
          do isb = 1,nlm
             write (16,10003) (sbar(isb,jsb,hitc), jsb = 1,nlm)
          end do
       end if
    end do
    !print *,' hitc = ', hitc
    return
    !
    ! ... Format Declarations ...
    !
    10000 format (" LMAX=",i2,"   Q=",3f10.6)
    10001 format (" NDEF=",i10)
    10002 format (" VECTOR=",3f12.6,"   ICLUS=",i5)
    10003 format (9f10.4)
  end subroutine shldch

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Calculation of the unscreened structure constants in real 
  !> space using the slater_koster table.
  !
  !> @param[in] w
  !> @param[in] dr
  !> @param[inout] sc
  !---------------------------------------------------------------------------
  subroutine canso(w,dr,sc)
    implicit none
    ! Input                     
    real(rp), intent(in) :: w
    real(rp), dimension(3), intent(in) :: dr
    ! Output
    real(rp), dimension(9,9), intent(out) :: sc
    ! Local variables   
    integer :: i,j,l,ll
    real(rp) :: el,el2,elem,elen,em,em2,emen,en,en2,r1,r2,r3,rr,s2,s3,s4,s5, &
                    sbyr,sq3,sq5
    integer, dimension(9) :: ip
    real(rp), dimension(9,9) :: s
    ! Intrinsic Functions   
    intrinsic SQRT
    !.. Data Declarations ..
    ! original and correct
    !     S=1,X=2,Y=3,Z=4,XY=5,YZ=6,ZX=7,X**2-Y**2=8,3Z*Z-R*R=9
    data ip/1,2,3,4,5,6,7,8,9/
    !data ip/1,2,3,4,5,7,6,8,9/
    !data ip/1,2,3,4,5,6,7,9,8/
    ! testing (old convention)
    !data ip /  1,  4,2,3,   5,6,8,9,7  /
    ! testing reversing d-orbitals (bstr convention)
    !data ip/1,2,3,4,9,8,7,6,5/
    !     S=1,X=2,Y=3,Z=4,XY=5,YZ=6,ZX=7,X**2-Y**2=8,3Z*Z-R*R=9
    r1 = dr(1)
    r2 = dr(2)
    r3 = dr(3)
    rr = SQRT(r1*r1+r2*r2+r3*r3)
    do i = 1,9
      do j = 1,9
        sc(i,j) = 0.0d0
      end do
    end do
    if (rr/w <= 0.30d0) return
    sbyr = w / rr
    s2 = sbyr * sbyr
    s3 = s2 * sbyr
    s4 = s3 * sbyr
    s5 = s4 * sbyr
    sq3 = SQRT(3.d0)
    sq5 = SQRT(5.d0)
    el = r1 / rr
    em = r2 / rr
    en = r3 / rr
    el2 = el * el
    em2 = em * em
    en2 = en * en
    elem = el * em
    elen = el * en
    emen = em * en
    !
    sc(1,1) = -2.0d0*sbyr
    sc(1,2) = el * s2 * 2.d0 * sq3
    sc(1,3) = em * s2 * 2.d0 * sq3
    sc(1,4) = en * s2 * 2.d0 * sq3
    sc(1,5) = -2.d0*sq3*sq5*elem*s3
    sc(1,6) = -2.d0*sq3*sq5*emen*s3
    sc(1,7) = -2.d0*sq3*sq5*elen*s3
    sc(1,8) = -sq3*sq5*s3*(el2-em2)
    sc(1,9) = sq5 * s3 * (1.d0-3.d0*en2)
    !-----------------------------------------------------------------------
    sc(2,2) = (3.0d0*el2-1.0) * 6.d0 * s3
    sc(2,3) = 18.0d0 * s3 * elem
    sc(2,4) = 18.0d0 * s3 * elen
    sc(2,5) = 6.d0 * sq5 * s4 * em * (1.0d0-5.0d0*el2)
    sc(2,6) = -30.0d0*sq5*s4*elem*en
    sc(2,7) = 6.d0 * sq5 * s4 * en * (1.0d0-5.0d0*el2)
    sc(2,8) = 6.0d0 * sq5 * s4 * el * (1.0d0-2.5d0*el2+2.5d0*em2)
    sc(2,9) = 3.d0 * sq3 * sq5 * s4 * el * (1.0d0-5.0d0*en2)
    !-----------------------------------------------------------------------
    sc(3,3) = 6.d0 * s3 * (3.0d0*em2-1.d0)
    sc(3,4) = 18.d0 * s3 * emen
    sc(3,5) = 6.0d0 * sq5 * s4 * el * (1.d0-5.d0*em2)
    sc(3,6) = 6.d0 * sq5 * s4 * en * (1.d0-5.d0*em2)
    sc(3,7) = sc(2,6)
    sc(3,8) = -6.d0*sq5*s4*em*(1.d0-2.5d0*em2+2.5d0*el2)
    sc(3,9) = 3.d0 * sq3 * sq5 * s4 * em * (1.d0-5.d0*en2)
    !----------------------------------------------------------------------
    sc(4,4) = 6.d0 * s3 * (3.d0*en2-1.d0)
    sc(4,5) = sc(2,6)
    sc(4,6) = 6.d0 * sq5 * s4 * em * (1.d0-5.d0*en2)
    sc(4,7) = 6.d0 * sq5 * s4 * el * (1.d0-5.d0*en2)
    sc(4,8) = -15.d0*sq5*s4*en*(el*el-em2)
    sc(4,9) = 3.d0 * sq3 * sq5 * s4 * en * (3.0d0-5.0d0*en2)
    !-----------------------------------------------------------------------
    sc(5,5) = 10.d0 * s5 * (-35.d0*el2*em2-5.d0*en2+4.d0)
    sc(5,6) = -50.d0*s5*elen*(7.d0*em2-1.d0)
    sc(5,7) = -50.d0*s5*emen*(7.d0*el2-1.d0)
    sc(5,8) = -175.d0*s5*elem*(el2-em2)
    sc(5,9) = -25.d0*sq3*s5*elem*(7.d0*en2-1.d0)
    !----------------------------------------------------------------------
    sc(6,6) = 10.d0 * s5 * (-35.d0*em2*en2-5.d0*el2+4.d0)
    sc(6,7) = -50.d0*s5*elem*(7.d0*en2-1.d0)
    sc(6,8) = 50.d0 * s5 * emen * (3.5d0*em2-3.5d0*el2-1.d0)
    sc(6,9) = -25.d0*sq3*s5*emen*(7.d0*en2-3.d0)
    !-----------------------------------------------------------------------
    sc(7,7) = 10.d0 * s5 * (-35.d0*el2*en2-5.d0*em2+4.d0)
    sc(7,8) = -50.d0*s5*elen*(3.5d0*el2-3.5d0*em2-1.0d0)
    sc(7,9) = -25.d0*sq3*s5*elen*(7.d0*en2-3.d0)
    !-----------------------------------------------------------------------
    sc(8,8) = 10.d0 * s5 * (-8.75d0*(el2-em2)**2-5.d0*en2+4.d0)
    sc(8,9) = -12.5d0*sq3*s5*(7.d0*en2-1.d0)*(el2-em2)
    !-----------------------------------------------------------------------
    sc(9,9) = -7.5d0*s5*(35.d0*en2*en2-30.d0*en2+3.d0)
    !-----------------------------------------------------------------------
    do l = 2,9
      ll = l - 1
      do j = 1,ll
        sc(l,j) = sc(j,l)
      end do
    end do
    do l = 1,3
      sc(l+1,1) = -sc(l+1,1)
    end do
    do l = 5,9
      do j = 2,4
        sc(l,j) = -sc(l,j)
      end do
    end do
    ! ------ THIS PART CHANGES THE YLM ORDER AND MULTIPLIES BY -0.5 ----
    do i = 1,9
      do j = 1,9
        s(ip(j),ip(i)) = -0.5d0*sc(j,i)
      end do
    end do
    do i = 1,9
      do j = 1,9
        sc(j,i) = s(j,i)
      end do
    end do
  end subroutine canso

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Cholesky factorization of (N X N) REAL*8 matrix c.
  !> returns as soon as matrix becomes indefinite, rest is unchanged
  !> 
  !
  !> @param[inout] C 
  !> @param[in] NA
  !> @param[inout] W Workspace of dimension N
  !> @param[in] N
  !> @param[out] NDEF
  !---------------------------------------------------------------------------
  subroutine chlr2f(c,na,w,n,ndef)
    implicit none
    ! Input                       
    integer, intent(in) :: n,na
    ! Output
    integer, intent(out) :: ndef
    real(rp), dimension(n), intent(inout) :: w
    real(rp), dimension(na), intent(inout) :: c
    ! Local variables    
    integer :: i,ic0,j,jc0,k
    real(rp) :: csum
    ! Intrinsic Functions    
    intrinsic SQRT
 
    write(17,*)'chol',na,n
    ic0 = 0
    do i = 1,n
      jc0 = 0
      do j = 1,i-1
        csum = c(ic0+j)
        do k = 1,j-1
          csum = csum - w(k)*c(jc0+k)
        end do
        jc0 = jc0 + j
        w(j) = csum / c(jc0)
      end do
      csum = c(ic0+i)
      do k = 1,i-1
        csum = csum - w(k)*w(k)
      end do
      if (csum <= 0.d0) then
        ndef = i - 1
        goto 1000
      else
        do j = 1,i-1
          c(ic0+j) = w(j)
        end do
        !
        ic0 = ic0 + i
        c(ic0) = SQRT(csum)
      end if
    end do
    ndef = n
    1000 if (ndef < n) then
           write (17,10000) n, ndef
         end if
    return
    10000 format (" CHLR2F:    N=",i4,"    NDEF=",i4)
  end subroutine chlr2f

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Back substitution after cholesky factorization
  !> solution overwrites RHS in V
  !
  !> @param[in] C
  !> @param[in] NA
  !> @param[inout] V
  !> @param[in] N
  !> @param[in] M Number of RHS's
  !---------------------------------------------------------------------------
  subroutine chlr2s(c,na,v,n,m)
    implicit none
    ! Input                    
    integer, intent(in) :: m,n,na
    real(rp), dimension(na), intent(in) :: c
    ! Output
    real(rp), dimension(n,m), intent(inout) :: v
    ! Local Variables     
    integer :: i,ic0,ind1,ind2,k,mm
    real(rp) :: csum

    do mm = 1,m
      ic0 = 0
      do i = 1,n
        csum = v(i,mm)
        do k = 1,i-1
          csum = csum - v(k,mm)*c(ic0+k)
        end do
        ic0 = ic0 + i
        v(i,mm) = csum / c(ic0)
      end do
      !
      do i = n,1,-1
        ind2 = (i*(i+1)) / 2
        csum = v(i,mm)
        ind1 = ind2 + i
        do k = i+1,n
          csum = csum - v(k,mm)*c(ind1)
          ind1 = ind1 + k
        end do
        v(i,mm) = csum / c(ind2)
      end do
    end do
  end subroutine chlr2s

  function LL(ilm)
    implicit none
    ! Input                       
    integer, intent(in) :: ilm
    ! Function Declaration 
    integer :: LL
    ! Local variables    
    integer, dimension(100) :: lla
    ! Data Declarations   
    data lla/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9/

    LL = lla(ilm)
  end function LL

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Reorders map according to eq vectors
  !
  !> @param[in] crd
  !> @param[in] no
  !> @param[in] iu
  !> @param[inout] nn
  !> @param[in] nat
  !> @param[in] ntot
  !> @param[in] nomx
  !> @param[in] ndi
  !> @param[in] nnmx
  !> @param[inout] set
  !> @param[inout] idnn
  !> @param[out] ret
  !---------------------------------------------------------------------------
  subroutine remd(crd,no,iu,nn,nat,ntot,nomx,ndi,nnmx,set,idnn,ret)
    implicit none
    ! Inputs                                 
    integer, intent(in) :: nat,ndi,nnmx,nomx,ntot
    integer, dimension(ndi), intent(in) :: no
    real(rp), dimension(3,ndi), intent(in) :: crd
    integer, dimension(nomx), intent(in) :: iu
    ! Output
    integer, dimension(nnmx), intent(inout) :: idnn
    integer, dimension(ndi,nnmx), intent(inout) :: nn
    real(rp), dimension(3), intent(out) :: ret
    real(rp), dimension(3,nomx,nnmx), intent(inout) :: set
    ! Local variables    
    integer :: i,ii,iii,imax,inn,ino,j,jj,jnn,jsz,k,la,lk,lm,m,n
    real(rp) :: a1,a2,a3,aaa,eps
    !-BUILDS VECTORS SET(3,NOMX,NNMX) CONNECTING NEIGHBORS OF EACH TYPE NO-
    do i = 1,ntot
      la = iu(i)
      jsz = nn(la,1)
      do j = 2,jsz
        jj = nn(la,j)
        do m = 1,3
          set(m,i,j) = crd(m,la) - crd(m,jj)
        end do
      end do
    end do
    do i = 1,nat
      n = no(i)
      do lk = 1,ntot
        ino = iu(lk)
        lm = no(ino)
        if (n == lm) goto 1000
      end do
      write (17,10003)
      1000 imax = nn(ino,1)
      do iii = 1,imax
        idnn(iii) = 0
      end do
      jsz = nn(i,1)
      do j = 2,jsz
        jj = nn(i,j)
        do m = 1,3
          ret(m) = crd(m,i) - crd(m,jj)
        end do
        !----------FINDS EQUIVALENT VECTOR------------------------
        eps = .0001
        do ii = 2,imax
          a1 = ret(1) - set(1,n,ii)
          a2 = ret(2) - set(2,n,ii)
          a3 = ret(3) - set(3,n,ii)
          aaa = a1**2 + a2**2 + a3**2
          if (aaa < eps) goto 1100
        end do
        goto 1200
        1100 k = ii
        idnn(k) = jj
      end do
      nn(i,1) = imax
      do j = 2,imax
        nn(i,j) = idnn(j)
      end do
    end do
    do inn = 1,nat
      write (12) (nn(inn,jnn), jnn = 1,nn(inn,1))
    end do
    return
    1200 write (17,10002)
    print *,'atom:',i
    stop
    !  9  WRITE(6,223)I,SET(1,I,J),SET(2,I,J),SET(3,I,J)
    10000 format (i5,3f9.4)
    !--REORDERS NEIGHBORS FOR I LARGER THAN NTOT ACCORDING TO TYPICAL NO--
    10001 format (3i5)
    10002 format (" VECTOR   NOT FOUND ")
    10003 format (" TYPE NO  NOT FOUND ")
  end subroutine remd

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> TODO
  !
  !> @param[in] IM
  !> @param[in] IZP
  !> @param[in] NN
  !> @param[in] NO
  !> @param[in] ND
  !> @param[in] NM
  !> @param[in] NTOT
  !---------------------------------------------------------------------------
  subroutine outmap(IM,IZP,NN,NO,ND,NM,NTOT)
    implicit none
    ! Input                  
    integer, intent(in) :: IM,ND,NM,NTOT
    integer, dimension(ND), intent(in) :: NO
    integer, dimension(ND) :: IZP
    integer, dimension(ND,NM), intent(in) :: NN
    ! Local variables   
    integer :: I,ID,IDM,J,K
    ! Intrinsic functions     
    intrinsic MIN

    write (IM,10000)
    do I = 1,NTOT
!!    K = NO(I)
      K = IZP(I)
      ID = NN(I,1)
!!    IDM = MIN(ID,21)
      IDM = MIN(ID,16)
      write (IM,10001) I, K, ID, (NN(I,J), J = 2,IDM)
      if (IDM /= NN(I,1)) then
        !     ID=NN(I,1)
        write (IM,10002) (NN(I,J), J = 17,ID)
      end if
    end do
    return

    10000 format ( &
          " NEAREST NEIGHBOUR MAP"/,"      ATOM   TYPE  CONNECTIVITY",5x, &
          "NEIGHBOURS")
    10001 format (1x,3i4,2x,20i5)
    10002 format (29x,16i5)
  end subroutine outmap

  integer function mapa(I,J,R2,DD,CT)
    implicit none
    ! Input                            
    integer, intent(in) :: I,J
    real(rp), intent(in) :: R2
    real(rp) :: DD
    real(rp), dimension(50), intent(in) :: CT
    ! Local variables   
    real(rp) :: CTM,CTSM

    CTM = (CT(1)+CT(1)) / 2.
    CTSM = CTM**2
    if (R2 >= CTSM) then
      MAPA = 0
    else
      MAPA = 1
    end if
  end function mapa


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Recursion library
  !
  !> @param[inout] ct
  !> @param[in] crd
  !> @param[in] ndim
  !> @param[in] nat
  !> @param[in] izp
  !> @param[inout] nn
  !> @param[in] nd
  !> @param[inout] nm
  !> @param[in] ngbr
  !> @param[in] ntot
  !---------------------------------------------------------------------------
  subroutine nncal(ct,crd,ndim,nat,izp,nn,nd,nm,ngbr,ntot)
    implicit none
    ! Input
    integer, intent(in) :: NAT,ND,NDIM,NTOT
    integer, dimension(NAT), intent(in) :: IZP
    real(rp), dimension(NDIM,NAT), intent(in) :: CRD
    ! Output
    integer, intent(inout) :: NM
    integer, dimension(ND,NM), intent(inout) :: NN
    real(rp), dimension(50), intent(inout) :: CT
    ! External function
    integer, external :: NGBR
    ! Intrinsic function
    intrinsic ABS, MAX
    ! Local variables
    integer :: I,IADD,ID,II,IIP,ILJ,J,JJP,L,NNMAX
    real(rp) :: R2
    real(rp), dimension(3) :: DDUM
    real(rp), dimension(NM) :: DUM
    
    NNMAX = 0
    IADD = 1
    if (IZP(1) < 0) then
      NN(1,1) = 1
    end if
    do II = 1,NAT
      if (IZP(II) < 0) goto 1000
    end do
    NN(1,1) = 1
    do I = 1,NAT
      do J = 2,NM
        NN(I,J) = 0
      end do
    end do
    IADD = 0
    II = 2
    1000 if (II <= 1) then
          II = 2
        end if
    do I = II,NAT
      if (IZP(I)*IADD <= 0) then
        NN(I,1) = 1
        IIP = IZP(I)
        IIP = ABS(IIP)
        ILJ = I - 1
        do J = 1,ILJ
          JJP = IZP(J)
          JJP = ABS(JJP)
          R2 = 0.0
          do L = 1,3
            DDUM(L) = CRD(L,I) - CRD(L,J)
            R2 = R2 + DDUM(L)*DDUM(L)
          end do
          ID = NGBR(IIP,JJP,R2,DUM,CT)
  !       if (ID /= 0 .and. ( IZP(I) > NTOT  .or. IZP(J) > NTOT) ) then
          if (ID /= 0) then
            ID = NN(I,1) + 1
            NN(I,1) = ID
  !         NN(I,ID) = IZP(J)
            NN(I,ID) = J
            NNMAX = MAX(NNMAX,ID)
            ID = NN(J,1) + 1
            NN(J,1) = ID
  !         NN(J,ID) = IZP(I)
            NN(J,ID) = I
            NNMAX = MAX(NNMAX,ID)
            if (NNMAX > NM) then
              write (6,10001)
              write (6,10002) I
              write(6,*) NNMAX,ID,NM
              stop
            end if
          end if
        end do
      end if
    end do
    nm = nnmax
    return

    10000 format ("FROM NNCAL ")
    10001 format (" TOO MANY NEIGHBOURS")
    10002 format (" NEIGHBOUR MAP AS FAR AS",i6,"TH SITE")
  end subroutine
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Read CR of clust file
  !
  !> @param[in] alat
  !> @param[in] nndim
  !> @param[inout] cr
  !> @param[inout] iz
  !> @param[inout] n
  !> @param[in] ip Unit of opened file 
  !---------------------------------------------------------------------------
  subroutine leia(alat,nndim,cr,iz,n,ip)
    implicit none
    ! Input
    integer, intent(in) :: ip,nndim
    real(rp), intent(in) :: alat
    ! Output
    integer, dimension(nndim+10), intent(inout) :: iz,n
    real(rp), dimension(3,nndim+10), intent(inout) :: cr
    ! Local variables    
    integer :: i,j

    read (IP,*)    
    do I = 1,nndim,2
      read (IP,*) &
          CR(1,I), CR(2,I), CR(3,I), IZ(I), N(I), CR(1,I+1), CR(2,I+1), &
          CR(3,I+1), IZ(I+1), N(I+1)
    end do
    do I = 1,nndim
      do J = 1,3
        CR(J,I) = ALAT * CR(J,I)
      end do
    end do
    return
  end subroutine
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Subrotina que ordena um determinado vetor 'M' em ordem decrescente
  !>     segundo o metodo da bolha (BUBBLE)
  !> ...M(1) >= M(2)...  ; onde:
  !>     NL,ndim   -  input
  !>     M         -  input/output (na saida M contem os elementos dispostos
  !>     em ordem decrescente em modulo.
  !>
  !> @param[in] nl
  !> @param[in] ndim
  !> @param[inout] m
  !> @param[in] nd
  !> @param[in] nt
  !---------------------------------------------------------------------------
  subroutine bubble(nl,ndim,m,nd,nt)
    implicit none
    ! Inputs
    integer, intent(in) :: ndim, nl, nd, nt
    ! Output
    real(rp), dimension(ndim,nd), intent(inout) :: m
    ! Local variables
    integer :: ind, inic, j, k
    real(rp) :: fim, z

    IND = 1
    INIC = 2
    FIM = NL
    do while ((IND==1) .and. (INIC<=FIM))
      IND = 0
      do J = INT(FIM),INIC,-1
        if (M(J,NT) < M(J-1,NT)) then
        do K = 1,ND
          Z = M(J,K)
          M(J,K) = M(J-1,K)
          M(J-1,K) = Z
        end do
        end if
      end do
      FIM = FIM - 1
      do J = INIC,INT(FIM)
        if (M(J+1,NT) < M(J,NT)) then
        do K = 1,ND
          Z = M(J+1,K)
          M(J+1,K) = M(J,K)
          M(J,K) = Z
        end do
          IND = 1
        end if
      end do
      INIC = INIC + 1
    end do
  end subroutine bubble

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> CUT:subprogram of the builds package
  !> Cuts a sphere of square raius RS around every atom in the unit
  !> cell. Best for cells with large numbers of atoms...  hould
  !> Carefull with monoatomic systems!!! In this case should cut
  !> around several sites, since spherical clusters should be avoided.
  !>
  !> @param[in] i
  !> @param[in] l
  !> @param[in] ndim
  !> @param[in] crd
  !> @param[inout] cr
  !> @param[in] izp
  !> @param[inout] iz
  !> @param[inout] num
  !> @param[in] no
  !> @param[in] rs
  !> @param[in] ii
  !---------------------------------------------------------------------------
  subroutine cut(i,l,ndim,crd,cr,izp,iz,num,no,rs,ii)
    ! Inputs
    real(rp),intent(in) :: rs
    integer,intent(in) :: i, l, ndim
    real(rp),dimension(3,ndim),intent(in) :: crd
    integer, dimension(ndim),intent(in) :: izp, no
    ! Output
    integer,intent(out) :: ii           
    real(rp),dimension(3,ndim),intent(out) :: cr     
    integer,dimension(ndim),intent(out) :: iz, num       
    ! Local variables
    real(rp) :: r2
    real(rp), dimension(3) :: dum
    integer :: na, j

    do na=1,l
      r2=0.0d0
      do j=1,3
        dum(j)=(crd(j,na)-crd(j,i))**2
        r2=r2+dum(j)
      end do
      if(r2.le.rs)then 
        ii=ii+1
        cr(1,ii)=crd(1,i)
        cr(2,ii)=crd(2,i)
        cr(3,ii)=crd(3,i)
        iz(ii)=izp(i)
        num(ii)=no(i)
        return
      end if
    end do
    return
  end subroutine cut

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Check if object is well fulfilled
  !---------------------------------------------------------------------------
  subroutine check_all(this)
    implicit none
    class(lattice) :: this

    if(this%crystal_sym /= 'bcc' &
    .and. this%crystal_sym /= 'fcc' &
    .and. this%crystal_sym /= 'hcp'  & 
    .and. this%crystal_sym /= 'fcc2' &
    .and. this%crystal_sym /= 'fcc3') then
      write(error_unit,'("[",A,":",I0,"]: lattice%crystal_sym must be one of: ''bcc'', ''fcc'', ''hcp''")') __FILE__,__LINE__
     error stop
    end if
  end subroutine check_all

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Calculate reduced_acr and reduced_nbas
  !---------------------------------------------------------------------------
  subroutine calculate_nbas(this)
    implicit none
    class(lattice) :: this
    integer :: size_iz
    integer, dimension(:) :: atype(this%nbas), amount(this%nbas)
    integer :: i,j

    amount(:) = 0
    atype(:) = 0

    j = 1

    amount(j) = 1
    atype(j) = this%iz(j)

    do i=2,this%nbas
      if(this%iz(i) .eq. this%iz(i-1)) then
        amount(j) = amount(j) + 1
      else
        j = j + 1
        atype(j) = this%iz(i)
        amount(j) = 1
      endif
    enddo

    this%reduced_nbas = 0

    do i=1,size(amount)
      if(amount(i) .eq. 0) exit
      this%reduced_nbas = this%reduced_nbas + 1
    enddo

    if(allocated(this%reduced_acr)) deallocate(this%reduced_acr)
    allocate(this%reduced_acr(this%reduced_nbas))

    do i=1,size(this%reduced_acr)
      this%reduced_acr(i) = atype(i)
    enddo

  end subroutine calculate_nbas

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Print class members values in namelist format 
  !>
  !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
  !> @param[in] unit File unit used to write namelist
  !> @param[in] file File name used to write namelist
  !---------------------------------------------------------------------------
  subroutine print_state_full(this,unit,file)
    implicit none
    class(lattice), intent(in) :: this

    integer,intent(in),optional :: unit
    character(len=*),intent(in),optional :: file
    integer :: newunit

    real(rp) :: zmin, zmax, zstep, wav, vol, rc, r2, celldm, alat
    integer :: reduced_nbas, ntype, ntot, nrec, nmax, nlay, ndim, npe, nclu, nbulk_bulk, nbulk, nbas, kk, dx, dy, dz, dw
    character(len=4) :: crystal_sym
    character(len=10) :: surftype
    real(rp), dimension(3,3) :: a
    real(rp), dimension(:), allocatable :: z, ct
    integer, dimension(:), allocatable :: reduced_acr, num, no, izpsurf, izsurf, nosurf, izpo, izp, iz, iu, irec, ib
    real(rp), dimension(:,:), allocatable :: primcell, inclu, crsurf, crd, cr, acr

    namelist /lattice/ a, z, ct, primcell, inclu, crsurf, crd, cr, &
    acr, zmin, zmax, zstep, wav, vol, rc, r2, celldm, alat, reduced_acr,&
    num, no, izpsurf, izsurf, nosurf, izpo, izp, iz, iu, irec, ib, &
    reduced_nbas, ntype, ntot, nrec, nmax, nlay, ndim, npe, nclu, &
    nbulk_bulk, nbulk, nbas, kk, dx, dy, dz, dw, crystal_sym, surftype 


    ! scalar

    zmin = this%zmin
    zmax = this%zmax
    zstep = this%zstep
    wav = this%wav
    vol = this%vol
    rc = this%rc
    r2 = this%r2
    celldm = this%celldm
    alat = this%alat
    reduced_nbas = this%reduced_nbas
    ntype = this%ntype
    ntot = this%ntot
    nrec = this%nrec
    nmax = this%nmax
    nlay = this%nlay
    ndim = this%ndim
    npe = this%npe
    nclu = this%nclu
    nbulk_bulk = this%nbulk_bulk
    nbulk = this%nbulk
    nbas = this%nbas
    kk = this%kk
    dx = this%dx
    dy = this%dy
    dz = this%dz
    dw = this%dw
    crystal_sym = this%crystal_sym
    surftype = this%surftype
    a = this%a

    ! one dimensional allocatables

    if(allocated(this%z)) then
      allocate(z,mold=this%z)
      z = this%z
    else
      allocate(z(0))
    endif
    if(allocated(this%ct)) then
      allocate(ct,mold=this%ct)
      ct = this%ct
    else
      allocate(ct(0))
    endif
    if(allocated(this%reduced_acr)) then
      allocate(reduced_acr,mold=this%reduced_acr)
      reduced_acr = this%reduced_acr
    else
      allocate(reduced_acr(0))
    endif
    if(allocated(this%num)) then
      allocate(num,mold=this%num)
      num = this%num
    else
      allocate(num(0))
    endif
    if(allocated(this%no)) then
      allocate(no,mold=this%no)
      no = this%no
    else
      allocate(no(0))
    endif
    if(allocated(this%izpsurf)) then
      allocate(izpsurf,mold=this%izpsurf)
      izpsurf = this%izpsurf
    else
      allocate(izpsurf(0))
    endif
    if(allocated(this%izsurf)) then
      allocate(izsurf,mold=this%izsurf)
      izsurf = this%izsurf
    else
      allocate(izsurf(0))
    endif
    if(allocated(this%nosurf)) then
      allocate(nosurf,mold=this%nosurf)
      nosurf = this%nosurf
    else
      allocate(nosurf(0))
    endif
    if(allocated(this%izpo)) then
      allocate(izpo,mold=this%izpo)
      izpo = this%izpo
    else
      allocate(izpo(0))
    endif
    if(allocated(this%izp)) then
      allocate(izp,mold=this%izp)
      izp = this%izp
    else
      allocate(izp(0))
    endif
    if(allocated(this%iz)) then
      allocate(iz,mold=this%iz)
      iz = this%iz
    else
      allocate(iz(0))
    endif
    if(allocated(this%iu)) then
      allocate(iu,mold=this%iu)
      iu = this%iu
    else
      allocate(iu(0))
    endif
    if(allocated(this%irec)) then
      allocate(irec,mold=this%irec)
      irec = this%irec
    else
      allocate(irec(0))
    endif
    if(allocated(this%ib)) then
      allocate(ib,mold=this%ib)
      ib = this%ib
    else
      allocate(ib(0))
    endif

    ! two dimensional allocatables

    if(allocated(this%primcell)) then
      allocate(primcell,mold=this%primcell)
      primcell = this%primcell
    else
      allocate(primcell(0,0))
    endif
    if(allocated(this%inclu)) then
      allocate(inclu,mold=this%inclu)
      inclu = this%inclu
    else
      allocate(inclu(0,0))
    endif
    if(allocated(this%crsurf)) then
      allocate(crsurf,mold=this%crsurf)
      crsurf = this%crsurf
    else
      allocate(crsurf(0,0))
    endif
    if(allocated(this%crd)) then
      allocate(crd,mold=this%crd)
      crd = this%crd
    else
      allocate(crd(0,0))
    endif
    if(allocated(this%cr)) then
      allocate(cr,mold=this%cr)
      cr = this%cr
    else
      allocate(cr(0,0))
    endif
    if(allocated(this%acr)) then
      allocate(acr,mold=this%acr)
      acr = this%acr
    else
      allocate(acr(0,0))
    endif

    if(present(unit) .and. present(file)) then
      write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
      error stop
    else if(present(unit)) then
        write(unit,nml=lattice)
    else if(present(file)) then
        open(unit=newunit,file=file)
        write(newunit,nml=lattice)
        close(newunit)
    else
        write(*,nml=lattice)
    endif
    
  end subroutine print_state_full

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
    class(lattice), intent(in) :: this

    integer,intent(in),optional :: unit
    character(len=*),intent(in),optional :: file
    integer :: newunit

    real(rp) :: wav, rc, celldm, alat
    integer :: nlay, ndim, npe, nclu
    character(len=4) :: crystal_sym
    character(len=10) :: surftype
    real(rp), dimension(3,3) :: a
    integer, dimension(:), allocatable :: no, izp
    real(rp), dimension(:,:), allocatable :: inclu, crd

    namelist /lattice/ crystal_sym, rc, a, izp, no, crd, ndim, npe,& ! Bulk parameters
                       surftype, nlay,&
                       inclu, nclu, & ! Impurity varables 
                       alat, celldm, wav        ! General parameters 

    ! scalar

    wav = this%wav
    rc = this%rc
    celldm = this%celldm
    alat = this%alat
    nlay = this%nlay
    ndim = this%ndim
    npe = this%npe
    nclu = this%nclu
    crystal_sym = this%crystal_sym
    surftype = this%surftype
    a = this%a

    ! ! one dimensional allocatables

    if(allocated(this%no)) then
      allocate(no,mold=this%no)
      no = this%no
    else
      allocate(no(0))
    endif
    if(allocated(this%izp)) then
      allocate(izp,mold=this%izp)
      izp = this%izp
    else
      allocate(izp(0))
    endif

    ! ! two dimensional allocatables

    if(allocated(this%inclu)) then
      allocate(inclu,mold=this%inclu)
      inclu = this%inclu
    else
      allocate(inclu(0,0))
    endif
    if(allocated(this%crd)) then
      allocate(crd,mold=this%crd)
      crd = this%crd
    else
      allocate(crd(0,0))
    endif

    if(present(unit) .and. present(file)) then
      write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
      error stop
    else if(present(unit)) then
        write(unit,nml=lattice)
    else if(present(file)) then
        open(unit=newunit,file=file,action='write')
        write(newunit,nml=lattice)
        close(newunit)
    else
        write(*,nml=lattice)
    endif
    
  end subroutine print_state

end module lattice_mod