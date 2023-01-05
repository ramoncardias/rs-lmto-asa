!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Hamiltonian
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
!> Module to handle procedures related to the Hamiltonian
!------------------------------------------------------------------------------


module hamiltonian_mod
  
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use symbolic_atom_mod
  use element_mod
  use potential_mod
  use lattice_mod
  use charge_mod
  use precision_mod, only: rp 
  use math_mod
  implicit none

  private
 
  !> Module's main procedure
  type, public :: hamiltonian
    !> Spin-orbit coupling Hamiltonian
    complex(rp), dimension(:,:,:), allocatable :: lsham 
    !> Bulk Hamiltonian
    complex(rp), dimension(:,:,:,:), allocatable :: ee
    !> Local Hamiltonian
    complex(rp), dimension(:,:,:,:), allocatable :: hall
    !> Hamiltonian built in chbar_nc (description to be improved)
    complex(rp), dimension(:,:,:,:), allocatable :: hmag 
    !> Hamiltonian built in ham0m_nc (description to be improved
    complex(rp), dimension(:,:,:), allocatable :: hhmag
    !> Charge 
    class(charge), pointer :: charge
    !> Lattice
    class(lattice), pointer :: lattice
  contains
    procedure :: build_lsham
    procedure :: build_bulkham
    procedure :: build_locham
    procedure :: chbar_nc
    procedure :: ham0m_nc
    procedure :: hmfind
    procedure :: restore_to_default
  end type hamiltonian
  
  interface hamiltonian
    procedure :: constructor
  end interface

contains

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] charge_obj Variable that holds charge_mod properties
  !> @return type(hamiltonian)
  !---------------------------------------------------------------------------
  function constructor(charge_obj) result(obj)
    type(hamiltonian) :: obj
    type(charge), target, intent(in) :: charge_obj

    obj%charge => charge_obj
    obj%lattice => charge_obj%lattice

    call obj%restore_to_default()
  end function constructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine destructor(this)
    type(hamiltonian) :: this
    if(allocated(this%lsham)) deallocate(this%lsham)
    if(allocated(this%ee)) deallocate(this%ee)
    if(allocated(this%hmag)) deallocate(this%hmag)
    if(allocated(this%hhmag)) deallocate(this%hhmag)
    if(allocated(this%hall)) deallocate(this%hall)
  end subroutine destructor

  ! Member functions
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Reset all members to default
  !---------------------------------------------------------------------------
  subroutine restore_to_default(this)
    class(hamiltonian) :: this

    allocate(this%lsham(18,18,this%charge%lattice%ntype))
    allocate(this%hhmag(9,9,4),this%hmag(9,9,this%charge%lattice%kk,4))
    allocate(this%ee(18,18,(this%charge%lattice%nn(1,1)+1),this%charge%lattice%kk))
    allocate(this%hall(18,18,(this%charge%lattice%nn(1,1)+1),this%charge%lattice%nmax))

    this%hall(:,:,:,:) = 0.0d0
    this%ee(:,:,:,:) = 0.0d0
    this%lsham(:,:,:) = 0.0d0
    this%hhmag(:,:,:) = 0.0d0
    this%hmag(:,:,:,:) = 0.0d0
  end subroutine restore_to_default

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Build the spin-orbit coupling hamiltonian according to Wu/Freeman PRB 54, 61 (1996)
  !---------------------------------------------------------------------------
  subroutine build_lsham(this)
    class(hamiltonian),intent(inout) :: this
    ! Local variables
    integer :: i, j, k 
    complex(rp) :: sg
    complex(rp), dimension(9,9) :: Lx, Ly, Lz

    !  Getting the angular momentum operators from the math_mod that are in cartesian coordinates
    Lx(:,:) = L_x(:,:) 
    Ly(:,:) = L_y(:,:) 
    Lz(:,:) = L_z(:,:) 

    ! Transforming them into the spherical harmonics coordinates
    call hcpx(Lx,'cart2sph')
    call hcpx(Ly,'cart2sph')  
    call hcpx(Lz,'cart2sph')  
    
    ! Missing the xi of the p and d orbital to be added here

    ! Writing the L.S hamiltonian
    do k=1,this%charge%lattice%ntype
      do i=1,9
        do j=1,9
          this%lsham(j  ,i  ,k)=this%lsham(j  ,i  ,k)+Lz(j,i) ! H11
          this%lsham(j  ,i+9,k)=this%lsham(j  ,i+9,k)+Lx(j,i)-i_unit*Ly(j,i) ! H12
          this%lsham(j+9,i  ,k)=this%lsham(j+9,i  ,k)+Lx(j,i)+i_unit*Ly(j,i) ! H21
          this%lsham(j+9,i+9,k)=this%lsham(j+9,i+9,k)-Lz(j,i) ! H22
        end do
      end do    
    end do

    this%lsham(:,:,:) = 0.09*this%lsham(:,:,:)
  end subroutine build_lsham 

  subroutine build_bulkham(this)
    class(hamiltonian), intent(inout) :: this
    ! Local variables
    integer :: i, j, k, l, m, n, itype, ino, ja, jo, ji, nr, ia
    integer :: ntype

    do ntype=1,this%charge%lattice%ntype
      ia = this%charge%lattice%atlist(ntype) ! Atom number in clust
      ino = this%charge%lattice%num(ia) ! Atom bravais type of ia
      nr = this%charge%lattice%nn(ia,1) ! Number of neighbours considered
      !write(123,*)'bulkham' 
      call this%chbar_nc(ia,nr,ino,ntype)
      do m=1,nr
        do i=1,9
          do j=1,9
            this%ee(j  ,i  ,m,ntype) = this%hmag(j,i,m,4)+this%hmag(j,i,m,3)        ! H0+Hz
            this%ee(j+9,i+9,m,ntype) = this%hmag(j,i,m,4)-this%hmag(j,i,m,3)        ! H0-Hz
            this%ee(j  ,i+9,m,ntype) = this%hmag(j,i,m,1)-i_unit*this%hmag(j,i,m,2) ! Hx-iHy
            this%ee(j+9,i  ,m,ntype) = this%hmag(j,i,m,1)+i_unit*this%hmag(j,i,m,2) ! Hx+iHy
          end do ! end of orbital j loop
        end do ! end of orbital i loop
      !write(128,*) 'm=',m
      !write(128,'(18f10.6)') real(this%ee(:,:,m,ntype))
      end do ! end of neighbour number
    end do ! end of atom type number
  end subroutine build_bulkham

  subroutine build_locham(this)
    class(hamiltonian), intent(inout) :: this
    ! Local variables
    integer :: it, ino, nr, nlim, m, i, j


    do nlim = 1,this%charge%lattice%nmax
      nr = this%charge%lattice%nn(nlim,1) ! Number of neighbours considered
      ino = this%charge%lattice%num(nlim)
      call this%chbar_nc(nlim,nr,ino,nlim)
      do m=1,nr
        do i=1,9
          do j=1,9
            this%hall(j  ,i  ,m,nlim) = this%hmag(j,i,m,4)+this%hmag(j,i,m,3) ! H0+Hz
            this%hall(j+9,i+9,m,nlim) = this%hmag(j,i,m,4)-this%hmag(j,i,m,3) ! H0-Hz
            this%hall(j  ,i+9,m,nlim) = this%hmag(j,i,m,1)-i_unit*this%hmag(j,i,m,2) ! Hx-iHy
            this%hall(j+9,i  ,m,nlim) = this%hmag(j,i,m,1)+i_unit*this%hmag(j,i,m,2) ! Hx+iHy
          end do
        end do
      end do
    end do
  end subroutine build_locham

  subroutine ham0m_nc(this,it,jt,vet,hhh)
    class(hamiltonian), intent(inout) :: this
    ! Input
    integer, intent(in) :: it, jt ! Type of atom i and j
    real(rp), dimension(3), intent(in) :: vet
    real(rp), dimension(9,9), intent(in) :: hhh
    ! Local Variables
    integer :: i, j, ilm, jlm, m
    complex(rp), dimension(3) :: cross
    complex(rp), dimension(9,9) :: hhhc
    complex(rp), dimension(this%charge%lattice%ntype,3) :: momc
    complex(rp) :: dot
    real(rp) :: vv

    this%hhmag(:,:,:) = 0.0d0

    vv = norm2(vet)

    ! Real to complex
    dot = cmplx(dot_product(this%charge%lattice%symbolic_atoms(it)%potential%mom,this%charge%lattice%symbolic_atoms(jt)%potential%mom),kind=kind(0.0d0))
    do i=1,this%charge%lattice%ntype
      do j=1,3
        momc(i,j) = cmplx(this%charge%lattice%symbolic_atoms(i)%potential%mom(j),kind=kind(0.0d0))
      end do
    end do
    cross = cmplx(cross_product(this%charge%lattice%symbolic_atoms(it)%potential%mom,this%charge%lattice%symbolic_atoms(jt)%potential%mom),kind=kind(0.0d0))
    hhhc(:,:) = cmplx(hhh(:,:),kind=kind(0.0d0)) 

    do ilm = 1,9
      do jlm = 1,9
        this%hhmag(ilm,jlm,4) = &
        this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm,jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm) + &
        this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm,jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*dot   
      end do
    end do

!    do ilm=1,9
!      write(123,'(9f10.6)') (real(this%hhmag(ilm,jlm,4)),jlm=1,9)
!    end do

    if(vv<=0.01d0)then
      do ilm = 1,9
        this%hhmag(ilm,ilm,4) = this%hhmag(ilm,ilm,4)+this%charge%lattice%symbolic_atoms(it)%potential%cx0(ilm)
      end do
    end if

    do m=1,3
      do jlm = 1,9
        do ilm = 1,9
          this%hhmag(ilm,jlm,m) = &
          (this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm,jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx0(jlm))*momc(it,m)+&
          (this%charge%lattice%symbolic_atoms(it)%potential%wx0(ilm)*hhhc(ilm,jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm))*momc(jt,m)+&
          i_unit*this%charge%lattice%symbolic_atoms(it)%potential%wx1(ilm)*hhhc(ilm,jlm)*this%charge%lattice%symbolic_atoms(jt)%potential%wx1(jlm)*cross(m) 
        end do
      end do
    end do

    if(vv>0.01d0)return
    do m=1,3
      do ilm=1,9
        this%hhmag(ilm,ilm,m) = this%hhmag(ilm,ilm,m)+this%charge%lattice%symbolic_atoms(it)%potential%cx1(ilm)*momc(it,m)
      end do
    end do
    
!    do m=1,3
!      write(123,*)'m=',m
!      do ilm=1,9
!        write(123,'(9f10.6)') (real(this%hhmag(ilm,jlm,m)),jlm=1,9)
!      end do
!    end do
  end subroutine ham0m_nc

  subroutine chbar_nc(this,ia,nr,ino,ntype)
    class(hamiltonian), intent(inout) :: this
    ! Input
    integer, intent(in) :: ia ! Atom number in clust
    integer, intent(in) :: nr ! Number of neighbours considered
    integer, intent(in) :: ino ! Atom bravais type of ia 
    integer, intent(in) :: ntype ! Atom type 
    ! Local variables
    real(rp) :: r2 
    real(rp), dimension(3,size(this%charge%lattice%cr(1,:))) :: cralat ! Clust position times the lattice constant
    real(rp), dimension(3) :: vet
    real(rp), dimension(9,9) :: hhh
    integer :: i, j, k, l, m, n, it, jt, jj, dummy
    integer :: ni, mdir
    integer :: kk ! Clust size number

    this%hmag(:,:,:,:) = 0.0d0

    r2 = this%charge%lattice%r2
    cralat(:,:) = this%charge%lattice%cr(:,:)*this%charge%lattice%alat
    kk = this%charge%lattice%kk

    call this%charge%lattice%clusba(r2,cralat,ia,kk,kk,dummy)

    !do m=1,nr
    !  print '(9f10.6)', real(this%charge%lattice%sbar(:,:,m,ino))
    !end do
    it = this%charge%lattice%iz(ia)
    do m=1,nr
      jj = this%charge%lattice%nn(ia,m)
      !write(123,*)'ia,ii',ia,m, this%charge%lattice%nn(ia,m)
      if(m==1)then
        jj = ia
      end if
      if(jj/=0)then
        jt = this%charge%lattice%iz(jj)
        vet(:) = (this%charge%lattice%cr(:,jj)-this%charge%lattice%cr(:,ia))*this%charge%lattice%alat
        !write(123,'(3f10.6)') vet(:)
        !write(123,'(3f10.6)') this%charge%lattice%sbarvec(:,m)
        !write(123,'(a,3i4,3f10.6)') 'nn ', IA,m,JJ, VET(:)
        call this%hmfind(vet,nr,hhh,m,ia,m,ni,ntype)
        if(ni==0)then
          this%charge%lattice%nn(ia,m) = 0
        end if
        call this%ham0m_nc(it,jt,vet,hhh)
        do mdir =1,4
          call hcpx(this%hhmag(:,:,mdir),'cart2sph')
          this%hmag(:,:,m,mdir) = this%hhmag(:,:,mdir)
        end do
      end if      
    end do
!    do m=1,nr
!      write(123,*)'m=',m  
!      do mdir=1,4
!        write(123,*)'mdir=',mdir
!        do i=1,9
!          write(123,'(9f10.4)')(real(this%hmag(i,j,m,mdir)),j=1,9)
!        end do
!      end do
!    end do
  end subroutine chbar_nc

  subroutine hmfind(this,vet,nr,hhh,m,ia,jn,ni,ntype)
    class(hamiltonian), intent(inout) :: this
    ! Input
    integer, intent(in) :: ntype ! Atom type
    integer, intent(in) :: m ! Number of the given neighbour
    integer, intent(in) :: ia ! Atom number in clust
    integer, intent(in) :: jn ! ?
    integer, intent(in) :: nr ! Number of neighbours
    real(rp), dimension(3), intent(in) :: vet
    ! Output
    integer, intent(out) :: ni 
    real(rp), dimension(9,9), intent(inout) :: hhh
    ! Local variables
    real(rp) :: a1, a2, a3, aaa, eps
    integer :: i, ilm, jlm

    eps = 0.0001d0
    ni = 1
    a1=0.0d0
    a2=0.0d0
    a3=0.0d0
    aaa=0.0d0
    do i=1,nr
      !write(123,'(a,i4,3f10.4)')'i',i, this%charge%lattice%sbarvec(:,i)
      a1 = (vet(1)-this%charge%lattice%sbarvec(1,i))
      a2 = (vet(2)-this%charge%lattice%sbarvec(2,i))
      a3 = (vet(3)-this%charge%lattice%sbarvec(3,i))
      aaa = a1**2 + a2**2 + a3**2
      if(aaa<eps) goto 1000
    end do
    write(*,'(1x,a,i4,a,i4,a,3f10.6)') ' Error in hamiltonian%hmfind: Neighbour vector not found for atom', ia, &
                                       ' neighbour', jn, 'vector', vet(:)

    ni=0
    1000 continue 
    do ilm = 1,9
      do jlm = 1,9
        hhh(ilm,jlm) = real(this%charge%lattice%sbar(jlm,ilm,m,this%charge%lattice%num(ia)))
      end do
    end do 

    !do ilm = 1,9
    !    write(123,'(9f10.6)')(hhh(ilm,jlm),jlm=1,9)
    !end do
  end subroutine hmfind
end module hamiltonian_mod
