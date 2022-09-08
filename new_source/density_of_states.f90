!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Green
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
!> Module to handle the calculation of the density of states
!------------------------------------------------------------------------------

module density_of_states_mod

    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use control_mod
    use lattice_mod
    use self_mod
    use symbolic_atom_mod
    use recursion_mod
    use precision_mod, only:rp
    use math_mod
    implicit none
  
    private

    !> Module's main structure
    type, public :: dos
      ! General variables
      real(rp), dimension(:,:,:), allocatable :: doscheb

      !> Recursion
      type(recursion), pointer :: recursion
      !> Self
      type(self), pointer :: self
      !> Symbolic atom
      type(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Lattice
      type(lattice), pointer :: lattice
      !> Control
      type(control), pointer :: control
    contains
      procedure :: density
      procedure :: chebyshev_dos
      procedure :: chebyshev_dos_full
      procedure :: restore_to_default
      final :: destructor
    end type

    interface dos
      procedure :: constructor
    end interface dos

contains

  function constructor(recursion_obj,self_obj) result(obj)
    type(dos) :: obj
    type(recursion), target, intent(in) :: recursion_obj
    type(self), target, intent(in) :: self_obj

    obj%recursion => recursion_obj
    obj%self => self_obj
    obj%symbolic_atom => recursion_obj%hamiltonian%symbolic_atom
    obj%lattice => recursion_obj%lattice
    obj%control => recursion_obj%lattice%control

    call obj%restore_to_default()
  end function constructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine destructor(this)
    type(dos) :: this
    if (allocated(this%doscheb)) deallocate(this%doscheb)
  end subroutine destructor


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Uses the moments form the Chebyshev recursions to calculate the DOS
  !---------------------------------------------------------------------------
  subroutine chebyshev_dos(this)
    class(dos), intent(inout) :: this
    ! Local variables
    real(rp), dimension(this%control%lld*2+2) :: kernel
    real(rp), dimension(this%self%channels_ldos+10,0:this%control%lld*2+2) :: polycheb
    real(rp), dimension(this%self%channels_ldos+10) :: w, wscale
    real(rp) :: wstep, eps, wmin, wmax, a, b
    !complex(rp), dimension(this%lattice%nrec,18,18,this%self%channels_ldos+10) :: green
    integer :: i, j, k, l, m, n
    eps = 0.0001d0

    ! Defining rescaling coeficients
    a = (this%self%energy_max-this%self%energy_min)/(2-0.3)
    b = (this%self%energy_max+this%self%energy_min)/2

    wscale(:) = (this%self%ene(:)-b)/a

    ! Calculating the Jackson Kernel
    call jackson_kernel((this%control%lld)*2+2,kernel)

    ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld,kernel,4.0d0)


    do n=1,this%lattice%nrec ! Loop on self-consistent atoms
      ! Multiply the moments with the kernel
      do i=1,18
        this%recursion%mu_ng(n,:,i,i) = this%recursion%mu_n(n,:,i,i)*kernel(:)
      end do
 
      this%recursion%mu_ng(n,2:size(kernel),:,:) = this%recursion%mu_ng(n,2:size(kernel),:,:)*2

      do i=1,size(kernel)
        write(400+n,*) i, sum(this%recursion%mu_n(n,i,1:18,1:18))
      end do

 
      ! Calculate the Chebyshev polynomials
      call t_polynomial(size(w),size(kernel),wscale(:),polycheb)
 
      ! Calculate the density of states
      do i=1,size(kernel)
        do l=1,18
          this%doscheb(n,l,:) = this%doscheb(n,l,:) + this%recursion%mu_ng(n,i,l,l)*polycheb(:,i-1)
          !do m=1,18
          !  green(n,l,m,:) = green(n,l,m,:) + (-i_unit*(this%recursion%mu_ng(n,i,l,m)*exp(-i_unit*(i-1)*acos(wscale(:)))))
          !end do
        end do
      end do
      do l=1,18
        this%doscheb(n,l,:) = this%doscheb(n,l,:)/((sqrt((a**2)-((this%self%ene(:)-b)**2)))*pi) 
        !do m=1,18
        !green(n,l,m,:) = green(n,l,m,:)/((sqrt((a**2)-((this%self%ene(:)-b)**2))))
        !end do
      end do
      do l=1,18
        do i=1,this%self%channels_ldos+10
          if (isnan(this%doscheb(n,l,i))) this%doscheb(n,l,i)=0.0d0
        end do 
      end do
 
      do i=1,this%self%channels_ldos+10
        write(200+n,'(8f16.4)') this%self%ene(i), (this%doscheb(n,1,i)), sum(this%doscheb(n,2:4,i)), sum(this%doscheb(n,5:9,i)),&
                                                  (this%doscheb(n,10,i)), sum(this%doscheb(n,11:13,i)), sum(this%doscheb(n,14:18,i)),&
                                                  sum(this%doscheb(n,1:18,i))
        !write(125+n,'(10f10.6)') this%self%ene(i), real(green(n,1,1,i)),real(green(n,2,2,i)),real(green(n,3,3,i)),real(green(n,4,4,i)),real(green(n,5,5,i)),&
        !                                           &real(green(n,6,6,i)),real(green(n,7,7,i)),real(green(n,8,8,i)),real(green(n,9,9,i))  
      end do
    end do  ! End loop on self-consistent atoms
  end subroutine chebyshev_dos

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Uses the moments form the Chebyshev recursions to calculate the DOS
  !---------------------------------------------------------------------------
  subroutine chebyshev_dos_full(this)
    class(dos), intent(inout) :: this
    ! Local variables
    real(rp), dimension(this%control%lld) :: kernel
    real(rp), dimension(this%self%channels_ldos+10,0:this%control%lld) :: polycheb
    real(rp), dimension(this%self%channels_ldos+10) :: w, wscale
    real(rp) :: wstep, eps, wmin, wmax, a, b
    real(rp), dimension(this%control%lld*2+2,5) :: mu_dum
    complex(rp), dimension(this%lattice%nrec,18,18,this%self%channels_ldos+10) :: green
    integer :: i, j, k, l, m, n
    eps = 0.0001d0

    ! Defining rescaling coeficients
    a = (this%self%energy_max-this%self%energy_min)/(2-0.3)
    b = (this%self%energy_max+this%self%energy_min)/2

    wscale(:) = (this%self%ene(:)-b)/a

    ! Calculating the Jackson Kernel
    call jackson_kernel(this%control%lld,kernel)

    ! Calculating the Lorentz Kernel
!    call lorentz_kernel(this%control%lld,kernel,4.0d0)


    do n=1,this%lattice%nrec ! Loop on self-consistent atoms
      ! Multiply the moments with the kernel
      do i=1,18
        this%recursion%mu_ng(n,:,i,i) = this%recursion%mu_n(n,:,i,i)*kernel(:)
      end do
 
      this%recursion%mu_ng(n,2:size(kernel),:,:) = this%recursion%mu_ng(n,2:size(kernel),:,:)*2

      do i=1,size(kernel)
        write(400+n,*) i, this%recursion%mu_n(n,i,5,6) 
      end do

 
      ! Calculate the Chebyshev polynomials
      call t_polynomial(size(w),size(kernel),wscale(:),polycheb)
 
      ! Calculate the density of states
      do i=1,size(kernel)
        do l=1,18
          this%doscheb(n,l,:) = this%doscheb(n,l,:) + this%recursion%mu_ng(n,i,l,l)*polycheb(:,i-1)
          do m=1,18
            green(n,l,m,:) = green(n,l,m,:) + (-i_unit*(this%recursion%mu_ng(n,i,l,m)*exp(-i_unit*(i-1)*acos(wscale(:)))))
          end do
        end do
      end do
      do l=1,18
        this%doscheb(n,l,:) = this%doscheb(n,l,:)/((sqrt((a**2)-((this%self%ene(:)-b)**2)))*pi) 
        do m=1,18
          green(n,l,m,:) = green(n,l,m,:)/((sqrt((a**2)-((this%self%ene(:)-b)**2))))
        end do
      end do
      do l=1,18
        do i=1,this%self%channels_ldos+10
          if (isnan(this%doscheb(n,l,i))) this%doscheb(n,l,i)=0.0d0
        end do 
      end do
 
      do i=1,this%self%channels_ldos+10
        write(200+n,'(8f16.4)') this%self%ene(i), (this%doscheb(n,1,i)), sum(this%doscheb(n,2:4,i)), sum(this%doscheb(n,5:9,i)),&
                                                  (this%doscheb(n,10,i)), sum(this%doscheb(n,11:13,i)), sum(this%doscheb(n,14:18,i)),&
                                                  sum(this%doscheb(n,1:18,i))
        write(125+n,'(10f10.6)') this%self%ene(i), real(green(n,1,1,i)),real(green(n,2,2,i)),real(green(n,3,3,i)),real(green(n,4,4,i)),real(green(n,5,5,i)),&
                                                   &real(green(n,5,5,i)),real(green(n,6,6,i)),real(green(n,8,8,i)),aimag(green(n,9,9,i))  
      end do
    end do  ! End loop on self-consistent atoms
  end subroutine chebyshev_dos_full


  ! Member functions
  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Calculates the density of states
  !---------------------------------------------------------------------------
  subroutine density(this,tdens,ia,mdir)
    class(dos) :: this
    ! Input
    integer, intent(in) :: mdir ! Direction index
    integer, intent(in) :: ia ! Atom type
    ! Output
    real(rp), dimension(18,this%self%channels_ldos+10), intent(out) :: tdens
    ! Local variables
    integer :: icode,iii,k,l,ll1,nb,nbp1,nl,nt,eidx, ne, nq, ll_in, ifail, npts, nw
    real(rp) :: a1,a2,dens,emax,emin,eps,err, e_shift, e_canon, prefac
    complex(rp) :: dens_i

    integer, dimension(18) :: nb2
    integer, dimension(18) :: jc
    integer, dimension(10*this%lattice%ntype*this%control%lld) :: iwk
    integer, dimension(18) :: ll_fail
    real(rp), dimension(this%control%lld) :: aa,am,bb,bm, sqbb
    real(rp), dimension(10) :: edge,weight,width
    real(rp), dimension(200,10) :: bwk
    real(rp), dimension(18,this%control%lld) :: am2,bm2
    real(rp), dimension(18,10) :: edge2,width2,weight2
    real(rp), dimension(10*this%lattice%ntype*this%control%lld,2,5) :: work
    real(rp), dimension(2) :: bandedges
    real(rp), dimension(18) :: alpha_inf, beta_inf

    eps = 1.0d-14
    err = 0.00001d0
    nbp1 = 2
    ll_in = this%control%lld
    ll_fail = ll_in
    nw = 10*this%lattice%ntype*this%control%lld
    npts = this%self%channels_ldos+10

    do nl = 1,18
      do l = 1,this%control%lld
        aa(l) = this%recursion%a(l,nl,ia,mdir)
        bb(l) = this%recursion%b2(l,nl,ia,mdir)!**2
        sqbb(l) = sqrt(bb(l))
      end do
      if(nl==1.or.nl==10) bb=1.025d0*bb
      ll1 = this%control%lld
      am=0.0d0;bm=0.0d0;edge=0.0d0;width=0.0d0;weight=0.0d0
      err = 0.0000001d0
      call this%recursion%bpOPT(this%control%lld,AA,sqBB,this%control%lld-1,am(1),bm(1),ifail)
      if(nl==1.or.nl==10) bm=1.01d0*bm
      alpha_inf(nl)=am(1)
      beta_inf(nl)=bm(1)
      nb=1
      edge(1)=am(1)-2.0d0*bm(1)
      width(1)=4.0d0*bm(1)
      weight(1)=1.0d0
      nb2(nl)=nb
      am=AA
      bm=BB
      do k = 1,nb
        a1 = edge(k)
        a2 = edge(k) + width(k)
        edge2(nl,k) = edge(k)
        width2(nl,k) = width(k)
        weight2(nl,k) = weight(k)
      end do
      if (nl == 1) then
        emin = a1
        emax = a2
      else
        emin = MIN(emin,a1)
        emax = MAX(emax,a2)
      end if
      if(nb>0) then
        do l = 1,this%control%lld
          am2(nl,l) = am(l)
          bm2(nl,l) = bm(l)
        end do
      else
        do l = 1,this%control%lld
          am2(nl,l) = 0.0d0
          bm2(nl,l) = 0.0d0
        end do
      end if
    end do


    eps = 1.0d-10
    tdens = 0.0d0

    do eidx=1,npts
      do nl = 1,18
        nb = nb2(nl)
        if(nb>0) then
          do l = 1,this%control%lld
            aa(l) = this%recursion%a(l,nl,ia,mdir)
            bb(l) = this%recursion%b2(l,nl,ia,mdir) !*0.5d0
            am(l) = am2(nl,l)
            bm(l) = bm2(nl,l) !*0.5d0
          end do
          edge=0.0d0;width=0.0d0;weight=0.0d0
          do k = 1,nb
            edge(k) = edge2(nl,k)
            width(k) = width2(nl,k)
            weight(k) = weight2(nl,k)
          end do
          e_shift=this%self%ene(eidx)/this%symbolic_atom(ia)%potential%dw_l(nl)-1.00d0*this%symbolic_atom(ia)%potential%cshi(nl)
          prefac=1.0d0
          bandedges(1)=edge(1)!-0.1
          bandedges(2)=edge(1)+width(1)!+0.1
          ! only terminator 5 implemented
          dens=bprldos(e_shift,AA,BB,this%control%lld,bandedges)
        else 
          dens = 0.0d0
        end if
        tdens(nl,eidx) = tdens(nl,eidx) + prefac*dens/this%symbolic_atom(ia)%potential%dw_l(nl)
      end do
    end do


   ! do eidx=1,npts
   !    write(300,*) this%self%ene(eidx), sum(tdens(1:9,eidx)), sum(tdens(10:18,eidx))
   ! end do
  end subroutine density

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Calculates the density of states
  !---------------------------------------------------------------------------
  real(rp) function bprldos(e,a,b2,ll,ei)
    ! Input
    integer, intent(in) :: LL
    real(rp), intent(in) :: E
    real(rp), intent(in) :: A(LL),B2(LL)
    real(rp), intent(in) :: ei(2)
    real(rp) :: DISC,BI2,AIT,RT,T1,T2
    complex(rp) :: etop, ebot,ea,eb,emid
    complex(rp) :: det, zoff, Qt
    integer :: l
     
    ebot=ei(1) !a(ll)-2*b2(ll)
    etop=ei(2) !a(ll)+2*b2(ll)
    !     ebot=a(ll)-2*b2(ll)
    !     etop=a(ll)+2*b2(ll)
    emid= 0.5d0*(etop+ebot)
    ea=e-etop
    eb=e-ebot
    det=ea*eb
!    zoff = (0.0d0,0.01d0)
    zoff=sqrt(det)
    Qt = (e -emid-zoff) * 0.5d0
    if (AIMAG(Qt) > 0.d0) then
       Qt = (e -emid+zoff) * 0.5d0
    end if
    do l=ll-1,1,-1
       Qt = B2(l) / (e-A(l)-Qt)
    end do
    bpRLDOS = -AIMAG(Qt)/PI
    RETURN
  end function bprldos

  subroutine restore_to_default(this)
    class(dos), intent(inout) :: this
    allocate(this%doscheb(this%lattice%nrec,18,this%self%channels_ldos+10))

    this%doscheb(:,:,:) = 0.0d0 
  end subroutine

end module density_of_states_mod
