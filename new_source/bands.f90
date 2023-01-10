!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Bands
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
!> Module to handle calculation of energy bands and parameters related
!with it
!------------------------------------------------------------------------------

module bands_mod

    use control_mod
    use energy_mod
    use green_mod
    use lattice_mod
    use symbolic_atom_mod
    use density_of_states_mod
    use recursion_mod
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use precision_mod, only:rp
    use logger_mod, only: g_logger
    use math_mod
    use string_mod
    implicit none
   
    private
   
    !> Module's main structure
    type, public :: bands
      !> Green
      class(green), pointer :: green
      !> Lattice
      class(lattice), pointer :: lattice
      !> Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Density of states
      class(dos), pointer :: dos
      !> Energy
      class(energy), pointer :: en
      !> Control
      class(control), pointer :: control
      !> Recursion
      class(recursion), pointer :: recursion

      ! General variables
      !> Energy variable
      real(rp) :: e1, e1cheb
      !> Total number of valence electrons
      real(rp) :: qqv      
      !> Band energy mesh
      integer :: nv1, ik1, nv1cheb, ik1cheb
      !> Total density of states
      real(rp), dimension(:), allocatable :: dtot, dtotcheb
      !> Projected DOS
      real(rp), dimension(:,:), allocatable :: dx, dy, dz, dnmag, dup, ddw
      !> Projected Green Function
      complex(rp), dimension(:,:,:,:), allocatable :: g0_x, g0_y, g0_z     
      !> Energy bands (?)
      real(rp), dimension(:,:,:), allocatable :: dspd 
    contains
      procedure :: calculate_projected_green
      procedure :: calculate_projected_dos
      procedure :: calculate_moments
      procedure :: calculate_pl
      procedure :: calculate_moments_chebgauss
      procedure :: calculate_magnetic_moments
      procedure :: calculate_fermi
      procedure :: fermi
      procedure :: restore_to_default
      final :: destructor
    end type 

    interface bands
      procedure :: constructor
    end interface bands

  contains

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Constructor
    !
    !> @param[in] fname Namelist file
    !> @return type(bands)
    !---------------------------------------------------------------------------
    function constructor(green_obj) result(obj)
      type(bands) :: obj
      type(green), target, intent(in) :: green_obj

      obj%green => green_obj
      obj%lattice => green_obj%dos%recursion%lattice
      obj%symbolic_atom => green_obj%dos%recursion%hamiltonian%charge%lattice%symbolic_atoms
      obj%dos => green_obj%dos
      obj%en => green_obj%dos%en
      obj%control => green_obj%dos%recursion%lattice%control 
      obj%recursion => green_obj%dos%recursion

      call obj%restore_to_default() 
    end function constructor

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Destructor
    !---------------------------------------------------------------------------
    subroutine destructor(this)
      type(bands) :: this
      if(allocated(this%dtot)) deallocate(this%dtot)
      if(allocated(this%dtotcheb)) deallocate(this%dtotcheb)
      if(allocated(this%dx)) deallocate(this%dx)
      if(allocated(this%dy)) deallocate(this%dy)
      if(allocated(this%dz)) deallocate(this%dz)
      if(allocated(this%dnmag)) deallocate(this%dnmag)
      if(allocated(this%dup)) deallocate(this%dup)
      if(allocated(this%ddw)) deallocate(this%ddw)
      if(allocated(this%g0_x)) deallocate(this%g0_x)
      if(allocated(this%g0_y)) deallocate(this%g0_y)
      if(allocated(this%g0_z)) deallocate(this%g0_z)
      if(allocated(this%dspd)) deallocate(this%dspd)
    end subroutine destructor

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
      class(bands) :: this

      allocate(this%dtot(this%en%channels_ldos+10))
      allocate(this%dtotcheb(this%en%channels_ldos+10))
      allocate(this%dx(this%en%channels_ldos+10,this%lattice%ntype))
      allocate(this%dy(this%en%channels_ldos+10,this%lattice%ntype))
      allocate(this%dz(this%en%channels_ldos+10,this%lattice%ntype))      
      allocate(this%g0_x(9,9,this%en%channels_ldos+10,this%lattice%ntype))
      allocate(this%g0_y(9,9,this%en%channels_ldos+10,this%lattice%ntype))
      allocate(this%g0_z(9,9,this%en%channels_ldos+10,this%lattice%ntype))
      allocate(this%dspd(6,this%en%channels_ldos+10,this%lattice%ntype))

      this%dtot(:) = 0.0d0
      this%dtotcheb(:) = 0.0d0
      this%dx(:,:) = 0.0d0
      this%dy(:,:) = 0.0d0
      this%dz(:,:) = 0.0d0
      this%g0_x(:,:,:,:) = 0.0d0
      this%g0_y(:,:,:,:) = 0.0d0
      this%g0_z(:,:,:,:) = 0.0d0
      this%dspd(:,:,:) = 0.0d0
    end subroutine restore_to_default

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Calculates the Fermi energy
    !---------------------------------------------------------------------------
    subroutine calculate_fermi(this)
      class(bands) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, ia, ik1_mag, ik1, nv1, ifail
      real(rp) :: e1_mag, ef_mag, e1
 
      e1_mag = 0.0d0
      ef_mag = 0.0d0
      ik1_mag = 0
      ik1 = 0
      this%dtot(:) = 0.0d0

      ik1 = this%en%ik1
      nv1 = this%en%nv1
     
      ! Define the valence electrons from the bulk parameters
      this%qqv = real(sum(this%symbolic_atom(1:this%lattice%nbulk_bulk)%element%valence))

      ! Calculate the total density of states
      do ia=1,this%lattice%nrec
        do i=1,this%en%channels_ldos+10
          do j=1,9
            this%dtot(i) = this%dtot(i) - aimag(this%green%g0(i,j,j,ia)+this%green%g0(i,j+9,j+9,ia))/pi
            this%dtotcheb(i) = this%dtotcheb(i) + (this%dos%doscheb(ia,j,i)+this%dos%doscheb(ia,j+9,i))
          end do
          write(125,'(2f10.5)') this%en%ene(i), this%dtot(i)
        end do
      end do 

      ! Calculate the Fermi enery
      ef_mag = this%en%fermi 
      this%en%chebfermi = this%en%fermi
      if(.not.(this%en%fix_fermi).and.this%control%calctype=='B')then
        e1_mag = ef_mag
        call this%fermi(ef_mag,this%en%edel,ik1_mag,this%en%energy_min,this%en%channels_ldos+10,this%dtot,ifail,this%qqv,e1_mag)
        e1 = this%en%fermi
        call this%fermi(this%en%fermi,this%en%edel,ik1,this%en%energy_min,this%en%channels_ldos+10,this%dtot,ifail,this%qqv,e1_mag)
        this%nv1 = ik1
        this%e1 = e1_mag
        call g_logger%info('Fermi energy: '//real2str(this%en%fermi),__FILE__,__LINE__)
        !For Chebyshev DOS
        ef_mag = this%en%fermi
        e1_mag = ef_mag
        call this%fermi(ef_mag,this%en%edel,ik1_mag,this%en%energy_min,this%en%channels_ldos+10,this%dtotcheb,ifail,this%qqv,e1_mag)
        e1 = this%en%chebfermi
        call this%fermi(this%en%chebfermi,this%en%edel,ik1,this%en%energy_min,this%en%channels_ldos+10,this%dtotcheb,ifail,this%qqv,e1_mag)
        this%nv1cheb = ik1
        this%e1cheb = e1_mag
      else if(this%control%calctype=='B'.and.(this%en%fix_fermi))then
        ik1 = nint((this%en%fermi-this%en%energy_min)/this%en%edel)
        e1 = this%en%energy_min+(ik1-1)*this%en%edel
        this%nv1 = ik1
        this%e1 = e1      
        this%nv1cheb = ik1
        this%e1cheb = e1
      end if
    end subroutine calculate_fermi

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Auxiliar routine to calculate the Fermi enery
    !---------------------------------------------------------------------------
    subroutine fermi(this,EF,H,IK1,AINF,NPTS,Y,IFAIL,QQV,E1)
      class(bands) :: this
      ! Input
      integer, intent(in) :: NPTS
      integer, intent(inout) :: IK1
      integer, intent(out) :: IFAIL
      real(rp), intent(in) :: AINF,H,QQV
      real(rp), intent(inout) :: EF
      real(rp), intent(out) :: E1
      real(rp), dimension(NPTS), intent(in) :: Y
      ! Local variables    
      integer :: i
      real(rp) :: AINT,AINT0,ALPHA

      IFAIL = 1
      AINT = 0.d0
      AINT0 = 0.d0
      do I = 2,NPTS-1,2
         AINT = AINT + H*(Y(I-1)+4.d0*Y(I)+Y(I+1))/3.d0 !*fermifun(AINF+(I)*H,Ef,kBT)
         if (AINT >= QQV) goto 1000
         AINT0 = AINT
      end do
      return
      1000 IFAIL = 0
      !print '(3x,a,2f12.6)','Fermil', QQV,AINT
      if (AINT == QQV) then
         IK1 = I + 1
         EF = AINF + H*I
         E1 = EF
      else
         ALPHA = (AINT-AINT0) / 2.d0 / H
         IK1 = I - 1
         E1 = AINF + H*(I-2)
         EF = ((QQV-AINT0)/ALPHA) + E1
      end if
      !print '(3x,a,4f12.6)','Fermil', QQV,AINT,EF,E1
    end subroutine fermi

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Calculates the moments m^(q), q = 0,1 and 2
    !---------------------------------------------------------------------------
    subroutine calculate_moments(this)
      class(bands) :: this
      ! Local variables
      integer :: i, j, l, m, o ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index
      integer :: isp, soff, jo, nsp
      real(rp) :: sgef, pmef, smef, isgn
      real(rp), dimension(this%en%channels_ldos+10) :: y 
      real(rp), dimension(18) :: chebmom(18), chebmom1(18), chebmom2(18)

      call this%calculate_magnetic_moments()

      this%dspd(:,:,:) = 0.0d0
      do na=1,this%lattice%nrec
        do isp=1,2  
           isgn=(-1.0d0)**(isp-1)
           soff=3*(isp-1)
           do l=1,3
              do m=1,2*l-1
                 o=(l-1)**2+m
                 do ie=1,this%en%channels_ldos
                    this%dspd(l+soff,ie,na)=this%dspd(l+soff,ie,na)-aimag(this%green%g0(ie,o,o,na)+this%green%g0(ie,o+9,o+9,na))-&
                       isgn*this%symbolic_atom(na)%potential%mom(3)*aimag(this%green%g0(ie,o,o,na)-this%green%g0(ie,o+9,o+9,na)) &
                      -isgn*this%symbolic_atom(na)%potential%mom(2)*aimag(i_unit*this%green%g0(ie,o,o+9,na)-i_unit*this%green%g0(ie,o+9,o,na)) &
                      -isgn*this%symbolic_atom(na)%potential%mom(1)*aimag(this%green%g0(ie,o,o+9,na)+this%green%g0(ie,o+9,o,na))
                 end do
              end do
           end do
        end do
      end do

      this%dspd(:,:,:) = this%dspd(:,:,:)*0.5d0/pi

!    do i=1,this%self%channels_ldos
!      write(301,*)this%self%ene(i), sum(this%dspd(:,i,1))
!    end do
!      do na=1,this%lattice%ntype
!         do i=1,this%self%channels_ldos
!            write(7000+na,'(f10.6,2x,18f12.6)') this%self%ene(i),(this%Dspd(jo,i,na),jo=1,3),(-this%Dspd(jo,i,na),jo=4,6)
!         end do
!         rewind(7000+na)
!      end do


      do na=1,this%lattice%nrec
        do i=1,6
          y(:) = 0.0d0
          if(i>3)then
            nsp = 2
          else
            nsp = 1
          end if  
          soff=3*(nsp-1)
          y(:) = this%dspd(i,:,na) 
          sgef = 0.0d0 ; pmef = 0.0d0 ; smef = 0.0d0
          call simpson_m(sgef,this%en%edel,this%en%fermi,this%nv1,y,this%e1,0,this%en%ene)
          call simpson_m(pmef,this%en%edel,this%en%fermi,this%nv1,y,this%e1,1,this%en%ene)
          call simpson_m(smef,this%en%edel,this%en%fermi,this%nv1,y,this%e1,2,this%en%ene)

          this%symbolic_atom(na)%potential%gravity_center(i-soff,nsp) = (pmef/sgef) - this%symbolic_atom(na)%potential%vmad
          this%symbolic_atom(na)%potential%ql(1,i-soff-1,nsp) = sgef
          this%symbolic_atom(na)%potential%ql(2,i-soff-1,nsp) = 0.0d0
          this%symbolic_atom(na)%potential%ql(3,i-soff-1,nsp) = smef - 2.0d0*(pmef/sgef)*pmef + ((pmef/sgef)**2)*sgef
        end do
      end do 
!      chebmom = 0.0d0
!      chebmom1 = 0.0d0
!      chebmom2 = 0.0d0
!    do na=1,this%lattice%nrec
!      do i=1,18
!        call simpson_m(sgef,this%self%edel,this%self%chebfermi,this%nv1cheb,this%dos%doscheb(na,i,:),this%e1cheb,0,this%self%ene)
!        call simpson_m(pmef,this%self%edel,this%self%chebfermi,this%nv1cheb,this%dos%doscheb(na,i,:),this%e1cheb,1,this%self%ene)
!        call simpson_m(smef,this%self%edel,this%self%chebfermi,this%nv1cheb,this%dos%doscheb(na,i,:),this%e1cheb,2,this%self%ene)
!        chebmom(i) = sgef
!        chebmom1(i) = pmef
!        chebmom2(i) = smef !- 2.0d0*(pmef/sgef)*pmef + ((pmef/sgef)**2)*sgef)
!      end do

      !write(*,*) this%symbolic_atom(1)%potential%gravity_center(1,1)
!      write(102,*) sum(this%symbolic_atom(1)%element%q0(:,1))+sum(this%symbolic_atom(1)%element%q0(:,2)), &
!                   sum(this%symbolic_atom(1)%element%q1(:,1))+sum(this%symbolic_atom(1)%element%q1(:,2)), &
!                   sum(this%symbolic_atom(1)%element%q2(:,1))+sum(this%symbolic_atom(1)%element%q2(:,2))
!      write(101,*) chebmom(1), sum(chebmom(2:4)), sum(chebmom(5:9)), chebmom(10), sum(chebmom(11:13)), sum(chebmom(14:18))
!      write(101,*) this%qqv
!    end do
    end subroutine calculate_moments


    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Calculates the band moments using the Chebyshev-Gauss quadrature method
    !---------------------------------------------------------------------------
    subroutine calculate_moments_chebgauss(this)
      class(bands), intent(inout) :: this
      ! Local variables
      real(rp), dimension(this%control%lld) :: x_i, w_i, q_i 
      real(rp), dimension(this%control%lld) :: w, wscale
      real(rp), dimension(this%control%lld,0:this%control%lld) :: polycheb
      real(rp), dimension(this%lattice%ntype,18,this%control%lld) :: doscheb_i
      real(rp), dimension(this%lattice%ntype,18) :: q0l, q1l, q2l
      real(rp) :: a, b, mom0, mom1, mom2
      integer :: i, l, llplusone, n

      x_i(:) = 0.0d0
      w_i(:) = 0.0d0
      q_i(:) = 0.0d0

      call chebyshev_gauss_quadrature(this%control%lld,x_i(:),q_i(:)) 

      w_i(:) = ((0.5*(x_i(:)+1))*(this%en%energy_min))
      w_i(:) = w_i(:) + this%en%chebfermi
      q_i(:) = -((q_i(:)*(this%en%energy_min))/2)
      ! Defining rescaling coeficients
      a = (this%en%energy_max-this%en%energy_min)/(2-0.3)
      b = (this%en%energy_max+this%en%energy_min)/2
        
      wscale(:) = ((w_i(:)-b)/a)
      ! Calculate the Chebyshev polynomials
      call t_polynomial(this%control%lld,this%control%lld,wscale(:),polycheb)
 
      doscheb_i(:,:,:) = 0.0d0

      do n=1,this%lattice%nrec
        ! Calculate the density of states
        do l=1,18
          do i=1,this%control%lld
            doscheb_i(n,l,:) = doscheb_i(n,l,:) + this%recursion%mu_ng(n,i,l,l)*polycheb(:,i-1)
          end do
        end do
        do l=1,18
          doscheb_i(n,l,:) = doscheb_i(n,l,:)/((sqrt((a**2)-((w_i(:)-b)**2)))*pi)
        end do
 
        do i=1,this%control%lld
          write(128,*) w_i(i), sum(doscheb_i(n,1:9,i))
        end do
 
        do l=1,18
          do i=1,this%control%lld
            if (isnan(doscheb_i(n,l,i))) doscheb_i(n,l,i)=0.0d0
          end do
        end do
        
        do l=1,18
          q0l(n,l) = sum(doscheb_i(n,l,:)*q_i(:))
          q1l(n,l) = sum(doscheb_i(n,l,:)*w_i(:)*q_i(:))
          q2l(n,l) = sum(doscheb_i(n,l,:)*(w_i(:)**2)*q_i(:)) 
        end do
        write(100,*) sum(q0l(n,:)), sum(q1l(n,:)), sum(q2l(n,:))
      end do
    end subroutine calculate_moments_chebgauss

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Calculates the magnetic moments mx, my, mz and mtot
    !---------------------------------------------------------------------------
    subroutine calculate_magnetic_moments(this)
      class(bands) :: this
      ! Local variables
      integer :: i, j ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index

      call this%calculate_projected_dos() 
     
      do na=1,this%lattice%nrec
        call simpson_m(this%symbolic_atom(na)%potential%mx,this%en%edel,this%en%fermi,this%nv1,this%dx(:,na),this%e1,0,this%en%ene)
        call simpson_m(this%symbolic_atom(na)%potential%my,this%en%edel,this%en%fermi,this%nv1,this%dy(:,na),this%e1,0,this%en%ene)
        call simpson_m(this%symbolic_atom(na)%potential%mz,this%en%edel,this%en%fermi,this%nv1,this%dz(:,na),this%e1,0,this%en%ene)

        this%symbolic_atom(na)%potential%mtot = sqrt((this%symbolic_atom(na)%potential%mx**2) +&
                                                     (this%symbolic_atom(na)%potential%my**2) +&
                                                     (this%symbolic_atom(na)%potential%mz**2))
 
        this%symbolic_atom(na)%potential%mom(1) = this%symbolic_atom(na)%potential%mx/this%symbolic_atom(na)%potential%mtot
        this%symbolic_atom(na)%potential%mom(2) = this%symbolic_atom(na)%potential%my/this%symbolic_atom(na)%potential%mtot
        this%symbolic_atom(na)%potential%mom(3) = this%symbolic_atom(na)%potential%mz/this%symbolic_atom(na)%potential%mtot
      end do

    end subroutine calculate_magnetic_moments

    subroutine calculate_projected_dos(this)
      class(bands) :: this
      ! Local variables
      integer :: i, j ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index

      do na = 1,this%lattice%nrec
        do ie = 1,this%en%channels_ldos+10
          do i=1,9
            this%dz(ie,na) = this%dz(ie,na) - aimag(this%green%g0(ie,i,i,na)-this%green%g0(ie,i+9,i+9,na))/pi
            this%dy(ie,na) = this%dy(ie,na) - aimag(i_unit*this%green%g0(ie,i,i+9,na)-i_unit*this%green%g0(ie,i+9,i,na))/pi
            this%dx(ie,na) = this%dx(ie,na) - aimag(this%green%g0(ie,i,i+9,na)+this%green%g0(ie,i+9,i,na))/pi
         end do
        end do
      end do
    end subroutine calculate_projected_dos

    subroutine calculate_projected_green(this)
      class(bands) :: this
      ! Local variables
      integer :: i, j ! Orbital index
      integer :: na ! Atom index
      integer :: ie ! Energy channel index

      do na = 1,this%lattice%nrec
        do ie = 1,this%en%channels_ldos+10
          do i=1,9
            do j=1,9
              this%g0_z(i,j,ie,na) = this%g0_z(i,j,ie,na) + (this%green%g0(ie,i,i,na)-this%green%g0(ie,i+9,i+9,na))
              this%g0_y(i,j,ie,na) = this%g0_y(i,j,ie,na) + (i_unit*this%green%g0(ie,i,i+9,na)-i_unit*this%green%g0(ie,i+9,i,na))
              this%g0_x(i,j,ie,na) = this%g0_z(i,j,ie,na) + (this%green%g0(ie,i,i+9,na)+this%green%g0(ie,i+9,i,na))
            end do
          end do
        end do
      end do
    end subroutine calculate_projected_green

    subroutine calculate_pl(this)
      class(bands) :: this
      ! Local variables
      integer :: i ! orbital index    
      integer :: is ! spin channel index
      integer :: ia ! atom index
      real(rp) :: rq, dnu, pli, delta2 ! Local variables

      do ia=1,this%lattice%nrec
        do is=1,2
          do i=1,3
            rq = 1/this%symbolic_atom(ia)%potential%qpar(i-1,is)
            delta2 = this%symbolic_atom(ia)%potential%srdel(i-1,is)*this%symbolic_atom(ia)%potential%srdel(i-1,is)
            dnu = (i-1.) + &
                  (2.*(i-1)+1.)/ & 
                  (rq*(this%symbolic_atom(ia)%potential%c(i-1,is)-this%symbolic_atom(ia)%potential%gravity_center(i,is))/2./(2*(i-1)+1.)/ &
                  (this%symbolic_atom(ia)%potential%c(i-1,is)-this%symbolic_atom(ia)%potential%gravity_center(i,is)-delta2*rq)-1.)
            pli = -atan(dnu)/pi + 0.5d0 + INT(this%symbolic_atom(ia)%potential%pl(i-1,is))
            this%symbolic_atom(ia)%potential%pl(i-1,is) = pli
          end do
        end do
      end do  
    end subroutine calculate_pl
end module bands_mod
