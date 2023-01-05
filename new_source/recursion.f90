!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Recursion
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
!> Module to handle recursion calculations
!------------------------------------------------------------------------------


module recursion_mod

    use hamiltonian_mod
    use lattice_mod
    use energy_mod
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use precision_mod, only: rp
    use math_mod
    implicit none
   
    private
  
    !> Module's main structure
    type, public :: recursion
      !> Hamiltonian
      class(hamiltonian), pointer :: hamiltonian
      !> Lattice
      class(lattice), pointer :: lattice
      !> Energy
      class(energy), pointer :: en
      ! General variables

      !> Scalar recursion coefficients
      real(rp), dimension(:,:,:,:), allocatable :: a, b2
      real(rp), dimension(:), allocatable :: atemp, b2temp
      !> Atom list in hopping region
      integer, dimension(:), allocatable :: izero
      !> Atom list in hopping region as a function of recursion step
      integer, dimension(:,:), allocatable :: izeroll
      !> Wave functions for recursion hopping
      complex(rp), dimension(:,:), allocatable :: psi, pmn  
      !> Wave functions for recursion hopping (Chebyshev)
      complex(rp), dimension(:,:,:), allocatable :: psi0, psi1, psi2
      !> Chebyshev moments
      real(rp), dimension(:,:,:,:), allocatable :: mu_n, mu_ng
      !> Variable to save H|psi>
      complex(rp), dimension(:,:), allocatable :: v
    contains
      procedure :: hop
      procedure :: crecal
      procedure :: recur
      procedure :: bpopt
      procedure :: emami
      procedure :: create_ll_map
      procedure :: chebyshev_recur
      procedure :: chebyshev_recur_full
      procedure :: restore_to_default
      final :: destructor
    end type recursion
  
    interface recursion
      procedure :: constructor
    end interface recursion
  
  contains
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Constructor
    !
    !> @param[in] fname Namelist file
    !> @return type(recursion)
    !---------------------------------------------------------------------------
    function constructor(hamiltonian_obj,energy_obj) result(obj)
      type(recursion) :: obj
      type(hamiltonian), target, intent(in) :: hamiltonian_obj
      type(energy), target, intent(in) :: energy_obj 
      
      obj%hamiltonian => hamiltonian_obj
      obj%lattice => hamiltonian_obj%charge%lattice
      obj%en => energy_obj
      
      call obj%restore_to_default()
    end function constructor
    
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Destructor
    !---------------------------------------------------------------------------
    subroutine destructor(this)
      type(recursion) :: this
      if(allocated(this%a)) deallocate(this%a)
      if(allocated(this%b2)) deallocate(this%b2)
      if(allocated(this%izero)) deallocate(this%izero)
      if(allocated(this%izeroll)) deallocate(this%izeroll)
      if(allocated(this%psi)) deallocate(this%psi)
      if(allocated(this%psi1)) deallocate(this%psi1)
      if(allocated(this%psi2)) deallocate(this%psi2)
      if(allocated(this%pmn)) deallocate(this%pmn)
      if(allocated(this%v)) deallocate(this%v)
      if(allocated(this%mu_n)) deallocate(this%mu_n)
      if(allocated(this%mu_ng)) deallocate(this%mu_ng)
      if(allocated(this%atemp)) deallocate(this%atemp)
      if(allocated(this%b2temp)) deallocate(this%b2temp)
    end subroutine destructor
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Recursion method using Chebyshev moments
    !---------------------------------------------------------------------------
    subroutine chebyshev_recur(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: nb, ih, i, j, k, nr, ll, m, n, l, hblocksize, nat, nnmap, nlimplus1
      integer, dimension(0:this%lattice%kk) :: idumll
      complex(rp) :: cone = (1.0d0,0.0d0) 
      complex(rp), dimension(18,18) :: dum, dum1, dum2
      complex(rp), dimension(18,this%lattice%kk) :: v
      complex(rp), dimension(:,:,:,:), allocatable :: hcheb
      real(rp) :: a, b, start, finish
      ! External functions
      complex(rp), external :: zdotc

      hblocksize = 18
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1
      allocate(hcheb(18,18,(this%lattice%nn(1,1)+1),this%lattice%kk))
 
      a = (this%en%energy_max-this%en%energy_min)/(2-0.3)
      b = (this%en%energy_max+this%en%energy_min)/2

      hcheb(:,:,:,:) = this%hamiltonian%ee(:,:,:,:) !hcheb(:,:,:,:)!/cmplx(a,0.d0)

      do i=1,this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
        j = this%lattice%irec(i) ! Atom number in the clust file
        ! Initialize neighbouring map
        this%izeroll(:,:) = 0
        this%izeroll(j,1) = 1

        call cpu_time(start)
        call this%create_ll_map()
        call cpu_time(finish)
        print '("Mapping neighrbours time = ",f10.3," seconds.")',(finish-start)!/32

        ! Initializing wave functions
        this%psi0(:,:,:) = (0.0d0,0.0d0)
        this%psi1(:,:,:) = (0.0d0,0.0d0)
        this%psi2(:,:,:) = (0.0d0,0.0d0)

        do l=1,18
          !> Starting state for |phi_0>
          this%psi0(l,l,j) = (1.0d0,0.0d0)
        end do

        dum(:,:) = (0.0d0,0.0d0)
        ! Write the 0th moment
        call zgemm('c','n',18,18,18,cone,this%psi0(:,:,j),18,this%psi0(:,:,j),18,cone,dum,18)
        this%mu_n(i,1,:,:) = real(dum(:,:))

        ! Write |phi_1>=H|phi_0>
        do k = nlimplus1, this%lattice%kk ! Loop in the clust
          idumll(k) = this%izeroll(k,1)
          ih = this%lattice%iz(k)
          nr = this%lattice%nn(k,1)
          if (this%izeroll(k,1)/=0) then
            call zgemm('n','n',18,18,18,cone,hcheb(1,1,1,ih),18,this%psi0(:,:,k),18,cone,this%psi1(:,:,k),18)
            call zgemm('n','n',18,18,18,cone,this%hamiltonian%lsham(:,:,ih),18,this%psi0(:,:,k),18,cone,this%psi1(:,:,k),18)
          end if
          if(nr>=2)then
            do nb = 2,nr ! Loop in the neighbouring
              nnmap = this%lattice%nn(k,nb)
              if(nnmap/=0.and.this%izeroll(nnmap,1)/=0)then
                call zgemm('n','n',18,18,18,cone,hcheb(1,1,nb,ih),18,this%psi0(:,:,nnmap),18,cone,this%psi1(:,:,k),18)
              end if
            end do ! End of the loop in the neighbouring
          end if
          ! Do the scaling and shifting
          this%psi1(:,:,k) = this%psi1(:,:,k) - b*this%psi0(:,:,k)
          this%psi1(:,:,k) = this%psi1(:,:,k)/a
        end do ! End loop in the clust

        ! Write the 1st moment
        dum(:,:) = (0.0d0,0.0d0)
        do n=1,this%lattice%kk 
          if(this%izeroll(n,1)/=0)then
            call zgemm('c','n',18,18,18,cone,this%psi1(:,:,n),18,this%psi0(:,:,n),18,cone,dum,18)
          end if
        end do
        this%mu_n(i,2,:,:) = real(dum(:,:))

        ! Start the recursion
        do ll=1,this%lattice%control%lld
          write(180,*) 'll=', ll
          ! Write H*|phi_1>
          do k = nlimplus1, this%lattice%kk ! Loop in the clust
            idumll(k) = this%izeroll(k,ll+1)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k,1)
            if (this%izeroll(k,ll+1)/=0) then
              call zgemm('n','n',18,18,18,cone,hcheb(1,1,1,ih),18,this%psi1(:,:,k),18,cone,this%psi2(:,:,k),18)
              call zgemm('n','n',18,18,18,cone,this%hamiltonian%lsham(:,:,ih),18,this%psi1(:,:,k),18,cone,this%psi2(:,:,k),18)
            end if
            if(nr>=2)then
              do nb = 2,nr ! Loop in the neighbouring
                nnmap = this%lattice%nn(k,nb)
                if(nnmap/=0.and.this%izeroll(nnmap,ll+1)/=0)then
                  call zgemm('n','n',18,18,18,cone,hcheb(1,1,nb,ih),18,this%psi1(:,:,nnmap),18,cone,this%psi2(:,:,k),18)
                end if
              end do ! End of the loop in the neighbouring
            end if
            ! Do the scaling and shifting
            this%psi2(:,:,k) = this%psi2(:,:,k) - b*this%psi1(:,:,k)
            this%psi2(:,:,k) = this%psi2(:,:,k)/a
            ! Write 2*H*|phi_1>
            this%psi2(:,:,k) = 2*this%psi2(:,:,k)
          end do ! End loop in the clust

          ! Write |phi_2>=2*H*|phi_1> - |phi_0>
          this%psi2(:,:,:) = this%psi2(:,:,:) - this%psi0(:,:,:)

          ! Write Moments
          dum1(:,:) = (0.0d0,0.0d0)
          dum2(:,:) = (0.0d0,0.0d0)
          do n=1,this%lattice%kk
            if(this%izeroll(n,ll+1)/=0)then
               call zgemm('c','n',18,18,18,cone,this%psi1(:,:,n),18,this%psi1(:,:,n),18,cone,dum1,18) 
               call zgemm('c','n',18,18,18,cone,this%psi2(:,:,n),18,this%psi1(:,:,n),18,cone,dum2,18)
            end if
          end do
          this%mu_n(i,2*ll+1,:,:) = 2*real(dum1(:,:)) - this%mu_n(i,1,:,:)
          this%mu_n(i,2*ll+2,:,:) = 2*real(dum2(:,:)) - this%mu_n(i,2,:,:) 

          this%psi0(:,:,:) = this%psi1(:,:,:)
          this%psi1(:,:,:) = this%psi2(:,:,:)
          this%psi2(:,:,:) = (0.0d0,0.0d0)
        end do ! End loop in the recursion steps 
      end do ! End loop on the number of atoms to be treat self-consistently

      deallocate(hcheb)
    end subroutine chebyshev_recur


    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Recursion method using to find the Chebyshev moments for diagonal and
    !off-diagonal terms. Here, one cannot use the 'double-trick'. See 
    ! Rev. Mod. Phys. 78, 275.
    !---------------------------------------------------------------------------
    subroutine chebyshev_recur_full(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: nb, ih, i, j, k, nr, ll, m, n, l, hblocksize, nat, nnmap, nlimplus1
      integer, dimension(0:this%lattice%kk) :: idumll
      complex(rp) :: cone = (1.0d0,0.0d0) 
      complex(rp), dimension(18,18) :: dum, dum1, dum2
      complex(rp), dimension(18,this%lattice%kk) :: v
      complex(rp), dimension(:,:,:,:), allocatable :: hcheb
      complex(rp), dimension(:,:,:), allocatable :: psiref
      real(rp) :: a, b, start, finish
      ! External functions
      complex(rp), external :: zdotc

      hblocksize = 18
      nat = this%lattice%kk
      nlimplus1 = this%lattice%nmax + 1
      allocate(hcheb(18,18,(this%lattice%nn(1,1)+1),this%lattice%kk),psiref(18,18,this%lattice%kk))
 
      a = (this%en%energy_max-this%en%energy_min)/(2-0.3)
      b = (this%en%energy_max+this%en%energy_min)/2

      hcheb(:,:,:,:) = this%hamiltonian%ee(:,:,:,:) !hcheb(:,:,:,:)!/cmplx(a,0.d0)

      do i=1,this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
        j = this%lattice%irec(i) ! Atom number in the clust file
        ! Initialize neighbouring map
        this%izeroll(:,:) = 0
        this%izeroll(j,1) = 1

        call cpu_time(start)
        call this%create_ll_map()
        call cpu_time(finish)
        print '("Mapping neighrbours time = ",f10.3," seconds.")',(finish-start)!/32

        ! Initializing wave functions
        this%psi0(:,:,:) = (0.0d0,0.0d0)
        this%psi1(:,:,:) = (0.0d0,0.0d0)
        this%psi2(:,:,:) = (0.0d0,0.0d0)

        do l=1,18
          !> Starting state for |phi_0>
          this%psi0(l,l,j) = (1.0d0,0.0d0)
          psiref(l,l,j) = (1.0d0,0.0d0)
        end do

        dum(:,:) = (0.0d0,0.0d0)
        ! Write the 0th moment
        call zgemm('c','n',18,18,18,cone,this%psi0(:,:,j),18,this%psi0(:,:,j),18,cone,dum,18)
        this%mu_n(i,1,:,:) = real(dum(:,:))

        ! Write |phi_1>=H|phi_0>
        do k = nlimplus1, this%lattice%kk ! Loop in the clust
          idumll(k) = this%izeroll(k,1)
          ih = this%lattice%iz(k)
          nr = this%lattice%nn(k,1)
          if (this%izeroll(k,1)/=0) then
            call zgemm('n','n',18,18,18,cone,hcheb(:,:,1,ih),18,this%psi0(:,:,k),18,cone,this%psi1(:,:,k),18)
          end if
          if(nr>=2)then
            do nb = 2,nr ! Loop in the neighbouring
              nnmap = this%lattice%nn(k,nb)
              if(nnmap/=0.and.this%izeroll(nnmap,1)/=0)then
                call zgemm('n','n',18,18,18,cone,hcheb(:,:,nb,ih),18,this%psi0(:,:,nnmap),18,cone,this%psi1(:,:,k),18)
              end if
            end do ! End of the loop in the neighbouring
          end if
          ! Do the scaling and shifting
          this%psi1(:,:,k) = this%psi1(:,:,k) - b*this%psi0(:,:,k)
          this%psi1(:,:,k) = this%psi1(:,:,k)/a
        end do ! End loop in the clust

        ! Write the 1st moment
        dum(:,:) = (0.0d0,0.0d0)
        do n=1,this%lattice%kk 
          if(this%izeroll(n,1)/=0)then
            call zgemm('c','n',18,18,18,cone,psiref(:,:,n),18,this%psi1(:,:,n),18,cone,dum,18)
          end if
        end do
        this%mu_n(i,2,:,:) = real(dum(:,:))

        ! Start the recursion
        do ll=1,this%lattice%control%lld
          ! Write H*|phi_1>
          do k = nlimplus1, this%lattice%kk ! Loop in the clust
            idumll(k) = this%izeroll(k,ll+1)
            ih = this%lattice%iz(k)
            nr = this%lattice%nn(k,1)
            if (this%izeroll(k,ll+1)/=0) then
              call zgemm('n','n',18,18,18,cone,hcheb(:,:,1,ih),18,this%psi1(:,:,k),18,cone,this%psi2(:,:,k),18)
            end if
            if(nr>=2)then
              do nb = 2,nr ! Loop in the neighbouring
                nnmap = this%lattice%nn(k,nb)
                if(nnmap/=0.and.this%izeroll(nnmap,ll+1)/=0)then
                  call zgemm('n','n',18,18,18,cone,hcheb(:,:,nb,ih),18,this%psi1(:,:,nnmap),18,cone,this%psi2(:,:,k),18)
                end if
              end do ! End of the loop in the neighbouring
            end if
            ! Do the scaling and shifting
            this%psi2(:,:,k) = this%psi2(:,:,k) - b*this%psi1(:,:,k)
            this%psi2(:,:,k) = this%psi2(:,:,k)/a
            ! Write 2*H*|phi_1>
            this%psi2(:,:,k) = 2*this%psi2(:,:,k)
          end do ! End loop in the clust

          ! Write |phi_2>=2*H*|phi_1> - |phi_0>
          this%psi2(:,:,:) = this%psi2(:,:,:) - this%psi0(:,:,:)

          ! Write Moments
          dum1(:,:) = (0.0d0,0.0d0)
          do n=1,this%lattice%kk
            if(this%izeroll(n,ll+1)/=0)then
               call zgemm('c','n',18,18,18,cone,psiref(:,:,n),18,this%psi2(:,:,n),18,cone,dum1,18) 
            end if
          end do
          this%mu_n(i,ll+2,:,:) = real(dum1(:,:))

          this%psi0(:,:,:) = this%psi1(:,:,:)
          this%psi1(:,:,:) = this%psi2(:,:,:)
          this%psi2(:,:,:) = (0.0d0,0.0d0)
        end do ! End loop in the recursion steps 
      end do ! End loop on the number of atoms to be treat self-consistently

      deallocate(hcheb)
    end subroutine chebyshev_recur_full

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Calculates the map of neighbouring map as a function of the recursion
    !steps
    !---------------------------------------------------------------------------
    subroutine create_ll_map(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, nr, nnmap, ll
      integer, dimension(0:this%lattice%kk) :: idumll

      idumll(:) = 0

      do ll=1,this%lattice%control%lld
        idumll(:) = 0
        do i=1,this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
          idumll(i) = this%izeroll(i,ll)
          nr = this%lattice%nn(i,1) ! Number of neighbours
          if(nr>=2)then
            do j=2,nr ! Loop on the neighbouring
              nnmap = this%lattice%nn(i,j)
              if(nnmap/=0)then
                if(this%izeroll(nnmap,ll)/=0)then
                  idumll(i) = 1
                end if
              end if
            end do ! End of loop in the neighbouring
          end if
        end do
        this%izeroll(:,ll+1) = idumll(:)
      end do  
    end subroutine create_ll_map

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Calculates the recursion coefficients A.
    !---------------------------------------------------------------------------
    subroutine hop(this,ll)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m, n, nr, nnmap, nlimplus1
      integer :: ll ! Recursion step
      integer :: ino ! Atom type
      integer :: ih ! Atom number in the clust
      integer, dimension(0:this%lattice%kk) :: idum
      complex(rp), dimension(18) :: dum
      real(rp) :: summ, start, finish
 
      summ = 0.0d0
      idum(:) = 0
       
      nlimplus1 = this%lattice%nmax + 1
 
      select case(this%lattice%control%nsp)

      case (1) ! Scalar relativistic case

        if(this%lattice%nmax/=0)then ! In case of impurities using the local hamiltonian
          do i=1,this%lattice%nmax ! Loop in the neighbouring 
            idum(i) = this%izero(i)
            dum(:) = (0.0d0,0.0d0)
            nr = this%lattice%nn(i,1) ! Number of neighbours of atom i
            if(this%izero(i)/=0)then
              do m=1,9 ! Loop on the orbital m
                do l=1,9 ! Loop on the orbital l
                  dum(l)   = dum(l)   + this%hamiltonian%hall(l,m,1,i)*this%psi(m,i)
                  dum(l+9) = dum(l+9) + this%hamiltonian%hall(l+9,m+9,1,i)*this%psi(m+9,i)
                end do ! End of loop on orbital m
              end do ! End of loop on orbital l    
            end if
            if(nr>=2)then
              do j=2,nr ! Loop on the neighbouring
                nnmap = this%lattice%nn(i,j)
                if(nnmap/=0)then 
                  if(this%izero(nnmap)/=0)then
                    do m=1,9 ! Loop in the orbital m
                      do l=1,9 ! Loop in the orbital l
                        dum(l)   = dum(l)   + this%hamiltonian%hall(l  ,m  ,j,i)*this%psi(m,nnmap)
                        dum(l+9) = dum(l+9) + this%hamiltonian%hall(l+9,m+9,j,i)*this%psi(m+9,nnmap)
                      end do ! End of loop in orbital m
                    end do ! End of loop in orbital l
                    idum(i) = 1
                  end if 
                end if
              end do ! End of loop in the neighbouring
            end if
            do l=1,18
              this%v(l,i) = dum(l)
            end do
          end do ! End of loop in the neighbouring   
        end if ! End of local Hamiltonian loop

        do i=nlimplus1,this%lattice%kk ! Loop to find the bulk atoms using the bulk Hamiltonian
          idum(i) = this%izero(i) 
          ih = this%lattice%iz(i) ! Atom type
          dum(:) = (0.0d0,0.0d0)
          nr = this%lattice%nn(i,1) ! Number of neighbours
          if(this%izero(i)/=0)then
            !write(125,*) 'i is ', i
            do m=1,9 ! Loop on the orbital m
              do l=1,9 ! Loop on the orbital l
                dum(l)   = dum(l)   + this%hamiltonian%ee(l  ,m  ,1,ih)*this%psi(m  ,i)
                dum(l+9) = dum(l+9) + this%hamiltonian%ee(l+9,m+9,1,ih)*this%psi(m+9,i)
              end do ! End of the loop on the orbital l
            end do ! End of loop on the orbital m
          end if
          if(nr>=2)then
            do j=2,nr ! Loop on the neighbouring
              nnmap = this%lattice%nn(i,j)
              if(nnmap/=0)then
                if(this%izero(nnmap)/=0)then
                  !write(125,*) i,j,this%lattice%nn(i,j)
                  do m=1,9 ! Loop on the orbital m
                    do l=1,9 ! Loop on orbital l
                      dum(l)   = dum(l)   + this%hamiltonian%ee(l  ,m  ,j,ih)*this%psi(m  ,nnmap)
                      dum(l+9) = dum(l+9) + this%hamiltonian%ee(l+9,m+9,j,ih)*this%psi(m+9,nnmap)
                    end do ! End of loop on orbital l
                  end do ! End of loop on the orbital m
                  idum(i) = 1
                end if
              end if
            end do ! End of loop in the neighbouring
          end if
          do l=1,18
            this%v(l,i) = dum(l)
          end do
        end do

        ! Redefines idum for all clust atoms
        do i=1,this%lattice%kk
          this%izero(i) = idum(i)
          do l=1,18
            dum(l) = this%v(l,i)
            summ = summ + real(dum(l)*conjg(this%psi(l,i)))
            this%pmn(l,i) = dum(l) + this%pmn(l,i)
          end do
        end do
        this%atemp(ll) = summ
      end select
    end subroutine hop 
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Calculates the recursion coefficients B2.
    !---------------------------------------------------------------------------
    subroutine crecal(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, k, l, m
      integer :: ll, llmax ! Recursion step
      integer :: nm1 ! LL-1
      integer :: nat ! Clust size
      integer :: hblocksize ! Hamiltonian size (18)
      real(rp) :: s, summ, start, finish
      complex(rp) :: dum, ajc
      complex(rp), dimension(18,this%lattice%kk) :: thpsi 
      character(len=1) :: transa
      character, dimension(6) :: matdescra
      ! External functions
      complex(rp), external :: zdotc
     
      summ = this%b2temp(1)
      thpsi(:,:) = (0.0d0,0.0d0)

      hblocksize = 18
      nat = this%lattice%kk
      nm1 = this%lattice%control%lld - 1
      llmax = this%lattice%control%lld
      do ll=1,nm1
        call cpu_time(start)
        call this%hop(ll)
        call cpu_time(finish)
        this%b2temp(ll) = summ
        ajc = -this%atemp(ll)

        call zaxpy(nat*hblocksize,ajc,this%psi,1,this%pmn,1)



        summ = 0.0d0
      
        ! In case zdotc function lapack is not working
        do i=1,nat
          do k=1,18
           summ = summ + real(conjg(this%pmn(k,i))*this%pmn(k,i))
          end do
        end do

!        summ = real(zdotc(nat*hblocksize,this%pmn,1,this%pmn,1))

        s = 1.0d0/sqrt(summ)

        thpsi(:,:) = this%pmn(:,:)*s
        this%pmn(:,:) = this%psi(:,:)
        this%psi(:,:) = thpsi(:,:)

        s = sqrt(summ)
      
        this%pmn(:,:) = -this%pmn(:,:)*s
      end do  

      this%b2temp(llmax) = summ
    end subroutine crecal


    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Calculates the recursion coefficients A and B2.
    !---------------------------------------------------------------------------
    subroutine recur(this)
      class(recursion), intent(inout) :: this
      ! Local variables
      integer :: i, j, l, ll
      integer :: llmax ! Recursion steps

      llmax = this%lattice%control%lld

      do i=1,this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
        j = this%lattice%irec(i) ! Atom number in the clust file
        do l = 1,18 ! Loop on the orbitals
        ! Clear list of atoms in hopping region

        this%izero(:) = 0
        ! Initializing wave functions
        this%psi(:,:) = (0.0d0,0.0d0)
        this%pmn(:,:) = (0.0d0,0.0d0)

        this%psi(l,j) = (1.0d0,0.0d0)
        this%izero(j) = 1
        this%atemp(:) = 0.0d0 ; this%b2temp(:) = 0.0d0

        this%b2temp(1) = 1.0d0
        this%atemp(llmax) = 0.0d0
        !write(125,*) 'orbital ', l
        call this%crecal()

        do ll=1,llmax ! Loop on the recursion steps
          this%a(ll,l,i,1) = this%atemp(ll)
          this%b2(ll,l,i,1) = this%b2temp(ll)
        end do ! End of the loop on the recursion steps
      end do ! End of the loop on the orbitals
    end do ! End of the loop on the nrec

     ! For debug purposes
     do i=1,this%lattice%nrec ! Loop on the number of atoms to be treat self-consistently
       do l=1,9  ! Loop on the orbital l
         write(123,*) 'orbital', l
         do ll=1,llmax ! Loop on the recursion steps
           write(123,*) this%a(ll,l,i,1), this%b2(ll,l,i,1)
         end do
       end do
     end do
    end subroutine


    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Obtain the optmal values of the terminators ainf and binf by Pettifor's
    !termination
    !---------------------------------------------------------------------------
    subroutine bpopt(this,ll,a,rb,n,ainf,rbinf,ifail)
      class(recursion), intent(inout) :: this
      ! Input
      integer, intent(in) :: n, ll 
      real(rp), dimension(ll), intent(in) :: A,RB
      ! Output
      integer, intent(out) :: ifail
      real(rp), intent(out) :: ainf,rbinf
      ! Local variables    
      integer :: I,JITER,NDIME
      real(rp) :: BM,BMAX,BMIN,EPS
      real(rp), dimension(ll) :: AZ,RBZ

      NDIME=ll  
      IFAIL = 0
      EPS = 1.0d-05
      JITER = 0
      AINF = A(N)
      do
         JITER = JITER + 1
         AZ(1) = 0.5d0 * (A(1)-AINF)
         do I = 2,N-1
            AZ(I) = 0.5d0 * (A(I)-AINF)
            RBZ(I) = 0.5d0 * RB(I)
         end do
         AZ(N) = A(N) - AINF
         RBZ(N) = 1.0d0 / SQRT(2.0d0) * RB(N)
         call this%EMAMI(NDIME,AZ,RBZ,N,BMAX,BMIN)
         BM = BMAX + BMIN
         BM = ABS(BM)
         AINF = AINF + (BMAX+BMIN)
         if (BM <= EPS) then
            exit
            !    elseif (JITER > 30) then
         elseif (JITER > 300) then
            IFAIL = 1
            exit
         end if
      end do
      RBINF = (BMAX-BMIN) / 2.0d0
      ! write(700,*) AINF,RBINF
    end subroutine bpopt

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Obtain the max. and the min. eigenvalues of symmetric tridiagonal matrix
    !by bisection method
    !---------------------------------------------------------------------------
    subroutine emami(this,nl,as,bs,n,emax,emin)
      class(recursion), intent(inout) :: this
      ! Input
      integer, intent(in) :: N, NL
      real(rp), dimension(NL), intent(in) :: AS,BS
      ! Output
      real(rp), intent(out) :: EMAX,EMIN
      ! Local Variables
      integer :: I,ISTOP,NUM
      real(rp) :: DELE,E,E1,E2,EMAX0,EMIN0,EPS,P,RELFEH,X1,X2
      real(rp), dimension(NL) :: A,B

      EMAX0 = -1.0d6
      EMIN0 = +1.0d6
      do I = 1,N
         A(I) = AS(I)
         B(I) = BS(I)
      end do
      B(1) = 0.0d0
      B(N+1) = 0.0d0
      do I = 1,N
         X1 = A(I) + ABS(B(I)) + ABS(B(I+1))
         X2 = A(I) - ABS(B(I)) - ABS(B(I+1))
         if (EMAX0 <= X1) then
            EMAX0 = X1
         end if
         if (EMIN0 > X2) then
            EMIN0 = X2
         end if
      end do
      RELFEH = 2.d0**(-39)
      EPS = 1.0d-6
      !C....CALCULATION OF EMAX                                               
      ISTOP = 0
      EMAX = EMAX0
      EMIN = EMIN0
      do
         E = (EMAX+EMIN) / 2.0d0
         ISTOP = ISTOP + 1
         if (ISTOP > 50) goto 1000
         NUM = 0
         P = A(1) - E
         if (P < 0.d0) then
            NUM = NUM + 1
         end if
         do I = 2,N
            if (P == 0.0) then
               P = (A(I)-E) - ABS(B(I))/RELFEH
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            else
               P = (A(I)-E) - B(I)**2/P
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            end if
         end do
         if (NUM == N) then
            EMAX = E
         end if
         if (NUM < N) then
            EMIN = E
         end if
         DELE = (EMAX-EMIN) / ((EMAX+EMIN)/2.d0)
         DELE = ABS(DELE)
         if (DELE <= EPS) exit
      end do
      E1 = E
      !.....CALCULATION ON EMIN                                               
      ISTOP = 0
      EMAX = E1
      EMIN = EMIN0
      do
         E = (EMAX+EMIN) / 2.d0
         ISTOP = ISTOP + 1
         if (ISTOP > 50) goto 1000
         NUM = 0
         P = A(1) - E
         if (P < 0.d0) then
            NUM = NUM + 1
         end if
         do I = 2,N
            if (P == 0.d0) then
               P = (A(I)-E) - ABS(B(I))/RELFEH
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            else
               P = (A(I)-E) - B(I)**2/P
               if (P < 0.d0) then
                  NUM = NUM + 1
               end if
            end if
         end do
         if (NUM == 0) then
            EMIN = E
         end if
         if (NUM > 0) then
            EMAX = E
         end if
         DELE = (EMAX-EMIN) / ((EMAX+EMIN)/2.d0)
         DELE = ABS(DELE)
         if (DELE <= EPS) exit
      end do
      E2 = E
      !....                                                                   
      EMAX = E1
      EMIN = E2
      return
      1000 continue
      !1000 write (6,10000)
      return
      ! 
      ! ... Format Declarations ...
      ! 
      10000 format (" ","NON-CONVERGE IN EMAMI")
    end subroutine emami

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this,full)
      class(recursion) :: this
      logical, intent(in), optional :: full
      allocate(this%a(max(this%lattice%control%llsp,this%lattice%control%lld),&
                         &18,this%lattice%nrec,3))
      allocate(this%atemp(max(this%lattice%control%llsp,this%lattice%control%lld)))
      allocate(this%b2temp(max(this%lattice%control%llsp,this%lattice%control%lld)))              
      allocate(this%b2(max(this%lattice%control%llsp,this%lattice%control%lld),&
                         &18,this%lattice%nrec,3))
      allocate(this%izero(0:this%lattice%kk),this%izeroll(0:this%lattice%kk,this%lattice%control%lld+1))

      allocate(this%psi(18,this%lattice%kk),this%pmn(18,this%lattice%kk))
      allocate(this%psi1(18,18,this%lattice%kk),this%psi2(18,18,this%lattice%kk),this%psi0(18,18,this%lattice%kk))
      allocate(this%v(18,this%lattice%kk))
      allocate(this%mu_n(this%lattice%nrec,(2*this%lattice%control%lld)+2,18,18),&
               this%mu_ng(this%lattice%nrec,(2*this%lattice%control%lld)+2,18,18))

      this%v(:,:) = 0.0d0
      this%psi(:,:) = 0.0d0
      this%psi1(:,:,:) = 0.0d0
      this%psi2(:,:,:) = 0.0d0
      this%pmn(:,:) = 0.0d0
      this%mu_n(:,:,:,:) = 0.0d0
      this%mu_ng(:,:,:,:) = 0.0d0
      this%izero(:) = 0
      this%izeroll(:,:) = 0
      this%a(:,:,:,:) = 0.0d0
      this%b2(:,:,:,:) = 0.0d0
      this%atemp(:) = 0.0d0
      this%b2temp(:) = 0.0d0

      if (present(full) .and. full) then
        if (associated(this%hamiltonian)) call this%hamiltonian%restore_to_default()
        if (associated(this%lattice)) call this%lattice%restore_to_default()
        if (associated(this%en)) call this%en%restore_to_default()
      endif

    end subroutine restore_to_default

end module recursion_mod
  
