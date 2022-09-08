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
!> Module to handle Greens functions calculation and related routines
!------------------------------------------------------------------------------

module green_mod

    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use self_mod
    use control_mod
    use lattice_mod
    use symbolic_atom_mod
    use recursion_mod
    use density_of_states_mod
    use precision_mod, only : rp
    use math_mod
    implicit none
  
    private

    !> Module's main structure
    type, public :: green
      ! General variables
      !> Recursion
      type(recursion), pointer :: recursion
      !> Lattice
      type(lattice), pointer :: lattice
      !> Self
      type(self), pointer :: self
      !> Symbolic atom
      type(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Density of states
      type(dos), pointer :: dos
      !> Control
      type(control), pointer :: control
      !> Onsite Green Function
      complex(rp), dimension(:,:,:,:), allocatable :: g0
    contains
      procedure :: sgreen
      procedure :: restore_to_default
      final :: destructor
    end type
   
    interface green
      procedure :: constructor
    end interface green

contains

  ! DESCRIPTION:
  !> @brief
  !> Constructor
  !
  !> @param[in] fname Namelist file
  !> @return type(green)
  !---------------------------------------------------------------------------
  function constructor(dos_obj) result(obj)
    type(green) :: obj
    type(dos), target, intent(in) :: dos_obj

    obj%dos => dos_obj
    obj%recursion => dos_obj%recursion
    obj%self => dos_obj%self
    obj%symbolic_atom => dos_obj%recursion%hamiltonian%symbolic_atom
    obj%lattice => dos_obj%recursion%lattice
    obj%control => dos_obj%recursion%lattice%control

    call obj%restore_to_default()
  end function constructor

  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Destructor
  !---------------------------------------------------------------------------
  subroutine destructor(this)
    type(green) :: this
    if(allocated(this%g0)) deallocate(this%g0)
  end subroutine destructor

  ! Member functions
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
      class(green) :: this
      allocate(this%g0(this%self%channels_ldos+10,18,18,this%lattice%ntype))

      this%g0(:,:,:,:) = (0.0d0,0.0d0)
    end subroutine restore_to_default


  !---------------------------------------------------------------------------
  ! DESCRIPTION:
  !> @brief
  !> Calculates the Greens function 
  !---------------------------------------------------------------------------
  subroutine sgreen(this)
!this,bdos,g0,e,nv,ldim,ll,na,nmdir,mom,ef)
    class(green) :: this
    ! Input
    !integer, intent(in) :: nv, ldim, ll, na, nmdir
    !real(rp), intent(in) :: ef
    ! Output
    !real(rp), dimension(nv), intent(inout) :: e
    !real(rp), dimension(nv,ldim,ldim,na), intent(out) :: bdos
    !complex(rp), dimension(nv,ldim,ldim,na), intent(out) :: g0
    !real(rp), dimension(na,3), intent(inout) :: mom
    ! Local variables
    integer :: ia,ja,mdir,nw, ll_t, ie, j, i
    real(rp), dimension(this%self%channels_ldos+10,this%lattice%ntype) :: dx,dy,dz
    real(rp), dimension(18,this%self%channels_ldos+10) :: doso
    real(rp), dimension(18,this%self%channels_ldos+10) :: dmag,dnmag
    complex(rp) :: dfac, sfac, impi
    complex(rp),dimension(4) :: gspinor
    complex(rp), dimension(3) :: gmask, lmask
    complex(rp), dimension(2,3) :: gfac
    integer, dimension(4,3) :: goff
    real(rp),dimension(this%control%lld,this%control%lld)  :: Sm
    real(rp),dimension(18)  :: q_int
    real(rp) :: e_start,e_stop

    impi=(pi,0.0d0)
    gfac(1,1)=1.0d0; gfac(1,2)=-i_unit; gfac(1,3)=1.0d0
    gfac(2,1)=1.0d0; gfac(2,2)= i_unit; gfac(2,3)=-1.0d0
    gmask(1)=0.0d0;gmask(2)=0.0d0;gmask(3)=1.0d0
!   if(lrot) then
!      lmask(1)=0.0d0;lmask(2)=0.0d0;lmask(3)=1.0d0
!   else
      lmask(1)=1.0d0/3.0d0;lmask(2)=1.0d0/3.0d0;lmask(3)=1.0d0/3.0d0
!   end if

    goff(1,1)=0;goff(2,1)=9;goff(1,2)=0;goff(2,2)=9;goff(1,3)=0;goff(2,3)=0
    goff(3,1)=9;goff(4,1)=0;goff(3,2)=9;goff(4,2)=0;goff(3,3)=9;goff(4,3)=9
    dfac= i_unit*impi/(2.0d0,0.0d0)
    dx=0;dy=0;dz=0
    nw=10*this%lattice%ntype*this%control%lld
    ll_t=this%control%lld

!    bdos=0.0d0
!    g0=0.0d0
    do ia = 1, this%lattice%nrec
      do mdir = 1, this%control%nmdir
        doso = 0.0d0
        call this%dos%density(doso,ia,mdir)
        if(this%control%nmdir==1) then
          do ie=1,this%self%channels_ldos+10
            do j=1,18
              this%g0(ie,j,j,ia)=-i_unit*doso(j,ie)*impi
            end do
            write(300+ia,*) this%self%ene(ie), sum(doso(1:9,ie)), sum(doso(10:18,ie))
          end do
        else  
          do ie=1,this%self%channels_ldos+10
            do j=1,9
              ! Charge, from main direction.. (not z-component)
              this%g0(ie,j,j,ia)    =this%g0(ie,j,j,ia)    -(doso(j,ie)+doso(j+9,ie))*dfac*lmask(mdir)!*1.0d0 /3.0d0 !mom(ja,mdir)**2
              this%g0(ie,j+9,j+9,ia)=this%g0(ie,j+9,j+9,ia)-(doso(j,ie)+doso(j+9,ie))*dfac*lmask(mdir)!*1.0d0 /3.0d0 !mom(ja,mdir)**2
              ! Spin dependent part
              this%g0(ie,j+goff(1,mdir),j+goff(2,mdir),ia) = &
                 this%g0(ie,j+goff(1,mdir),j+goff(2,mdir),ia)-(doso(j,ie)-doso(j+9,ie))*gfac(1,mdir)*dfac
  
              this%g0(ie,j+goff(3,mdir),j+goff(4,mdir),ia) = &
                 this%g0(ie,j+goff(3,mdir),j+goff(4,mdir),ia)-(doso(j,ie)-doso(j+9,ie))*gfac(2,mdir)*dfac
            end do
          end do
        end if
      end do
    end do
  end subroutine sgreen
end module
