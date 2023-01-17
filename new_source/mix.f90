!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Mix
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
!> lorem ipsum
!------------------------------------------------------------------------------


module mix_mod
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use control_mod
    use lattice_mod
    use charge_mod
    use symbolic_atom_mod, only: symbolic_atom
    use precision_mod, only: rp
    use string_mod
    use logger_mod, only: g_logger
    implicit none
  
    private
  
    !> Module's main structure
    type, public :: mix
      ! Lattice
      class(lattice), pointer :: lattice
      ! Charge
      class(charge), pointer :: charge
      ! Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Description
      !> 
      !> Description
      integer :: var
      !> Variables to save parameters for mixing
      real(rp), dimension(:,:), allocatable :: qia, qia_new, qia_old, qiaprev
      !> Mixing parameter
      real(rp) :: beta
      !> Difference between interactions
      real(rp) :: delta
    contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state
      procedure :: save_to
      procedure :: mixpq
      final :: destructor
    end type mix
  
    interface mix
      procedure :: constructor
    end interface mix
  
  contains
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Constructor
    !
    !> @param[in] lattice_obj pointer
    !> @return type(mix)
    !---------------------------------------------------------------------------
    function constructor(lattice_obj,charge_obj) result(obj)
      type(mix) :: obj
      type(lattice), target, intent(in) :: lattice_obj
      type(charge), target, intent(in) :: charge_obj 

      obj%lattice => lattice_obj
      obj%charge => charge_obj
      obj%symbolic_atom => lattice_obj%symbolic_atoms

      call obj%restore_to_default()
      call obj%build_from_file()
    end function constructor
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Destructor
    !---------------------------------------------------------------------------
    subroutine destructor(this)
      type(mix) :: this
    end subroutine destructor
  
    ! Member functions
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Read parameters from input file
    !
    !> @param[in] fname Namelist file
    !---------------------------------------------------------------------------
    subroutine build_from_file(this)
      class(mix),intent(inout) :: this
      character(len=sl) :: fname 
      integer :: var
      real(rp) :: beta
      ! variables associated with the reading processes
      integer :: iostatus, funit



      namelist /mix/ var, beta
  
      fname = this%lattice%control%fname 
      var = this%var
      beta = this%beta

      open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
      if(iostatus /= 0) then
        write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
        error stop
      endif
  
      read(funit,nml=mix,iostat=iostatus)
      if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
        write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
        write(error_unit,'("iostatus = ",I0)') iostatus
      endif
      close(funit)
  
      this%var = 0
      this%beta = beta  
    end subroutine build_from_file
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
      implicit none
      class(mix), intent(inout) :: this
  
      allocate(this%qia(this%lattice%nrec,18),this%qia_new(this%lattice%nrec,18),this%qia_old(this%lattice%nrec,18))
      allocate(this%qiaprev(this%lattice%nrec,18))

      this%var = 0
      this%beta = 0.1d0 
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
      class(mix), intent(in) :: this
  
      integer,intent(in),optional :: unit
      character(len=*),intent(in),optional :: file
      integer :: newunit
  
      integer :: var
  
      namelist /mix/ var
  
      var = this%var
  
      if(present(unit) .and. present(file)) then
        call g_logger%fatal('Argument error: both unit and file are present',__FILE__,__LINE__)
      else if(present(unit)) then
          write(unit,nml=mix)
      else if(present(file)) then
          open(unit=newunit,file=file)
          write(newunit,nml=mix)
          close(newunit)
      else
          write(*,nml=mix)
      endif
      close(newunit)
    end subroutine print_state

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Save potential parameters and band moments.
    !>
    !> Save potential parameters and band moments. Save the parameters read from the &par namelist into qia variables. To be used to mix later.
    !> @param[in] dummy variable of type mix
    !> @return type(mix)
    !---------------------------------------------------------------------------
    subroutine save_to(this,whereto)
      class(mix), intent(inout) :: this
      character(len=*), intent(in) :: whereto
      integer :: i, it


      select case(whereto)

      case('old')
        do it=1,this%lattice%nrec
          do I = 1,3
            this%qia_old(IT,I) = this%symbolic_atom(it)%potential%ql(1,I-1,1)
            this%qia_old(IT,I+6) = this%symbolic_atom(it)%potential%ql(3,I-1,1)
            this%qia_old(IT,I+12) = this%symbolic_atom(it)%potential%pl(I-1,1) !- (0.5d0 + INT(PL(I,1)))
            !QI_OLD(IT,I+12) = ENU(I,1)
            this%qia_old(IT,I+3) = this%symbolic_atom(it)%potential%ql(1,I-1,2)
            this%qia_old(IT,I+9) = this%symbolic_atom(it)%potential%ql(3,I-1,2)
            this%qia_old(IT,I+15) = this%symbolic_atom(it)%potential%pl(I-1,2) !- (0.5d0 + INT(PL(I,2)))
            !QI_OLD(IT,I+15) = ENU(I,2)
          end do
        end do
      case('new')
        do it=1,this%lattice%nrec
          do I = 1,3
            this%qia_new(it,i) = this%symbolic_atom(it)%potential%ql(1,i-1,1)
            this%qia_new(it,i+6) = this%symbolic_atom(it)%potential%ql(3,i-1,1)
            this%qia_new(it,i+12) = this%symbolic_atom(it)%potential%pl(i-1,1) !- (0.5d0 + INT(PL(I,1)))
            !QI_OLD(IT,I+12) = ENU(I,1)
            this%qia_new(it,i+3) = this%symbolic_atom(it)%potential%ql(1,i-1,2)
            this%qia_new(it,i+9) = this%symbolic_atom(it)%potential%ql(3,i-1,2)
            this%qia_new(it,i+15) = this%symbolic_atom(it)%potential%pl(i-1,2) !- (0.5d0 + INT(PL(I,2)))
            !QI_OLD(IT,I+15) = ENU(I,2)
          end do
        end do

      case('prev')
        do it=1,this%lattice%nrec
          do I = 1,3
            this%qiaprev(it,i) = this%symbolic_atom(it)%potential%ql(1,i-1,1)
            this%qiaprev(it,i+6) = this%symbolic_atom(it)%potential%ql(3,i-1,1)
            this%qiaprev(it,i+12) = this%symbolic_atom(it)%potential%pl(i-1,1) !- (0.5d0 + INT(PL(I,1)))
            !QI_OLD(IT,I+12) = ENU(I,1)
            this%qiaprev(it,i+3) = this%symbolic_atom(it)%potential%ql(1,i-1,2)
            this%qiaprev(it,i+9) = this%symbolic_atom(it)%potential%ql(3,i-1,2)
            this%qiaprev(it,i+15) = this%symbolic_atom(it)%potential%pl(i-1,2) !- (0.5d0 + INT(PL(I,2)))
            !QI_OLD(IT,I+15) = ENU(I,2)
          end do
        end do

      case('current')
        do it=1,this%lattice%nrec
          do I = 1,3
            this%symbolic_atom(it)%potential%ql(1,i-1,1) = this%qia(it,i)
            this%symbolic_atom(it)%potential%ql(3,i-1,1) = this%qia(it,i+6)
            this%symbolic_atom(it)%potential%pl(i-1,1) = this%qia(it,i+12)
            this%symbolic_atom(it)%potential%ql(1,i-1,2) = this%qia(it,i+3)
            this%symbolic_atom(it)%potential%ql(3,i-1,2) = this%qia(it,i+9) 
            this%symbolic_atom(it)%potential%pl(i-1,2) = this%qia(it,i+15) 
          end do
        end do
      end select
    end subroutine save_to 

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Mix the parameters QL and PL. Can be linear, andersen or broyden.
    !>
    !> Mix the parameters QL and PL. Can be linear, andersen or broyden.
    !> @param[in] dummy variable of type mix
    !> @param[in] old and new potential variables qia_old and qia_new, respectively.
    !> @return type(mix): qia, the mix potential variable.
    !---------------------------------------------------------------------------
    subroutine mixpq(this,qia_old,qia_new)
      class(mix), intent(inout) :: this
      real(rp), dimension(this%lattice%nrec,18), intent(in) :: qia_old, qia_new
      integer :: ia ! Atom index

      this%qia(:,:) = (1-this%beta)*qia_old(:,:) + this%beta*qia_new(:,:)

      this%charge%dq(:) = 0.0d0
      do ia=1,this%lattice%nrec
        this%charge%trq = 0.0d0 
        this%charge%cht = 0.0d0
        this%charge%cht = sum(this%qia_new(ia,1:6)) - this%symbolic_atom(ia)%element%valence
        this%charge%dq(ia) = sum(this%qia(ia,1:6)) - this%symbolic_atom(ia)%element%valence
        this%charge%trq = this%charge%dq(ia)
        call g_logger%info('charge transfer at atom '//int2str(ia),__FILE__,__LINE__)
        call g_logger%info(real2str(this%charge%cht)//' '//real2str(this%charge%trq)//' '//real2str(this%charge%cht-this%charge%trq),__FILE__,__LINE__)
      end do
 
      this%delta = sqrt(sum((this%qia_old(:,:)-this%qia_new(:,:))**2))/6.0d0/this%lattice%nrec
    end subroutine mixpq
end module mix_mod
