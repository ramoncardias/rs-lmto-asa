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
    use symbolic_atom_mod, only: symbolic_atom
    use precision_mod, only: rp
    use string_mod
    implicit none
  
    private
  
    !> Module's main structure
    type, public :: mix
      ! Lattice
      class(lattice), pointer :: lattice
      ! Symbolic atom
      class(symbolic_atom), dimension(:), pointer :: symbolic_atom
      !> Description
      !> 
      !> Description
      integer :: var
      !> Variables to save parameters for mixing
      real(rp), dimension(:,:), allocatable :: qia, qia_new, qia_old
    contains
      procedure :: build_from_file
      procedure :: restore_to_default
      procedure :: print_state
      procedure :: save_to_old
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
    function constructor(lattice_obj) result(obj)
      type(mix) :: obj
      type(lattice), target, intent(in) :: lattice_obj
  

      obj%lattice => lattice_obj
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

      ! variables associated with the reading processes
      integer :: iostatus, funit



      namelist /mix/ var
  
      fname = this%lattice%control%fname 
      var = this%var
  
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
  
    end subroutine build_from_file
  
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
      implicit none
      class(mix), intent(inout) :: this
  
      write(*,*) this%lattice%nrec
      allocate(this%qia(this%lattice%nrec,18),this%qia_new(this%lattice%nrec,18),this%qia_old(this%lattice%nrec,18))

      this%var = 0
  

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
        write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
        error stop
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
    !> Save potential parameters and band moments. Save the parameters read from the &par namelist into qia_old variable. To be used to mix later.
    !> @param[in] dummy variable of type mix
    !> @return type(mix)
    !---------------------------------------------------------------------------
    subroutine save_to_old(this)
      class(mix), intent(inout) :: this
      integer :: i, it

      do it=1,this%lattice%nrec
        do I = 1,3
          this%QIA_OLD(IT,I) = this%symbolic_atom(it)%potential%QL(1,I-1,1)
          this%QIA_OLD(IT,I+6) = this%symbolic_atom(it)%potential%QL(3,I-1,1)
          this%QIA_OLD(IT,I+12) = this%symbolic_atom(it)%potential%PL(I-1,1) !- (0.5d0 + INT(PL(I,1)))
          !QI_OLD(IT,I+12) = ENU(I,1)
          this%QIA_OLD(IT,I+3) = this%symbolic_atom(it)%potential%QL(1,I-1,2)
          this%QIA_OLD(IT,I+9) = this%symbolic_atom(it)%potential%QL(3,I-1,2)
          this%QIA_OLD(IT,I+15) = this%symbolic_atom(it)%potential%PL(I-1,2) !- (0.5d0 + INT(PL(I,2)))
          !QI_OLD(IT,I+15) = ENU(I,2)
        end do
      end do
    end subroutine save_to_old 
end module mix_mod
