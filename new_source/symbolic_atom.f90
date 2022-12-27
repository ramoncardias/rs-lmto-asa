!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Symbolic Atom
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
!> Module to hold symbolic_atom's parameters
!------------------------------------------------------------------------------

module symbolic_atom_mod
    use, intrinsic :: iso_fortran_env, only: error_unit
    use element_mod, only: element
    use potential_mod, only: potential
    use string_mod, only: sl, replace, endswith
    use precision_mod, only: rp
    implicit none

    private

    ! public functions
    public :: array_of_symbolic_atoms
    public :: print_state, print_state_full
    public :: save_state
    public :: load_state

    type, public :: symbolic_atom
        !> Mixture occupation in self-consistent calculation. Default: 0.01.
        !> 
        !> Mixture occupation in self-consistent calculation.
        !> 
        !> Default: 0.01.
        real(rp) :: mix

        !> Spin-(up/down) occupation in self-consistent calculation. Default: 0.05.
        !> 
        !> Spin-(up/down) occupation in self-consistent calculation.
        !> 
        !> Default: 0.05.
        real(rp) :: mixmag

        !> Specify the number of rigid band calculation used.
        !>
        !> Specify the number of rigid band calculation used
        integer :: rb
        ! TODO
        ! From common_cnstr
        real(rp), dimension(3) :: mag_cfield, mag_cfield_diff
        real(rp) :: chg_cfield, chg_cfield_diff
        real(rp) :: chg_con_val, mag_con_val
        ! TODO
        real(rp) :: a
        ! TODO (Charge?)
        real(rp) :: dq, qc, qv
        
        type(element) :: element
        type(potential) :: potential
    contains
        procedure :: build_from_file
        procedure :: restore_to_default
        procedure :: print_state => symbolic_atom_print_state
        procedure :: print_state_full => symbolic_atom_print_state_full
        procedure :: mesh_grid_size
        procedure :: B
        procedure :: rho0
    end type symbolic_atom

    interface symbolic_atom
        procedure :: constructor
    end interface symbolic_atom

    interface array_of_symbolic_atoms
        procedure :: array_of_symbolic_atoms_from_memory
        procedure :: array_of_symbolic_atoms_from_file
    end interface array_of_symbolic_atoms

    integer, parameter :: min_mesh_grid_size = 25

contains

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Constructor
    !
    !> @param[in] symbolic_atom Namelist file in database
    !> @param[in] database Directory to database files with 'symbolic_atom' namelist
    !> @return type(control)
    !---------------------------------------------------------------------------
    function constructor(label,database) result(obj)
        type(symbolic_atom) :: obj
        character(len=*), intent(in) :: label
        character(len=*), intent(in), optional :: database
        integer :: i

        call obj%restore_to_default()
        obj%element = element(label,database)
        obj%potential = potential(label,database)

        ! Initial guess for the Pl's potential
        !  do i=1,2 ! loop on spin channel
        !    obj%potential%pl(1,i) = obj%element%num_quant_s + 0.5d0
        !    obj%potential%pl(2,i) = obj%element%num_quant_p + 0.5d0
        !    obj%potential%pl(3,i) = obj%element%num_quant_d + 0.5d0
        !  end do
    end function constructor
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Read parameters from input file
    !
    !> @param[in] fname Namelist file
    !---------------------------------------------------------------------------
    subroutine build_from_file(this, fname)
        class(symbolic_atom), intent(out) :: this
        character(len=*), intent(in) :: fname
        call this%element%build_from_file(fname)
        call this%potential%build_from_file(fname)
    end subroutine build_from_file
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
        implicit none
        class(symbolic_atom), intent(out) :: this
        this%mix = -1
        this%mixmag = -1
        this%rb = -1
        this%a = 0.02d0

        call this%element%restore_to_default()
        call this%potential%restore_to_default()
    end subroutine restore_to_default

    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Print class members values in namelist format 
    !>
    !> Print class members values in namelist format 
    !---------------------------------------------------------------------------
    subroutine symbolic_atom_print_state(this,unit,file)
        class(symbolic_atom) :: this
        integer,optional :: unit
        character(len=*),optional :: file
        integer :: newunit
        if(present(unit) .and. present(file)) then
            write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
            error stop
        else if(present(unit)) then
            call this%element%print_state(unit=unit)
            call this%potential%print_state(unit=unit)
        else if(present(file)) then
            open(newunit=newunit,file=file,action='write')
            call this%element%print_state(unit=newunit)
            call this%potential%print_state(unit=newunit)
            close(newunit)
        else
            call this%element%print_state()
            call this%potential%print_state()
        endif
    end subroutine symbolic_atom_print_state

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Print class members values in namelist format 
    !>
    !> Print class members values in namelist format 
    !---------------------------------------------------------------------------
    subroutine symbolic_atom_print_state_full(this,unit,file)
        class(symbolic_atom) :: this
        integer,optional :: unit
        character(len=*),optional :: file
        integer :: newunit
        if(present(unit) .and. present(file)) then
            write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
            error stop
        else if(present(unit)) then
            call this%element%print_state_full(unit=unit)
            call this%potential%print_state_full(unit=unit)
        else if(present(file)) then
            open(newunit=newunit,file=file,action='write')
            call this%element%print_state_full(unit=newunit)
            call this%potential%print_state_full(unit=newunit)
            close(newunit)
        else
            call this%element%print_state_full()
            call this%potential%print_state_full()
        endif
    end subroutine symbolic_atom_print_state_full

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> TODO
    !>
    !> TODO
    !---------------------------------------------------------------------------
    function mesh_grid_size(this)
        class(symbolic_atom), intent(in) :: this
        integer :: mesh_grid_size
        real(rp) :: B, Z
        intrinsic EXP, LOG, MAX

        Z = this%element%atomic_number
        B = 1.d0 / (Z+Z+1.d0)
        mesh_grid_size = MAX(min_mesh_grid_size,int(((.5d0+LOG(1.d0+this%potential%ws_r/B)/this%A) * 2.d0-1)/2)*2 + 1)
    end function mesh_grid_size

    !---------------------------------------------------------------------------
    !> DESCRIPTION:
    !> @breif
    !> TODO Calculates the quantity B whatever it is
    !---------------------------------------------------------------------------
    function B(this)
        class(symbolic_atom), intent(in) :: this
        real(rp) :: B
        B = this%potential%ws_r / (EXP(this%A*this%mesh_grid_size()-this%A)-1.d0)
    end function B

    !---------------------------------------------------------------------------
    !> DESCRIPTION:
    !> @breif
    !> TODO Calculates the initial guess for charge density
    !---------------------------------------------------------------------------
    function rho0(this,nsp)
        class(symbolic_atom), intent(in) :: this
        integer, intent(in), optional :: nsp
        real(rp), dimension(:,:), allocatable :: rho0
        real(rp) :: b, ea, sum, rpb, r, ro, fac
        integer :: ir, nr
        integer :: nsp_
        nsp_ = merge(nsp,1,present(nsp))
        nr = this%mesh_grid_size()
        allocate(rho0(nr,nsp_))
        b = this%B()
        ea = exp(this%a)
        sum = 0.d0
        rpb = b
        do ir = 1,nr
            r = rpb - b
            ro = exp(-5.*r) * r * r
            rho0(ir,1) = ro
            sum = sum + this%a*rpb*ro
            rpb = rpb * ea
        end do
        fac = this%element%atomic_number / (sum*2)
        do ir = 1,nr
            rho0(ir,1) = rho0(ir,1) * fac
            rho0(ir,2) = rho0(ir,1)
            write(101,*) ir, rho0(ir,1)
            write(102,*) ir, rho0(ir,2)
        end do
    end function rho0

    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief Build an array of symbolic_atoms
    !
    !> @param[in] symbolic_atom List of labels in database to build a symbolic_atom
    !> @param[in] database Directory with files used as database
    !> @return type(symbolic_atom), dimension(:), allocatable
    !---------------------------------------------------------------------------
    function array_of_symbolic_atoms_from_memory(symbolic_atoms,database)
        type(symbolic_atom), dimension(:), allocatable :: array_of_symbolic_atoms_from_memory
        character(len=*), dimension(:), intent(in) :: symbolic_atoms
        character(len=*), intent(in), optional :: database
        integer :: i, j
        allocate(array_of_symbolic_atoms_from_memory(size(symbolic_atoms)))
        if(present(database)) then
            do i=1, size(symbolic_atoms)
                array_of_symbolic_atoms_from_memory(i) = symbolic_atom(symbolic_atoms(i),database)
            enddo
        else
            do i=1, size(symbolic_atoms)
                array_of_symbolic_atoms_from_memory(i) = symbolic_atom(symbolic_atoms(i))
            enddo
        endif
    end function array_of_symbolic_atoms_from_memory

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief Build an array of symbolic_atoms from namelist 'atoms' in fname
    !
    !> @param[in] fname File containg the namelist 'atoms'
    !> @param[in] size Size (integer) or array
    !> @return type(symbolic_atom), dimension(:), allocatable
    !---------------------------------------------------------------------------
    function array_of_symbolic_atoms_from_file(fname,size)
        use string_mod, only: sl
        type(symbolic_atom), dimension(:), allocatable :: array_of_symbolic_atoms_from_file
        integer, intent(in) :: size
        character(len=*), intent(in) :: fname

        ! Readable variables
        character(len=sl), dimension(:), allocatable :: label
        ! integer, dimension(:), allocatable :: label
        character(len=sl) :: database = './'

        ! variables associated with the reading processes
        integer :: iostatus, funit
        
        namelist /atoms/ database, label

        allocate(label(size))

        open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
        if(iostatus /= 0) then
            write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
            error stop
        endif
        
        read(funit,nml=atoms,iostat=iostatus)
        if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
            write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
            write(error_unit,'("iostatus = ",I0)') iostatus
              endif
        close(funit)
        
        array_of_symbolic_atoms_from_file = array_of_symbolic_atoms(label,database)

    end function array_of_symbolic_atoms_from_file

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Print class members values in namelist format 
    !>
    !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
    !> @param[in] unit File unit used to write namelist
    !> @param[in] file File name used to write namelist
    !> @param[in] suffix If provided, prints the state to file with the same symbol added to this suffix
    !---------------------------------------------------------------------------
    subroutine print_state(array,unit,file,suffix)
        type(symbolic_atom), dimension(:), intent(in) :: array
        integer,intent(in),optional :: unit
        character(len=*),intent(in),optional :: file, suffix
        character(len=sl) :: new_file
        integer :: i,newunit
        if((present(unit) .and. present(file)) .or. \
            (present(unit) .and. present(suffix)) .or. \
            (present(file) .and. present(suffix))) then
            write(error_unit,'("[",A,":",I0,"]: Argument error: both unit, file and suffix are present")') __FILE__,__LINE__
            error stop
        else if(present(file)) then
            open(newunit=newunit,file=file)
            do i=1, size(array)
                call array(i)%print_state(unit=newunit)
            enddo
            close(newunit)
        else if(present(suffix)) then
            do i=1, size(array)
                new_file = trim(array(i)%element%symbol)
                if(endswith(new_file,".nml")) then
                    call replace(new_file,".nml","")
                endif
                new_file = trim(new_file)//trim(suffix)//".nml"
                call array(i)%print_state(file=new_file)
            enddo
        else
            do i=1, size(array)
                call array(i)%print_state(unit=unit)
            enddo
        endif
    end subroutine print_state

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Print class members values in namelist format 
    !>
    !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
    !> @param[in] unit File unit used to write namelist
    !> @param[in] file File name used to write namelist
    !> @param[in] suffix If provided, prints the state to file with the same symbol added to this suffix
    !---------------------------------------------------------------------------
    subroutine print_state_full(array,unit,file,suffix)
        type(symbolic_atom), dimension(:), intent(in) :: array
        integer,intent(in),optional :: unit
        character(len=*),intent(in),optional :: file, suffix
        character(len=sl) :: new_file
        integer :: i,newunit
        if((present(unit) .and. present(file)) .or. \
            (present(unit) .and. present(suffix)) .or. \
            (present(file) .and. present(suffix))) then
            write(error_unit,'("[",A,":",I0,"]: Argument error: both unit, file and suffix are present")') __FILE__,__LINE__
            error stop
        else if(present(file)) then
            open(newunit=newunit,file=file,action='write')
            do i=1, size(array)
                call array(i)%print_state_full(unit=newunit)
            enddo
            close(newunit)
        else if(present(suffix)) then
            do i=1, size(array)
                new_file = array(i)%element%symbol
                if(endswith(new_file,".nml")) then
                    call replace(new_file,".nml","")
                endif
                new_file = trim(new_file)//trim(suffix)//".nml"
                call array(i)%print_state_full(file=new_file)
            enddo
        else
            do i=1, size(array)
                call array(i)%print_state_full(unit=unit)
            enddo
        endif
    end subroutine print_state_full

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Shortcut for print_state_full with suffix="_out" 
    !>
    !> Shortcut for print_state_full with suffix="_out" 
    !---------------------------------------------------------------------------
    subroutine save_state(array)
        type(symbolic_atom), dimension(:), intent(in) :: array
        call print_state(array,suffix="_out")
    end subroutine save_state

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Checks the existence of checkpoint files (those with fname = symbol//"_out.nml" )
    !>
    !> Checks the existence of checkpoint files (those with fname = symbol//"_out.nml" )
    !---------------------------------------------------------------------------
    subroutine load_state(array)
        type(symbolic_atom), dimension(:), intent(inout) :: array
        character(len=sl) :: fname
        integer :: i, array_size
        logical :: file_exists
        
        array_size = size(array)
        do i=1, array_size
            fname = trim(array(i)%element%symbol)//"_out.nml"
            INQUIRE(FILE=fname, EXIST=file_exists)
            if(file_exists) then
                call array(i)%build_from_file(fname)
            endif
        enddo
    end subroutine load_state

    
end module symbolic_atom_mod
