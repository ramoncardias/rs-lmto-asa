!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Element
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
!> Module to hold element's parameters
!------------------------------------------------------------------------------

module element_mod
    use, intrinsic :: iso_fortran_env, only: error_unit
    use precision_mod, only: rp
    use globals_mod, only: GLOBAL_DATABASE_FOLDER, GLOBAL_CHAR_SIZE
    use string_mod, only: path_join, sl
    use logger_mod, only: g_logger
    implicit none

    private

    ! public functions
    public :: array_of_elements

    type, public :: element
        character(len=10) :: symbol
        !> Element related variables
        integer :: f_core, num_quant_s, num_quant_p, num_quant_d
        real(rp) :: atomic_number, core, valence
    contains
        procedure :: build_from_file
        procedure :: restore_to_default
        procedure :: print_state
        procedure :: print_state_full
        final :: destructor
    end type element

    interface element
        procedure :: constructor
    end interface element

contains

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Constructor
    !
    !> @param[in] element Namelist file in database
    !> @param[in] database Directory to database files with 'element' namelist
    !> @return type(control)
    !---------------------------------------------------------------------------
    function constructor(label,database) result(obj)
        type(element) :: obj
        character(len=*), intent(in) :: label
        character(len=*), intent(in), optional :: database
        character(len=3*sl) :: path_to_file
        character(len=sl), dimension(2) :: lst_path_to_file

        call obj%restore_to_default()
        if(present(database)) then
            lst_path_to_file(1) = database
        else
            lst_path_to_file(1) = './'
        endif
        lst_path_to_file(2) = trim(label)//'.nml'
        path_to_file = path_join(lst_path_to_file)
        if(exists(path_to_file)) then
            call obj%build_from_file(path_to_file)
        else
            lst_path_to_file(1) = GLOBAL_DATABASE_FOLDER
            path_to_file = path_join(lst_path_to_file)
            if(exists(path_to_file)) then
                call obj%build_from_file(path_to_file)
            else
                call g_logger%fatal('Element '//trim(label)//' not found in any database',__FILE__,__LINE__)
            endif
        endif
    end function constructor
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Destructor
    !---------------------------------------------------------------------------
    subroutine destructor(this)
      type(element) :: this
    end subroutine destructor


    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Read parameters from input file
    !
    !> @param[in] fname Namelist file
    !---------------------------------------------------------------------------
    subroutine build_from_file(this, fname)
        class(element),intent(inout) :: this
        character(len=*), intent(in) :: fname

        ! Readable variables
        character(len=10) :: symbol
        integer :: f_core, num_quant_s, num_quant_p, num_quant_d
        real(rp) :: atomic_number, core, valence
        
        ! variables associated with the reading processes
        integer :: iostatus, funit

        namelist /element/ symbol, atomic_number, core, valence, f_core, num_quant_s, num_quant_p, num_quant_d
        ! Save previous values
        symbol = this%symbol
        atomic_number = this%atomic_number
        core = this%core
        valence = this%valence
        f_core = this%f_core
        num_quant_s = this%num_quant_s
        num_quant_p = this%num_quant_p
        num_quant_d = this%num_quant_d

        open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
        if(iostatus /= 0) then
            write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
            error stop
        endif
        
        read(funit,nml=element,iostat=iostatus)
        if(iostatus /= 0 .and. .not. IS_IOSTAT_END(iostatus)) then
            write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
            write(error_unit,'("iostatus = ",I0)') iostatus
        endif
        close(funit)

        ! Setting user values
        this%symbol = symbol
        this%atomic_number = atomic_number
        this%core = core
        this%valence = valence
        this%f_core = f_core
        this%num_quant_s = num_quant_s
        this%num_quant_p = num_quant_p
        this%num_quant_d = num_quant_d
    end subroutine build_from_file
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
        implicit none
        class(element), intent(out) :: this

        this%symbol = ''
        this%atomic_number = -1
        this%core = -1
        this%valence = -1
        this%f_core = -1
        this%num_quant_s = -1
        this%num_quant_p = -1
        this%num_quant_d = -1
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
        class(element), intent(in) :: this

        integer,intent(in),optional :: unit
        character(len=*),intent(in),optional :: file
        integer :: newunit
        
        character(len=10) :: symbol
        integer :: f_core, num_quant_s, num_quant_p, num_quant_d
        real(rp) :: atomic_number, core, valence

        namelist /element/ symbol, atomic_number, core, valence, f_core, num_quant_s, num_quant_p, num_quant_d

        symbol = this%symbol
        atomic_number = this%atomic_number
        core = this%core
        valence = this%valence
        f_core = this%f_core
        num_quant_s = this%num_quant_s
        num_quant_p = this%num_quant_p
        num_quant_d = this%num_quant_d

        if(present(unit) .and. present(file)) then
            write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
            error stop
        else if(present(unit)) then
            write(unit,nml=element)
        else if(present(file)) then
            open(unit=newunit,file=file)
            write(newunit,nml=element)
            close(newunit)
        else
            write(*,nml=element)
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
    !---------------------------------------------------------------------------
    subroutine print_state_full(this,unit,file)
        class(element), intent(in) :: this

        integer,intent(in),optional :: unit
        character(len=*),intent(in),optional :: file
        integer :: newunit

        character(len=10) :: symbol
        integer :: f_core, num_quant_s, num_quant_p, num_quant_d
        real(rp) :: atomic_number, core, valence

        namelist /element/ symbol, atomic_number, core, valence, f_core, num_quant_s, num_quant_p, num_quant_d

        symbol = this%symbol
        atomic_number = this%atomic_number
        core = this%core
        valence = this%valence
        f_core = this%f_core
        num_quant_s = this%num_quant_s
        num_quant_p = this%num_quant_p
        num_quant_d = this%num_quant_d

        if(present(unit) .and. present(file)) then
            write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
            error stop
        else if(present(unit)) then
            write(unit,nml=element)
        else if(present(file)) then
            open(unit=newunit,file=file)
            write(newunit,nml=element)
            close(newunit)
        else
            write(*,nml=element)
        endif

    end subroutine print_state_full

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief Build an array of elements
    !
    !> @param[in] element List of labels in database to build elements
    !> @param[in] database Directory to database files with 'element' namelist
    !> @return type(element), dimension(:), allocatable
    !---------------------------------------------------------------------------
    function array_of_elements(elements,database)
        type(element), dimension(:), allocatable :: array_of_elements
        character(len=*), dimension(:), intent(in) :: elements
        character(len=*), intent(in), optional :: database
        integer :: i
        allocate(array_of_elements(size(elements)))
        if(present(database)) then
            do i=1, size(elements)
                array_of_elements(i) = element(elements(i),database)
            enddo
        else
            do i=1, size(elements)
                array_of_elements(i) = element(elements(i))
            enddo
        endif
    end function array_of_elements

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Check whether 'fname' exists
    !
    !> @param[in] fname File to check
    !---------------------------------------------------------------------------
    function exists(fname)
        character(len=*), intent(in) :: fname
        logical :: exists
        exists = .False.
        INQUIRE(FILE=fname, EXIST=exists)
    end function exists
end module element_mod
