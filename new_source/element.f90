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
    implicit none

    private

    ! public functions
    public :: array_of_elements

    type, public :: element
        character(len=10) :: symbol
        !> Element related variables
        integer :: f_core, num_quant_s, num_quant_p, num_quant_d
        real(rp) :: atomic_number, core, valence
        !> Moments as defined in Eq. 48. of Phys. Rev. B 43, 9538 (1991).
        !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
        !> 2nd index 1 = spin-up, 2 = spin-dw
        real(rp), dimension(:,:), allocatable  :: q0, q1, q2
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
        call obj%restore_to_default()
        if(present(database)) then
            call obj%build_from_file(trim(database) // '/' // trim(label) // '.nml')
        else if(exists('./' // trim(label) // '.nml')) then
            call obj%build_from_file('./' // trim(label) // '.nml')
        else if(exists(GLOBAL_DATABASE_FOLDER // '/' // trim(label) // '.nml')) then
            call obj%build_from_file(GLOBAL_DATABASE_FOLDER // '/' // trim(label) // '.nml')
        else
            write(error_unit,'("[",A,":",I0,"]: Element ",A," not found in any database")') __FILE__,__LINE__,trim(label)
            error stop
        endif
    end function constructor
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Destructor
    !---------------------------------------------------------------------------
    subroutine destructor(this)
      type(element) :: this
      if(allocated(this%q0)) deallocate(this%q0)
      if(allocated(this%q1)) deallocate(this%q1)
      if(allocated(this%q2)) deallocate(this%q2)
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
        real(rp), dimension(:,:), allocatable  :: q0, q1, q2
        
        ! variables associated with the reading processes
        integer :: iostatus, funit

        namelist /element/ symbol, atomic_number, core, valence, f_core, num_quant_s, num_quant_p, num_quant_d,&
                           & q0, q1, q2
        ! Save previous values
        symbol = this%symbol
        atomic_number = this%atomic_number
        core = this%core
        valence = this%valence
        f_core = this%f_core
        num_quant_s = this%num_quant_s
        num_quant_p = this%num_quant_p
        num_quant_d = this%num_quant_d

        call move_alloc(this%q0,q0)
        call move_alloc(this%q1,q1)
        call move_alloc(this%q2,q2)

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
        
        call move_alloc(q0,this%q0)
        call move_alloc(q1,this%q1)
        call move_alloc(q2,this%q2)
    end subroutine build_from_file
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
        implicit none
        class(element), intent(out) :: this

        allocate(this%q0(3,2),this%q1(3,2),this%q2(3,2))

        this%q0(:,:) = 0.0d0
        this%q1(:,:) = 0.0d0
        this%q2(:,:) = 0.0d0
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
        real(rp), dimension(:,:), allocatable  :: q0, q1, q2

        namelist /element/ symbol, atomic_number, core, valence, f_core, num_quant_s, num_quant_p, num_quant_d,&
                           & q0, q1, q2

        symbol = this%symbol
        atomic_number = this%atomic_number
        core = this%core
        valence = this%valence
        f_core = this%f_core
        num_quant_s = this%num_quant_s
        num_quant_p = this%num_quant_p
        num_quant_d = this%num_quant_d

        if(allocated(this%q0)) then
            allocate(q0,mold=this%q0)
            q0 = this%q0
        else
            allocate(q0(0,0))
        endif
        if(allocated(this%q1)) then
            allocate(q1,mold=this%q1)
            q1 = this%q1
        else
            allocate(q1(0,0))
        endif
        if(allocated(this%q2)) then
            allocate(q2,mold=this%q0)
            q2 = this%q2
        else
            allocate(q2(0,0))
        endif

      
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
        real(rp), dimension(:,:), allocatable  :: q0, q1, q2

        namelist /element/ symbol, atomic_number, core, valence, f_core, num_quant_s, num_quant_p, num_quant_d,&
                           & q0, q1, q2

        symbol = this%symbol
        atomic_number = this%atomic_number
        core = this%core
        valence = this%valence
        f_core = this%f_core
        num_quant_s = this%num_quant_s
        num_quant_p = this%num_quant_p
        num_quant_d = this%num_quant_d

        if(allocated(this%q0)) then
            allocate(q0,mold=this%q0)
            q0 = this%q0
        else
            allocate(q0(0,0))
        endif
        if(allocated(this%q1)) then
            allocate(q1,mold=this%q1)
            q1 = this%q1
        else
            allocate(q1(0,0))
        endif
        if(allocated(this%q2)) then
            allocate(q2,mold=this%q2)
            q2 = this%q2
        else
            allocate(q2(0,0))
        endif

      
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
        integer :: i, j
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
