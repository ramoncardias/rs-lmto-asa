!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Safe_Alloc
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
!> Module to handle logging informations
!------------------------------------------------------------------------------


module safe_alloc_mod
    use string_mod, only: sl, int2str
    use logger_mod, only : g_logger
    use string_mod
    use report_mod
    implicit none

    ! private

    type alloc_infos
        character(len=sl) :: label
        character(len=9) :: type
        integer :: size
    contains
        procedure :: report => alloc_infos_report
    end type

    type :: safe_alloc
        type(alloc_infos), dimension(:), allocatable, private :: allocs_historic
    contains
        procedure, private :: push_new_info
        procedure, private :: alloc_int_1d
        procedure, private :: dealloc_int_1d
        procedure :: print_report => safe_alloc_report
        generic :: allocate => alloc_int_1d!, alloc_int_nd
        generic :: deallocate => dealloc_int_1d!, alloc_int_nd
        final :: destructor
    end type

    interface safe_alloc
        procedure :: constructor
    end interface safe_alloc

    interface free
        procedure :: destructor
    end interface free

    ! Global
    type(safe_alloc) :: g_safe_alloc

    ! Publics
    public :: g_safe_alloc, safe_alloc

contains

    function constructor() result(obj)
        type(safe_alloc) :: obj
        allocate(obj%allocs_historic(0))
    end function constructor

    subroutine destructor(this)
        type(safe_alloc) :: this
        character(len=sl) :: label
        integer :: i,j, allocated_size
        type(alloc_infos), dimension(:), allocatable :: new_allocs_historic
        do while (size(this%allocs_historic) > 1)
            label = this%allocs_historic(1)%label
            allocated_size = this%allocs_historic(1)%size
            ! list with new elements
            allocate(new_allocs_historic(size(this%allocs_historic)))
            j = 1
            do i = 2, size(this%allocs_historic)
                if (this%allocs_historic(i)%label == label) then
                    allocated_size = allocated_size + this%allocs_historic(i)%size
                else
                    new_allocs_historic(j) = this%allocs_historic(i)
                    j = j + 1
                endif
            enddo
            if (allocated_size > 0) then
                call g_logger%error('Possible memory leaky on list "' // trim(label) // '". Non-deallocated size: ' // int2str(allocated_size),__FILE__,__LINE__)
            else if (allocated_size < 0) then
                call g_logger%error('More deallocations than allocations registered for list "' // trim(label) // '". Allocation-deallication balance: ' // int2str(allocated_size),__FILE__,__LINE__)
            endif
            deallocate(this%allocs_historic)
            allocate(this%allocs_historic(j))
            this%allocs_historic = new_allocs_historic(:j)
        enddo
    end subroutine destructor

    subroutine push_new_info(this,new_alloc)
        class(safe_alloc), intent(inout) :: this
        type(alloc_infos), intent(in) :: new_alloc
        type(alloc_infos), dimension(:), allocatable :: aux_allocs_historic
        integer :: new_size
        new_size = size(this%allocs_historic) + 1
        allocate(aux_allocs_historic(new_size-1))
        aux_allocs_historic = this%allocs_historic
        deallocate(this%allocs_historic)
        allocate(this%allocs_historic(new_size))
        this%allocs_historic(:new_size-1) = aux_allocs_historic
        this%allocs_historic(new_size) = new_alloc
        deallocate(aux_allocs_historic)
        ! this%allocs_historic = [this%allocs_historic, new_alloc]
    end subroutine push_new_info
    
    subroutine alloc_int_1d(this,label,list,list_size)
        class(safe_alloc), intent(inout) :: this
        character(len=*), intent(in) :: label
        integer, intent(in) :: list_size
        integer, dimension(:), allocatable, intent(inout) :: list
        type(alloc_infos) :: new_alloc
        if (allocated(list)) then
            deallocate(list)
            call g_logger%warning('Reallocating list "' // trim(label) // '" of size ' // int2str(size(list)) // ' to new size ' // int2str(list_size),__FILE__,__LINE__)
        endif
        allocate(list(list_size))
        new_alloc%label = label
        new_alloc%type = 'integer'
        new_alloc%size = list_size
        call this%push_new_info(new_alloc)
    end subroutine alloc_int_1d

    subroutine dealloc_int_1d(this,label,list)
        class(safe_alloc), intent(inout) :: this
        character(len=*), intent(in) :: label
        integer, dimension(:), allocatable, intent(inout) :: list
        type(alloc_infos) :: new_alloc
        if (allocated(list)) then
            new_alloc%label = label
            new_alloc%type = 'integer'
            new_alloc%size = -size(list)
            call this%push_new_info(new_alloc)
            deallocate(list)
        else
            call g_logger%warning('Deallocating list "' // trim(label) // '" not allocated',__FILE__,__LINE__)
        endif
    end subroutine dealloc_int_1d

    subroutine safe_alloc_report(this)
        class(safe_alloc), intent(in) :: this
        type(report) :: rep
        integer :: i
        rep = report('total')
        do i=1, size(this%allocs_historic)
            if (this%allocs_historic(i)%size > 0) call rep%add_value(trim(rep%label)//'.'//this%allocs_historic(i)%label,this%allocs_historic(i)%size)
        enddo
        call rep%print_report('MEMORY REPORT')
    end subroutine safe_alloc_report

    subroutine alloc_infos_report(this)
        class(alloc_infos), intent(in) :: this
        character(len=11) :: action
        action = 'allocate'
        if (this%size < 0) action = 'deallocate'
        call g_logger%info('List ' // trim(this%label) // ' ' // trim(action) // ' ' // int2str(abs(this%size)) // ' of ' // trim(this%type), __FILE__, __LINE__)
    end subroutine alloc_infos_report

    ! subroutine alloc_int_nd(self,label,list,size)
    !     class(safe_alloc), intent(in) :: self
    !     character(len=*), intent(in) :: label
    !     integer, dimension(:), intent(in) :: size
    !     integer, dimension(:), allocatable, intent(inout) :: list
    !     type(alloc_infos) :: new_alloc
    !     allocate(list(size))
    !     new_alloc%label = label
    !     new_alloc%bytes = size * integer_bytes
    !     this%allocs = (this%allocs, new_alloc)
    ! end subroutine alloc_int_nd


end module safe_alloc_mod