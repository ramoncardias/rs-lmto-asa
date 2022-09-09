!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Report
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
!> Module to process and report information
!------------------------------------------------------------------------------

module report_mod
    use string_mod, only: sl, startswith, fcount_str, fjoin, split, str_contains, int2str, indent_str, upper
    use precision_mod, only: rp
    implicit none

    private

    type report_data
        character(len=sl) :: label
        real(rp) :: value
    end type report_data
        
    type report
        character(len=sl) :: label
        real(rp), dimension(:), allocatable, private :: values
        type(report), dimension(:), allocatable :: children
    contains
        procedure, private :: report_add_value, report_add_label_value
        procedure, private :: report_add_value_r4, report_add_value_i
        procedure, private :: report_add_label_value_r4, report_add_label_value_i
        procedure, private :: push_child => report_push_child
        procedure, private :: report_get_values, report_get_values_with_label
        procedure, private :: get_child_values_internal => report_get_child_values_private
        procedure, private :: report_get_child_values, report_get_child_values_label
        procedure, private :: report_get_num_children, report_get_num_children_label
        procedure, private :: report_get_size_values, report_get_size_values_label
        procedure :: assign => report_assign
        procedure :: print_report => report_print_report
        procedure :: add_child => report_add_child
        procedure :: get_labels => report_labels
        procedure :: get_child => report_get_child
        procedure :: get_valid_child => report_get_valid_child
        procedure :: get_parent => report_get_parent_label
        generic :: add_value => report_add_value, report_add_label_value, report_add_value_r4, report_add_value_i, report_add_label_value_r4, report_add_label_value_i
        generic :: get_values => report_get_values, report_get_values_with_label
        generic :: get_child_values => report_get_child_values, report_get_child_values_label
        generic :: get_num_children => report_get_num_children, report_get_num_children_label
        generic :: get_size_values => report_get_size_values, report_get_size_values_label
    end type report

    interface report
        procedure :: report_constructor
    end interface report

    public :: report
    
    character(len=sl), parameter :: default_fields = 'MEAN,MAX,MIN,TOTAL,TOTAL (%),COUNT'
    character(len=*), parameter :: report_category_label = 'CATEGORY'
    integer, parameter :: report_field_margin = 3
    integer, parameter :: report_field_size = 7

contains
    recursive function report_constructor(label) result(obj)
        character(len=*), intent(in) :: label
        character(len=sl), dimension(:), allocatable :: lst_labels
        type(report) :: obj
        allocate(obj%values(0))
        if (str_contains(label,'.')) then
            call split(label,'.',lst_labels)
            obj%label = lst_labels(1)
            allocate(obj%children(1))
            obj%children(1) = report(fjoin(lst_labels(2:),'.'))
        else
            allocate(obj%children(0))
            obj%label = label
        endif
    end function report_constructor

    !> Assigns obj to this
    recursive subroutine report_assign(this,obj)
        class(report), intent(inout) :: this
        type(report), intent(in) :: obj
        integer :: i
        this%label = obj%label
        if (allocated(obj%values)) then
            if(allocated(this%values)) deallocate(this%values)
            allocate(this%values(size(obj%values)))
            this%values = obj%values
        endif
        if (allocated(obj%children)) then
            if(allocated(this%children)) deallocate(this%children)
            allocate(this%children(size(obj%children)))
            do i=1, size(obj%children)
                call this%children(i)%assign(obj%children(i)) 
            enddo
        endif
    end subroutine report_assign
    
    !> Append child argument variable to this report
    subroutine report_add_child(this,child)
        class(report), intent(inout) :: this
        type(report), intent(inout) :: child
        call this%push_child(child)
    end subroutine report_add_child
    
    !> Private subroutine to push a new child into children list as it is
    subroutine report_push_child(this,child)
        class(report), intent(inout) :: this
        type(report), intent(in) :: child
        type(report), dimension(:), allocatable :: aux_children
        integer :: i, n_childrens
        n_childrens = size(this%children) + 1
        allocate(aux_children(n_childrens))
        do i=1,n_childrens-1
            call aux_children(i)%assign(this%children(i))
        enddo
        call aux_children(n_childrens)%assign(child)
        deallocate(this%children)
        !TODO: try to use move_alloc
        allocate(this%children(n_childrens))
        do i=1,n_childrens
            call this%children(i)%assign(aux_children(i))
        enddo
    end subroutine report_push_child

    subroutine push_value(values,value)
        real(rp), dimension(:), allocatable, intent(inout) :: values
        real(rp), intent(in) :: value
        real(rp), dimension(:), allocatable :: aux_values
        integer :: new_size
        if (.not.allocated(values)) allocate(values(0))
        new_size = size(values)+1
        allocate(aux_values(new_size))
        aux_values(:new_size-1) = values
        aux_values(new_size) = value
        call move_alloc(aux_values,values)
    end subroutine push_value

    subroutine push_values(dst,src)
        real(rp), dimension(:), allocatable, intent(inout) :: dst
        real(rp), dimension(:), allocatable, intent(in) :: src
        real(rp), dimension(:), allocatable :: aux_values
        integer :: dst_size, src_size
        dst_size = size(dst)
        src_size = size(src)
        if (.not.allocated(dst)) then
            allocate(dst(src_size))
            dst = src
        else
            allocate(aux_values(dst_size + src_size))
            aux_values(:dst_size) = dst
            aux_values(dst_size+1:) = src
            call move_alloc(aux_values,dst)
        endif
    end subroutine push_values
    
    recursive subroutine report_add_label_value(this,label,value)
        class(report), intent(inout) :: this
        character(len=*), intent(in) :: label
        real(rp), intent(in) :: value
        type(report), pointer :: new_child, found_child
        type(report), target :: t_new_child
        character(len=sl):: new_label
        if (this%get_child(label, found_child)) then
            call found_child%add_value(value)
        else if (this%get_valid_child(label, found_child, new_label)) then
            call found_child%add_value(new_label,value)
        else
            new_label = label
            if (startswith(trim(label),trim(this%label)//'.')) new_label = label(len(trim(this%label))+2:)
            t_new_child = report(new_label)
            if (t_new_child%get_child(new_label,new_child)) call new_child%add_value(value)
            call this%add_child(t_new_child)
        endif
    end subroutine report_add_label_value

    subroutine report_add_value_r4(this, value)
        class(report), intent(inout) :: this
        real(4) :: value    
        call this%add_value(real(value,rp))
    end subroutine report_add_value_r4

    subroutine report_add_value_i(this, value)
        class(report), intent(inout) :: this
        integer :: value    
        call this%add_value(real(value,rp))
    end subroutine report_add_value_i

    subroutine report_add_label_value_r4(this,label,value)
        class(report), intent(inout) :: this
        character(len=*), intent(in) :: label
        real(4), intent(in) :: value
        call this%add_value(label,real(value,rp))
    end subroutine report_add_label_value_r4

    subroutine report_add_label_value_i(this,label,value)
        class(report), intent(inout) :: this
        character(len=*), intent(in) :: label
        integer, intent(in) :: value
        call this%add_value(label,real(value,rp))
    end subroutine report_add_label_value_i

    subroutine report_add_value(this, value)
        class(report) :: this
        real(rp), intent(in) :: value
        call push_value(this%values,value)
    end subroutine report_add_value

    recursive function report_get_child(this,label, child) result(found)
        class(report), intent(in) :: this
        character(len=*), intent(in) :: label
        type(report), intent(out), pointer :: child
        character(len=sl), dimension(:), allocatable :: lst_labels
        logical :: found
        integer :: i
        found = .False.
        if (this%label == label) then
            child = this
            found = .True.
        else if (startswith(label, this%label) .and. allocated(this%children)) then
            call split(label,'.',lst_labels)
            do i=1, size(this%children)
                if (lst_labels(2) == this%children(i)%label) then
                    found = this%children(i)%get_child(fjoin(lst_labels(2:),'.'), child)
                endif
            enddo
        endif
    end function report_get_child

    function report_get_valid_child(this, label, child, new_label) result(found)
        class(report), intent(in) :: this
        character(len=*), intent(in) :: label
        character(len=sl), intent(out), optional :: new_label
        type(report), intent(out), pointer :: child
        character(len=sl), dimension(:), allocatable :: lst_labels
        logical :: found
        integer :: i, init
        found = .False.
        call split(label,'.',lst_labels)
        init = merge(2,1,lst_labels(1) == this%label)
        do i=size(lst_labels),init,-1
            if (this%get_child(fjoin(lst_labels(:i),'.'),child)) then
                found = .True.
                if(present(new_label)) then
                  new_label = fjoin(lst_labels(i+1:),'.')
                endif
                exit
            endif
        enddo 
    end function report_get_valid_child

    function report_get_num_children(this) result(num)
        class(report), intent(in) :: this
        integer :: num
        num = merge(size(this%children),0,allocated(this%children))
    end function report_get_num_children

    function report_get_num_children_label(this,label) result(num)
        class(report), intent(in) :: this
        character(len=*), intent(in) :: label
        type(report), pointer :: child
        integer :: num
        num = 0
        if (this%get_child(label,child)) num = child%get_num_children()
    end function report_get_num_children_label

    function report_get_size_values(this) result(num)
        class(report), intent(in) :: this
        integer :: num
        num = merge(size(this%values),0,allocated(this%values))
    end function report_get_size_values

    function report_get_size_values_label(this,label) result(num)
        class(report), intent(in) :: this
        character(len=*), intent(in) :: label
        type(report), pointer :: child
        type(report), target :: t_child
        integer :: num
        num = 0
        if (this%get_child(label,child)) num = child%get_size_values()
    end function report_get_size_values_label

    recursive subroutine report_labels(this, labels, root)
        class(report), intent(in) :: this
        character(len=sl), dimension(:), allocatable, intent(inout) :: labels
        character(len=*), optional, intent(in) :: root
        character(len=sl) :: root_
        integer :: i, children_size
        if (.not.allocated(labels)) allocate(labels(0))
        root_=this%label
        if (present(root)) root_=trim(root)//'.'//trim(this%label)
        call push_labels(labels,root_)
        if (allocated(this%children)) then
            children_size=size(this%children)
            do i=1,children_size
                call this%children(i)%get_labels(labels,trim(root_))
            enddo
        endif
    end subroutine report_labels
    
    subroutine push_labels(labels,label)
        character(len=sl), dimension(:), allocatable, intent(inout) :: labels
        character(len=*), intent(in) :: label
        character(len=sl), dimension(:), allocatable :: aux_labels
        integer :: labels_size
        labels_size = size(labels)+1
        allocate(aux_labels(labels_size))
        aux_labels(:labels_size-1) = labels
        aux_labels(labels_size) = label
        deallocate(labels)
        call move_alloc(aux_labels, labels)
    end subroutine push_labels

    !> Search for a child with 'label', if found returns .True. and fills 'values' array with child's values, otherwise returns .False.
    !> You should enter child levels in 'label', ex: 'calc1.subcalc2.subsubcalc3' for getting the values of 'subsubcalc3' 
    subroutine report_get_values(this, values, append)
        class(report), intent(in) :: this
        real(rp), dimension(:), allocatable, intent(inout) :: values
        logical, optional, intent(in) :: append
        call this%get_values(this%label, values, append)
    end subroutine report_get_values

    subroutine report_get_values_with_label(this,label, values, append)
        class(report), intent(in) :: this
        character(len=*), intent(in) :: label
        real(rp), dimension(:), allocatable, intent(inout) :: values
        logical, optional, intent(in) :: append
        type(report), pointer :: child
        if (.not. (present(append) .and. append)) then
            if(allocated(values)) deallocate(values)
            allocate(values(0))
        endif
        if (this%get_child(label,child) .and. allocated(child%values)) then
            call push_values(values,child%values)
        endif
    end subroutine report_get_values_with_label

    subroutine report_get_child_values(this, values, append)
        class(report), intent(in) :: this
        real(rp), dimension(:), allocatable, intent(inout) :: values
        logical, optional, intent(in) :: append
        if (.not. (present(append) .and. append)) then
            if(allocated(values)) deallocate(values)
            allocate(values(0))
        endif
        call this%get_child_values_internal(values)
    end subroutine report_get_child_values

    subroutine report_get_child_values_label(this,label, values, append)
        class(report), intent(in) :: this
        character(len=*), intent(in) :: label
        real(rp), dimension(:), allocatable, intent(inout) :: values
        logical, optional, intent(in) :: append
        type(report), pointer :: child
        if (.not. (present(append) .and. append)) then
            if(allocated(values)) deallocate(values)
            allocate(values(0))
        endif
        if (this%get_child(label,child)) call child%get_child_values_internal(values)
    end subroutine report_get_child_values_label
    
    recursive subroutine report_get_child_values_private(this,values)
        class(report), intent(in) :: this
        real(rp), dimension(:), allocatable, intent(inout) :: values
        integer :: i
        if(.not. allocated(values)) allocate(values(0))
        if (allocated(this%children)) then
            do i=1,size(this%children)
                call push_values(values,this%children(i)%values)
                call this%children(i)%get_child_values_internal(values)
            enddo
        endif
    end subroutine report_get_child_values_private

    function report_get_parent_label(this,label,parent) result(found)
        class(report), intent(in) :: this
        character(len=*), intent(in) :: label
        type(report), pointer, intent(out) :: parent
        character(len=sl) :: parent_label
        character(len=sl), dimension(:), allocatable :: lst_labels
        logical :: found
        call split(label,'.',lst_labels)
        parent_label = fjoin(lst_labels(:size(lst_labels)-1),'.')
        found = this%get_child(parent_label,parent)
    end function report_get_parent_label

    subroutine get_short_labels(labels,short_labels,levels)
        character(len=sl), dimension(:), allocatable, intent(in) :: labels
        character(len=sl), dimension(:), allocatable, intent(out) :: short_labels
        integer, dimension(:), allocatable, intent(out) :: levels
        character(len=sl), dimension(:), allocatable :: lst_labels
        integer :: i
        allocate(short_labels(size(labels)))
        allocate(levels(size(labels)))
        do i=1,size(labels)
            call split(labels(i),'.',lst_labels)
            levels(i) = size(lst_labels)
            short_labels(i) = lst_labels(size(lst_labels))
        enddo
    end subroutine get_short_labels

    subroutine report_print_report(this,title,fields)
        class(report), intent(in) :: this
        character(len=*), intent(in) :: title
        character(len=*), intent(in), optional :: fields
        character(len=sl) :: fields_, aux_string_max, aux_string_min, aux_string_mean, aux_string_total, aux_string_total_perc, aux_string_count
        character(len=sl), dimension(:), allocatable :: lst_fields
        character(len=sl), dimension(:), allocatable :: labels, short_labels
        type(report), pointer :: parent
        integer, dimension(:), allocatable :: levels
        real(rp), dimension(:), allocatable :: values, parent_values
        integer :: i, j, k, short_label_max_size, label_max_size, fields_max_size, report_width, n_prints
        logical :: should_check_increasing, should_print_comment
        real(rp) :: sum_values, parent_sum_values
        ! fields format
        if (present(fields)) then
            fields_ = fields
        else
            fields_ = default_fields
        endif
        call split(fields_,',',lst_fields)
        fields_max_size = 0
        do i=1,size(lst_fields)
            fields_max_size = max(fields_max_size,len(trim(lst_fields(i))))
        enddo
        ! get labels
        call this%get_labels(labels)
        call get_short_labels(labels,short_labels,levels)
        short_label_max_size = 0
        do i=1, size(short_labels)
            call indent_str(short_labels(i),levels(i))
            short_label_max_size = max(short_label_max_size,len(trim(short_labels(i))))
        enddo
        ! increasing fields_max_size to prevent **** on report
        should_check_increasing = .true.
        j = 1
        do while (should_check_increasing)
            should_check_increasing = .false.
            do i=j,size(labels)
            call this%get_child_values(labels(i),values)
                ! checking just with TOTAL field
                call print_max(values,aux_string_max)
                call print_min(values,aux_string_min)
                call print_mean(values,aux_string_mean)
                call print_total(values,aux_string_total)
                call print_total_perc(labels(i),values,aux_string_total_perc)
                call print_count(values,aux_string_count)
                if ( &
                    str_contains(aux_string_max,'*') .or. &
                    str_contains(aux_string_min,'*') .or. &
                    str_contains(aux_string_mean,'*') .or. &
                    str_contains(aux_string_total,'*') .or. &
                    str_contains(aux_string_total_perc,'*') .or. &
                    str_contains(aux_string_count,'*') &
                ) then
                    fields_max_size = fields_max_size + 1
                    should_check_increasing = .true.
                    exit
                endif
                j = j + 1
            enddo
        enddo
        report_width = short_label_max_size+report_field_margin+size(lst_fields)*(report_field_margin+fields_max_size)
#ifdef DEBUG
        write(*,*) ('-',i=1,(report_width-len(title))/2-1),' ',upper(title),' ',('-',i=1,(report_width-len(title))/2-1)
        label_max_size = 0
        do i=1, size(labels)
            label_max_size = max(label_max_size,len(trim(labels(i))))
        enddo
        do i=1, size(labels)
            call this%get_values(labels(i),values)
            write(*,'(A'//int2str(label_max_size)//'," ",'//int2str(label_max_size)//'F'//int2str(fields_max_size)//'.4)') adjustl(labels(i)),values
        enddo
#endif
        ! print report title
        write(*,*) ('-',i=1,(report_width-len(title))/2-1),' ',upper(title),' ',('-',i=1,(report_width-len(title))/2-1)
        ! print report
        write(*,'(A," ",'//int2str(short_label_max_size-len(report_category_label)+report_field_margin)//'X)',advance='no') trim(report_category_label)
        !            ^  this white space is to align with '*' or ' ' bellow
        do i=1, size(lst_fields)
            lst_fields(i) = upper(lst_fields(i))
            write(*,'(A'//int2str(fields_max_size)//','//int2str(report_field_margin)//'X)',advance='no') trim(lst_fields(i))
        enddo
        write(*,*) ! new line
        should_print_comment = .False.
        do i=1, size(labels)
            n_prints = merge(2,1,(this%get_size_values(labels(i)) > 0) .and. (this%get_num_children(labels(i)) > 0))
            if(n_prints == 2) should_print_comment = .True.
            do k=1,n_prints
                call this%get_values(labels(i),values)
                if (k == 1) call this%get_child_values(labels(i),values,append=.True.)
                write(*,'("'//trim(short_labels(i))//merge(' ','*',k == 1)//'",'//int2str(short_label_max_size-len(trim(short_labels(i)))+report_field_margin)//'X)',advance='no') 
                ! print fields
                do j=1, size(lst_fields)
                    if (lst_fields(j) == 'MAX') then
                        call print_max(values)
                    else if (lst_fields(j) == 'MIN') then
                        call print_min(values)
                    else if (lst_fields(j) == 'MEAN') then
                        call print_mean(values)
                    else if (lst_fields(j) == 'TOTAL') then
                        call print_total(values)
                    else if (lst_fields(j) == 'TOTAL (%)') then
                        call print_total_perc(labels(i),values)
                    else if (lst_fields(j) == 'COUNT') then
                        call print_count(values)
                    endif
                enddo
                write(*,*) ! new line
            enddo
        enddo
        write(*,*) ('-',i=1,report_width)
        if (should_print_comment) then
            write(*,*) '* exclusive values to this label'
            write(*,*) ('-',i=1,report_width)
        endif
    contains
        subroutine print_max(values,string)
            real(rp), dimension(:), allocatable, intent(in) :: values
            character(len=*), optional, intent(out) :: string
            real(rp) :: value
            value = merge(maxval(values),real(0,8),size(values) > 0)
            if (present(string)) then
                write(string,'(F'//int2str(fields_max_size)//'.4,'//int2str(report_field_margin)//'X)') value
            else
                write(*,'(F'//int2str(fields_max_size)//'.4,'//int2str(report_field_margin)//'X)',advance='no') value
            endif
        end subroutine print_max
        subroutine print_min(values,string)
            real(rp), dimension(:), allocatable, intent(in) :: values
            character(len=*), optional, intent(out) :: string
            real(rp) :: value
            value = merge(minval(values),real(0,8),size(values) > 0)
            if (present(string)) then
                write(string,'(F'//int2str(fields_max_size)//'.4,'//int2str(report_field_margin)//'X)') value
            else
                write(*,'(F'//int2str(fields_max_size)//'.4,'//int2str(report_field_margin)//'X)',advance='no') value
            endif
        end subroutine print_min
        subroutine print_mean(values,string)
            real(rp), dimension(:), allocatable, intent(in) :: values
            character(len=*), optional, intent(out) :: string
            real(rp) :: value
            value = sum(values)/size(values)
            if (present(string)) then
                write(string,'(F'//int2str(fields_max_size)//'.4,'//int2str(report_field_margin)//'X)') value
            else
                write(*,'(F'//int2str(fields_max_size)//'.4,'//int2str(report_field_margin)//'X)',advance='no') value
            endif
        end subroutine print_mean
        subroutine print_total(values,string)
            real(rp), dimension(:), allocatable, intent(in) :: values
            character(len=*), optional, intent(out) :: string
            real(rp) :: value
            value = sum(values)
            if (present(string)) then
                write(string,'(F'//int2str(fields_max_size)//'.4,'//int2str(report_field_margin)//'X)') value
            else
                write(*,'(F'//int2str(fields_max_size)//'.4,'//int2str(report_field_margin)//'X)',advance='no') value
            endif
        end subroutine print_total
        subroutine print_total_perc(label,values,string)
            real(rp), dimension(:), allocatable, intent(in) :: values
            character(len=*), intent(in) :: label
            type(report), pointer :: parent
            character(len=*), optional, intent(out) :: string
            real(rp), dimension(:), allocatable :: parent_values
            real(rp) :: value
            if (this%get_parent(label,parent)) then
                call parent%get_child_values(parent_values)
                call parent%get_values(parent_values,append=.true.)
                parent_sum_values = sum(parent_values)
                value = 100*sum(values)/parent_sum_values
                if (present(string)) then
                    write(string,'(F'//int2str(fields_max_size)//'.1,'//int2str(report_field_margin)//'X)') value
                else
                    write(*,'(F'//int2str(fields_max_size)//'.1,'//int2str(report_field_margin)//'X)',advance='no') value
                endif
            else
                value = 100
                if (present(string)) then
                    write(string,'(F'//int2str(fields_max_size)//'.1,'//int2str(report_field_margin)//'X)') value
                else
                    write(*,'(F'//int2str(fields_max_size)//'.1,'//int2str(report_field_margin)//'X)',advance='no') value
                endif
            endif

        end subroutine print_total_perc
        subroutine print_count(values,string)
            real(rp), dimension(:), allocatable, intent(in) :: values
            character(len=*), optional, intent(out) :: string
            integer :: value
            value = size(values)
            if (present(string)) then
                write(string,'(I'//int2str(fields_max_size)//','//int2str(report_field_margin)//'X)') value
            else
                write(*,'(I'//int2str(fields_max_size)//','//int2str(report_field_margin)//'X)',advance='no') value
            endif
        end subroutine print_count
    end subroutine report_print_report

end module report_mod