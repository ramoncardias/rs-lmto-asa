!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Potential
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
!> Module to hold potential's parameters
!------------------------------------------------------------------------------

module potential_mod
    use, intrinsic :: iso_fortran_env, only: error_unit
    use precision_mod, only: rp
    use globals_mod, only: GLOBAL_DATABASE_FOLDER, GLOBAL_CHAR_SIZE
    implicit none

    private

    ! public functions
    public :: array_of_potentials

    type, public :: potential
        !> Potential parameters \f$ C \f$ and \f$ \sqrt{\Delta} \f$  
        !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
        !> 2nd index 1 = spin-up, 2 = spin-dw
        real(rp), dimension(3,2) :: center_band
        real(rp), dimension(3,2) :: width_band
        real(rp), dimension(3,2) :: gravity_center 
        !> Potential parameters treated internally in the code
        !> cx -> center of the band
        !> wx -> width of the band
        !> cex -> center of the band - gravity center of the band
        !> obx -> o parameter for the 'hoh' calculation
        !> 1st index 1 = s-orbital, 2-4 = p-orbital, 5-9 = d-orbital
        !> 2nd index 1 = spin-up, 2 = spin-dw
        complex(rp), dimension(9,2) :: cx, wx, cex, obx

        !> cx0 -> cx-up + cx-down
        !> cx1 -> cx-up - cx-dwown
        !> wx0 -> wx-up + wx-down
        !> wx1 -> wx-up - wx-down
        complex(rp), dimension(9) :: cx0, cx1, wx0, wx1
        !> Potential parameters Pl's defined as \f$ P_l = 0.5 - \frac{1}{\pi}arctg(D_{l})\f$
        !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
        !> 2nd index 1 = spin-up, 2 = spin-dw
        real(rp), dimension(3,2) :: pl
        !> Normalized magnetic moments
        real(rp), dimension(3) :: mom 
        !> Magnetic moments
        real(rp) :: mx, my, mz, mtot
        !> Band variables
        real(rp), dimension(18) :: cshi, dw_l
    contains
        procedure :: build_from_file
        procedure :: restore_to_default
        procedure :: print_state
        procedure :: print_state_full
        final :: destructor
    end type potential

    interface potential
        procedure :: constructor
    end interface potential

contains

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Constructor
    !
    !> @param[in] potential Namelist file in database
    !> @param[in] database Directory to database files with 'potential' namelist
    !> @return type(control)
    !---------------------------------------------------------------------------
    function constructor(label,database) result(obj)
        type(potential) :: obj
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
            write(error_unit,'("[",A,":",I0,"]: Potential ",A," not found in any database")') __FILE__,__LINE__,trim(label)
            error stop
        endif
    end function constructor
    

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Destructor
    !---------------------------------------------------------------------------
    subroutine destructor(this)
      type(potential) :: this
    end subroutine destructor

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Read parameters from input file
    !
    !> @param[in] fname Namelist file
    !---------------------------------------------------------------------------
    subroutine build_from_file(this, fname)
        class(potential),intent(inout) :: this
        character(len=*), intent(in) :: fname

        ! Readable variables
        complex(rp), dimension(9,2) :: cex, obx
        real(rp), dimension(3,2) :: pl
        real(rp) :: &
            center_band_s_up, center_band_s_dw, &
            center_band_p_up, center_band_p_dw, &
            center_band_d_up, center_band_d_dw, &
            width_band_s_up, width_band_s_dw, &
            width_band_p_up, width_band_p_dw, &
            width_band_d_up, width_band_d_dw, &
            gravity_center_s_up, gravity_center_s_dw, &
            gravity_center_p_up, gravity_center_p_dw, &
            gravity_center_d_up, gravity_center_d_dw
        real(rp), dimension(3) :: mom    
        ! variables associated with the reading processes
        integer :: iostatus, funit

        ! Local variables
        integer :: i, j

        namelist /par/ &
            center_band_s_up, center_band_s_dw, &
            center_band_p_up, center_band_p_dw, &
            center_band_d_up, center_band_d_dw, &
            width_band_s_up, width_band_s_dw, &
            width_band_p_up, width_band_p_dw, &
            width_band_d_up, width_band_d_dw, &
            gravity_center_s_up, gravity_center_s_dw, &
            gravity_center_p_up, gravity_center_p_dw, &
            gravity_center_d_up, gravity_center_d_dw, &
            pl, mom, cex, obx

        ! Save previous values
        center_band_s_up = this%center_band(1,1)
        center_band_s_dw = this%center_band(1,2)
        center_band_p_up = this%center_band(2,1)
        center_band_p_dw = this%center_band(2,2)
        center_band_d_up = this%center_band(3,1)
        center_band_d_dw = this%center_band(3,2)

        width_band_s_up = this%width_band(1,1)
        width_band_s_dw = this%width_band(1,2)
        width_band_p_up = this%width_band(2,1)
        width_band_p_dw = this%width_band(2,2)
        width_band_d_up = this%width_band(3,1)
        width_band_d_dw = this%width_band(3,2)

        mom(:) = this%mom(:)
        
        pl = this%pl
        cex = this%cex
        obx = this%obx
        gravity_center_s_up = this%gravity_center(1,1)
        gravity_center_s_dw = this%gravity_center(1,2)
        gravity_center_p_up = this%gravity_center(2,1)
        gravity_center_p_dw = this%gravity_center(2,2)
        gravity_center_d_up = this%gravity_center(3,1)
        gravity_center_d_dw = this%gravity_center(3,2)

        open(newunit=funit,file=fname,action='read',iostat=iostatus,status='old')
        if(iostatus /= 0) then
            write(error_unit,'("[",A,":",I0,"]: file ",A," not found")') __FILE__,__LINE__,trim(fname)
            error stop
        endif
        
        read(funit,nml=par,iostat=iostatus)
        if( &
            iostatus /= 0 .and. &
            .not. IS_IOSTAT_END(iostatus) .and. &
            iostatus /= 5010 & ! LIBERROR_READ_VALUE (according to https://www.hep.manchester.ac.uk/u/samt/misc/gfortran_errors.html)
        ) then
            write(error_unit,'("[",A,":",I0,"]: Error while reading namelist")') __FILE__,__LINE__
            write(error_unit,'("iostatus = ",I0)') iostatus
            error stop
        endif
        close(funit)


        ! Setting user values
        this%center_band(1,1) = center_band_s_up
        this%center_band(1,2) = center_band_s_dw
        this%center_band(2,1) = center_band_p_up
        this%center_band(2,2) = center_band_p_dw
        this%center_band(3,1) = center_band_d_up
        this%center_band(3,2) = center_band_d_dw

        this%width_band(1,1) = width_band_s_up
        this%width_band(1,2) = width_band_s_dw
        this%width_band(2,1) = width_band_p_up
        this%width_band(2,2) = width_band_p_dw
        this%width_band(3,1) = width_band_d_up
        this%width_band(3,2) = width_band_d_dw

        this%mom(:) = mom(:)

        ! Setting the potential parameters 
        ! Imaginary part is set to 0
        ! kind is rp (same as real part)
        this%cx(1,:) = cmplx(this%center_band(1,:),0,rp)
        this%wx(1,:) = cmplx(this%width_band(1,:),0,rp)
        do i = 2,4
          this%cx(i,:) = cmplx(this%center_band(2,:),0,rp)
          this%wx(i,:) = cmplx(this%width_band(2,:),0,rp)
        end do
        do i = 5,9
          this%cx(i,:) = cmplx(this%center_band(3,:),0,rp)
          this%wx(i,:) = cmplx(this%width_band(3,:),0,rp)
        end do

        this%cx0(:) = 0.5d0*(this%cx(:,1)+this%cx(:,2))
        this%cx1(:) = 0.5d0*(this%cx(:,1)-this%cx(:,2))
        this%wx0(:) = 0.5d0*(this%wx(:,1)+this%wx(:,2))
        this%wx1(:) = 0.5d0*(this%wx(:,1)-this%wx(:,2))

        this%pl = pl
        this%cex = cex
        this%obx = obx

    end subroutine build_from_file
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
        implicit none
        class(potential), intent(out) :: this

        this%center_band(:,:) = 0.0d0
        this%width_band(:,:) = 0.0d0
        this%pl(:,:) = 0.0d0
        this%cx(:,:) = 0.0d0
        this%wx(:,:) = 0.0d0
        this%cex(:,:) = 0.0d0
        this%obx(:,:) = 0.0d0
        this%cx0(:) = 0.0d0
        this%cx1(:) = 0.0d0
        this%wx0(:) = 0.0d0
        this%wx1(:) = 0.0d0
        this%mom(:) = [0.0d0,0.0d0,1.0d0]
        this%cshi(:) = 0.0d0
        this%dw_l(:) = 1.0d0
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
        class(potential), intent(in) :: this

        integer,intent(in),optional :: unit
        character(len=*),intent(in),optional :: file
        integer :: newunit
        
        complex(rp), dimension(9,2) :: cex, obx
        real(rp), dimension(3,2) :: pl
        real(rp) :: &
            center_band_s_up, center_band_s_dw, &
            center_band_p_up, center_band_p_dw, &
            center_band_d_up, center_band_d_dw, &
            width_band_s_up, width_band_s_dw, &
            width_band_p_up, width_band_p_dw, &
            width_band_d_up, width_band_d_dw
        real(rp), dimension(3) :: mom    
        
        namelist /par/ &
            center_band_s_up, center_band_s_dw, &
            center_band_p_up, center_band_p_dw, &
            center_band_d_up, center_band_d_dw, &
            width_band_s_up, width_band_s_dw, &
            width_band_p_up, width_band_p_dw, &
            width_band_d_up, width_band_d_dw, &
            pl, mom, cex, obx
        
        center_band_s_up = this%center_band(1,1)
        center_band_s_dw = this%center_band(1,2)
        center_band_p_up = this%center_band(2,1)
        center_band_p_dw = this%center_band(2,2)
        center_band_d_up = this%center_band(3,1)
        center_band_d_dw = this%center_band(3,2)

        width_band_s_up = this%width_band(1,1)
        width_band_s_dw = this%width_band(1,2)
        width_band_p_up = this%width_band(2,1)
        width_band_p_dw = this%width_band(2,2)
        width_band_d_up = this%width_band(3,1)
        width_band_d_dw = this%width_band(3,2)

        mom = this%mom        

        pl = this%pl
        cex = this%cex
        obx = this%obx

        if(present(unit) .and. present(file)) then
            write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
            error stop
        else if(present(unit)) then
            write(unit,nml=par)
        else if(present(file)) then
            open(unit=newunit,file=file)
            write(newunit,nml=par)
            close(newunit)
        else
            write(*,nml=par)
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
        implicit none
        class(potential), intent(in) :: this

        integer,intent(in),optional :: unit
        character(len=*),intent(in),optional :: file
        integer :: newunit
        
        complex(rp), dimension(9,2) :: cx, wx, cex, obx
        complex(rp), dimension(9) :: cx0, cx1, wx0, wx1
        real(rp), dimension(3,2) :: pl
        real(rp), dimension(3) :: mom 
        real(rp), dimension(18) :: cshi, dw_l
        real(rp) :: &
            center_band_s_up, center_band_s_dw, &
            center_band_p_up, center_band_p_dw, &
            center_band_d_up, center_band_d_dw, &
            width_band_s_up, width_band_s_dw, &
            width_band_p_up, width_band_p_dw, &
            width_band_d_up, width_band_d_dw

        namelist /par/ &
            center_band_s_up, center_band_s_dw, &
            center_band_p_up, center_band_p_dw, &
            center_band_d_up, center_band_d_dw, &
            width_band_s_up, width_band_s_dw, &
            width_band_p_up, width_band_p_dw, &
            width_band_d_up, width_band_d_dw, &
            cx, wx, cx0, wx0, cx1, wx1,  &
            pl, mom, cshi, dw_l, cex, obx
            
        center_band_s_up = this%center_band(1,1)
        center_band_s_dw = this%center_band(1,2)
        center_band_p_up = this%center_band(2,1)
        center_band_p_dw = this%center_band(2,2)
        center_band_d_up = this%center_band(3,1)
        center_band_d_dw = this%center_band(3,2)

        width_band_s_up = this%width_band(1,1)
        width_band_s_dw = this%width_band(1,2)
        width_band_p_up = this%width_band(2,1)
        width_band_p_dw = this%width_band(2,2)
        width_band_d_up = this%width_band(3,1)
        width_band_d_dw = this%width_band(3,2)

        mom = this%mom
        
        cx = this%cx
        wx = this%wx
        cx0 = this%cx0
        wx0 = this%wx0
        cx1 = this%cx1
        wx1 = this%wx1

        pl = this%pl
        cshi = this%cshi
        dw_l = this%dw_l
        cex = this%cex
        obx = this%obx

        if(present(unit) .and. present(file)) then
            write(error_unit,'("[",A,":",I0,"]: Argument error: both unit and file are present")') __FILE__,__LINE__
            error stop
        else if(present(unit)) then
            write(unit,nml=par)
        else if(present(file)) then
            open(unit=newunit,file=file)
            write(newunit,nml=par)
            close(newunit)
        else
            write(*,nml=par)
        endif

    end subroutine print_state_full

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief Build an array of potentials
    !
    !> @param[in] potential List of labels in database to build potentials
    !> @param[in] database Directory to database files with 'potential' namelist
    !> @return type(potential), dimension(:), allocatable
    !---------------------------------------------------------------------------
    function array_of_potentials(potentials,database)
        type(potential), dimension(:), allocatable :: array_of_potentials
        character(len=*), dimension(:), intent(in) :: potentials
        character(len=*), intent(in), optional :: database
        integer :: i, j
        allocate(array_of_potentials(size(potentials)))
        if(present(database)) then
            do i=1, size(potentials)
                array_of_potentials(i) = potential(potentials(i),database)
            enddo
        else
            do i=1, size(potentials)
                array_of_potentials(i) = potential(potentials(i))
            enddo
        endif
    end function array_of_potentials

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
end module potential_mod
