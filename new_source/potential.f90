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
    use string_mod, only: path_join, sl
    use logger_mod, only: g_logger
    use namelist_generator_mod, only: namelist_generator
    implicit none

    private

    ! public functions
    public :: array_of_potentials
    
    type, public :: potential
        !> Orbital index. Determines the size of the Hamiltonian
        integer :: lmax
        !> Potential parameters \f$ C \f$ and \f$ \sqrt{\Delta} \f$  
        !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
        !> 2nd index 1 = spin-up, 2 = spin-dw
        real(rp), dimension(:,:), allocatable :: center_band
        real(rp), dimension(:,:), allocatable :: width_band
        real(rp), dimension(:,:), allocatable :: gravity_center 
        real(rp), dimension(:,:), allocatable :: shifted_band ! center_band - gravity_center
        real(rp), dimension(:,:), allocatable :: obar 
        !> Potential parameters treated internally in the code
        !> cx -> center of the band
        !> wx -> width of the band
        !> cex -> center of the band - gravity center of the band
        !> obx -> o parameter for the 'hoh' calculation
        !> 1st index 1 = s-orbital, 2-4 = p-orbital, 5-9 = d-orbital
        !> 2nd index 1 = spin-up, 2 = spin-dw
        complex(rp), dimension(:,:), allocatable :: cx, wx, cex, obx

        !> cx0 -> cx-up + cx-down
        !> cx1 -> cx-up - cx-dwown
        !> wx0 -> wx-up + wx-down
        !> wx1 -> wx-up - wx-down
        complex(rp), dimension(:), allocatable :: cx0, cx1, wx0, wx1
        !> Potential parameters Pl's defined as \f$ P_l = 0.5 - \frac{1}{\pi}arctg(D_{l})\f$
        !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
        !> 2nd index 1 = spin-up, 2 = spin-dw
        real(rp), dimension(:,:), allocatable :: pl
        !> Moments as defined in Eq. 48. of Phys. Rev. B 43, 9538 (1991).
        !> 1st index 1 = s-orbital, 2 = p-orbital, 3 = d-orbital
        !> 2nd index 1 = spin-up, 2 = spin-dw
        real(rp), dimension(:,:,:), allocatable :: ql
        !> Potential parameters on the orthogonal basis
        real(rp), dimension(:,:), allocatable :: c, enu, ppar, qpar, srdel, vl, pnu
        !> Normalized magnetic moments
        real(rp), dimension(:), allocatable :: mom 
        !> Magnetic moments
        real(rp) :: mx, my, mz, mtot
        !> Band variables
        real(rp), dimension(:), allocatable :: cshi, dw_l
        !> Wignzer Seitz Radius
        real(rp) :: ws_r
        ! Energy variables
        real(rp) :: sumec, sumev, etot, utot, ekin, rhoeps
        ! Madelung potential
        real(rp) :: vmad
    contains
        procedure :: build_from_file
        procedure :: restore_to_default
        procedure :: print_state
        procedure :: print_state_full
        procedure :: print_state_formatted
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
        character(len=3*sl) :: path_to_file
        character(len=sl), dimension(2) :: lst_path_to_file
 
        call obj%restore_to_default()
        if(present(database)) then
            lst_path_to_file(1) = database
        else
            lst_path_to_file(1) = './'
        endif
        lst_path_to_file(2) = trim(label) // '.nml'
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
        !complex(rp), dimension(:,:), allocatable :: cex, obx
        real(rp), dimension(:,:), allocatable :: pl
        real(rp), dimension(:,:,:), allocatable :: ql
        real(rp), dimension(:,:), allocatable :: center_band
        real(rp), dimension(:,:), allocatable :: width_band
        real(rp), dimension(:,:), allocatable :: gravity_center
        real(rp), dimension(:), allocatable :: mom    
        real(rp), dimension(:,:), allocatable :: c, enu, ppar, qpar, srdel, vl
        real(rp) :: ws_r
        real(rp) :: sumec, sumev, etot, utot, ekin, rhoeps
        real(rp) :: vmad
        integer :: lmax
        ! variables associated with the reading processes
        integer :: iostatus, funit

        ! Local variables
        integer :: i

        namelist /par/ &
            center_band, width_band, gravity_center, &
            sumec, sumev, etot, utot, ekin, rhoeps, &
            c, enu, ppar, qpar, srdel, vl, &
            pl, mom, ws_r, ql, lmax, vmad

        ! Save previous values
        ws_r = this%ws_r
        sumec = this%sumec
        sumev = this%sumev
        etot = this%etot
        utot = this%utot
        ekin = this%ekin
        rhoeps = this%rhoeps
        lmax = this%lmax

        call move_alloc(this%ql,ql)
        call move_alloc(this%mom,mom)
        call move_alloc(this%pl,pl)

        call move_alloc(this%center_band,center_band)
        call move_alloc(this%width_band,width_band)
        call move_alloc(this%gravity_center,gravity_center)

        call move_alloc(this%c,c)
        call move_alloc(this%enu,enu)
        call move_alloc(this%ppar,ppar)
        call move_alloc(this%qpar,qpar)
        call move_alloc(this%srdel,srdel)
        call move_alloc(this%vl,vl)

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
        this%ws_r = ws_r
        this%sumec = sumec
        this%sumev = sumev
        this%etot = etot
        this%utot = utot
        this%ekin = ekin
        this%rhoeps = rhoeps
        this%lmax = lmax

        call move_alloc(center_band,this%center_band)
        call move_alloc(width_band,this%width_band)
        call move_alloc(gravity_center,this%gravity_center)

        call move_alloc(mom,this%mom)
        call move_alloc(pl,this%pl)
        call move_alloc(ql,this%ql)

        call move_alloc(c,this%c)
        call move_alloc(enu,this%enu)
        call move_alloc(ppar,this%ppar)
        call move_alloc(qpar,this%qpar)
        call move_alloc(srdel,this%srdel)
        call move_alloc(vl,this%vl)

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
    end subroutine build_from_file
    
    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    !> @brief
    !> Reset all members to default
    !---------------------------------------------------------------------------
    subroutine restore_to_default(this)
        implicit none
        class(potential), intent(out) :: this

        this%lmax = 2

        allocate(this%center_band(this%lmax+1,2),this%width_band(this%lmax+1,2),this%pl(0:this%lmax,2))
        allocate(this%gravity_center(this%lmax+1,2),this%ql(3,0:this%lmax,2),this%cx((this%lmax+1)**2,2),this%wx((this%lmax+1)**2,2))
        allocate(this%cex((this%lmax+1)**2,2),this%obx((this%lmax+1)**2,2),this%cx0((this%lmax+1)**2))
        allocate(this%cx1((this%lmax+1)**2),this%wx0((this%lmax+1)**2),this%wx1((this%lmax+1)**2)) 
        allocate(this%mom(3),this%cshi((2*(this%lmax+1))**2),this%dw_l((2*(this%lmax+1))**2))
        allocate(this%c(0:this%lmax,2),this%enu(0:this%lmax,2),this%ppar(0:this%lmax,2),this%qpar(0:this%lmax,2),&
                 this%srdel(0:this%lmax,2),this%vl(0:this%lmax,2),this%pnu(0:this%lmax,2))

        allocate(this%shifted_band(this%lmax+1,2),this%obar(this%lmax+1,2))

        this%ws_r = 0.0d0
        this%center_band(:,:) = 0.0d0
        this%width_band(:,:) = 0.0d0
        this%shifted_band(:,:) = 0.0d0
        this%obar(:,:) = 0.0d0
        this%sumec = 0.0d0
        this%sumev = 0.0d0
        this%etot = 0.0d0
        this%utot = 0.0d0
        this%ekin = 0.0d0
        this%rhoeps = 0.0d0
        this%vmad = 0.0d0
        this%pl(:,:) = 0.0d0
        this%ql(:,:,:) = 0.0d0
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
        this%c(:,:) = 0.0d0
        this%enu(:,:) = 0.0d0
        this%ppar(:,:) = 0.0d0
        this%qpar(:,:) = 0.0d0
        this%srdel(:,:) = 0.0d0
        this%vl(:,:) = 0.0d0
        this%pnu(:,:) = 0.0d0
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
        
        !complex(rp), dimension(:,:), allocatable :: cex, obx
        real(rp), dimension(:,:), allocatable :: pl
        real(rp), dimension(:,:,:), allocatable :: ql
        real(rp), dimension(:,:), allocatable :: center_band
        real(rp), dimension(:,:), allocatable :: width_band
        real(rp), dimension(:,:), allocatable :: gravity_center
        real(rp), dimension(:,:), allocatable :: c, enu, ppar, qpar, srdel, vl
        real(rp), dimension(3) :: mom    
        real(rp) :: ws_r
        real(rp) :: sumec, sumev, etot, utot, ekin, rhoeps
        real(rp) :: vmad
        integer :: lmax
        
        namelist /par/ &
            center_band, width_band, gravity_center, &
            sumec, sumev, etot, utot, ekin, rhoeps, vmad, &
            c, enu, ppar, qpar, srdel, vl, &
            pl, mom, ql, lmax, ws_r
        
        mom = this%mom        
        pl = this%pl
        ql = this%ql
        center_band = this%center_band
        width_band = this%width_band
        gravity_center = this%gravity_center
        sumec = this%sumec
        sumev = this%sumev
        etot = this%etot
        utot = this%utot
        ekin = this%ekin
        rhoeps = this%rhoeps
        vmad = this%vmad 
        lmax = this%lmax
        ws_r = this%ws_r
        c = this%c
        enu = this%enu
        ppar = this%ppar
        qpar = this%qpar
        srdel = this%srdel
        vl = this%vl
       
        if(present(unit) .and. present(file)) then
            call g_logger%fatal('Argument error: both unit and file are present',__FILE__,__LINE__)
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
        real(rp), dimension(:,:), allocatable :: pl
        real(rp), dimension(:,:,:), allocatable :: ql
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
            pl, mom, cshi, dw_l, cex, obx, ql
            
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
            call g_logger%fatal('Argument error: both unit and file are present',__FILE__,__LINE__)
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
    !> @brief
    !> Print class members values in namelist format 
    !>
    !> Print class members values in namelist format. Either unit or file should be provided. If none of them are provided, then the program will write to standart output.
    !> @param[in] unit File unit used to write namelist
    !> @param[in] file File name used to write namelist
    !---------------------------------------------------------------------------
    subroutine print_state_formatted(this,unit,file)
        implicit none
        class(potential), intent(in) :: this

        integer,intent(in),optional :: unit
        character(len=*),intent(in),optional :: file
        integer :: newunit
        
        type(namelist_generator) :: nml
            
        nml = namelist_generator('par')
        
        call nml%add('center_band', this%center_band)
        call nml%add('width_band', this%width_band)
        call nml%add('gravity_center', this%gravity_center)
        call nml%add('sumec', this%sumec)
        call nml%add('sumev', this%sumev)
        call nml%add('etot', this%etot)
        call nml%add('utot', this%utot)
        call nml%add('ekin', this%ekin)
        call nml%add('rhoeps', this%rhoeps)
        call nml%add('c', this%c)
        call nml%add('enu', this%enu)
        call nml%add('ppar', this%ppar)
        call nml%add('qpar', this%qpar)
        call nml%add('srdel', this%srdel)
        call nml%add('vl', this%vl)
        call nml%add('pl', this%pl)
        call nml%add('mom', this%mom)
        call nml%add('ws_r', this%ws_r)
        call nml%add('ql', this%ql)
        call nml%add('lmax', this%lmax)
        call nml%add('vmad', this%vmad)
        
        if(present(unit) .and. present(file)) then
            call g_logger%fatal('Argument error: both unit and file are present',__FILE__,__LINE__)
          else if(present(unit)) then
            call nml%generate_namelist(unit=unit)
          else if(present(file)) then
            call nml%generate_namelist(file=file)
          else
            call nml%generate_namelist()
          endif
    end subroutine print_state_formatted
    
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
