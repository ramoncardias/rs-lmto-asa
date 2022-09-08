module globals_mod
    
    ! File/Folder character size
    integer, parameter :: GLOBAL_CHAR_SIZE = 300
    
    ! Project folder root
    character(len=*), parameter ::GLOBAL_ROOT_FOLDER = '/home/lucas/projetos/freelancer/rslmto_refac'
    
    ! Project database folder
    character(len=*), parameter :: GLOBAL_DATABASE_FOLDER = GLOBAL_ROOT_FOLDER // '/database/elements'

end module globals_mod