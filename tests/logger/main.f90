program Main
    use globals_mod, only: g_logger
    implicit none

    integer :: a, b, c
    logical :: q, q_hidden
    character(len=254) :: name
    real(4) :: x
    real(8) :: y, y_hidden

    integer :: verbose = 1

    call g_logger%init(verbose)
    
    a = 1
    b = 2

    call g_logger%debug("First debug message",__FILE__,__LINE__)
    call g_logger%info( "First info message",__FILE__,__LINE__)
    call g_logger%warning("First warning message",__FILE__,__LINE__)
    call g_logger%error("First error message",__FILE__,__LINE__)
    ! uncomment bellow to test
    ! call g_logger%fatal("First fatal message",__FILE__,__LINE__)

    
    call g_logger%info( "Logging variables",__FILE__,__LINE__)
    
    call g_logger%log("a", a, __FILE__, __LINE__)
    call g_logger%log("b", b, __FILE__, __LINE__)

    call g_logger%info( "c = a + b",__FILE__,__LINE__)
    c = a + b
    call g_logger%log("c", c, __FILE__, __LINE__)

    call g_logger%info("a = c * b",__FILE__,__LINE__)
    a = c * b
    call g_logger%log("a", a, __FILE__, __LINE__)

    call g_logger%info("testing logical",__FILE__,__LINE__)
    q = .true.
    q_hidden = .true.
    call g_logger%log("q", q, __FILE__, __LINE__)
    call g_logger%log("q_hidden", q_hidden, __FILE__, __LINE__)

    call g_logger%info("testing character",__FILE__,__LINE__)
    name = "James"
    call g_logger%log("name", name, __FILE__, __LINE__)

    call g_logger%info("testing real(4)",__FILE__,__LINE__)
    x = 3.1415
    call g_logger%log("x", x, __FILE__, __LINE__)

    call g_logger%info("testing real(8)",__FILE__,__LINE__)
    y = 3.1415
    call g_logger%log("y", y, __FILE__, __LINE__)
    y_hidden = 3.1415
    call g_logger%log("y_hidden", y_hidden, __FILE__, __LINE__)

    call g_logger%info("End of Run",__FILE__,__LINE__)


end program Main