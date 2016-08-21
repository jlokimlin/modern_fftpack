module auxiliary_routines

    use, intrinsic :: ISO_Fortran_env, only: &
        stderr => ERROR_UNIT

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    implicit none

    private
    public :: get_real_1d_saved_workspace_length
    public :: get_real_1d_workspace_length
    public :: get_complex_1d_saved_workspace_length
    public :: get_complex_1d_workspace_length
    public :: get_cost_1d_saved_workspace_length
    public :: get_cost_1d_workspace_length
    public :: get_real_2d_saved_workspace_length
    public :: get_real_2d_workspace_length
    public :: get_complex_2d_saved_workspace_length
    public :: get_complex_2d_workspace_length
    public :: get_real_nd_saved_workspace_length
    public :: get_real_nd_workspace_length
    public :: get_complex_nd_saved_workspace_length
    public :: get_complex_nd_workspace_length
    public :: get_1d_saved_workspace_length
    public :: get_1d_saved_workspace
    public :: get_1d_workspace_length
    public :: get_1d_workspace
    public :: get_2d_saved_workspace_length
    public :: get_2d_saved_workspace
    public :: get_2d_workspace_length
    public :: get_2d_workspace
    public :: get_1d_sin_workspace_length
    public :: get_1d_sin_workspace
    public :: w2r, r2w

    ! error handlers
    public :: fft_consistent
    public :: fft_error_handler

    interface
        !
        !==> error_handlers
        !
        pure module function fft_consistent(inc, jump, n, lot) result (return_value)
            !
            !  Purpose:
            !
            !  Checks integer arguments inc, jump, n and lot for consistency.
            !  More specifically, inc, jump, n and lot are "consistent" if,
            !  for any values i1 and i2 < n, and j1 and j2 < lot,
            !
            !  i1 * inc + j1 * jump = i2 * inc + j2 * jump
            !
            !  can only occur if i1 = i2 and j1 = j2.
            !
            !  For multiple FFT's to execute correctly, inc, jump, n and lot must
            !  be consistent, or else at least one array element will be
            !  transformed more than once.
            !
            !  INPUTS
            !
            !  integer inc, jump, n, lot, the parameters to check.
            !
            !  OUTPUT
            !
            !  logical return_value is .true. if the parameters are consistent and
            ! .false. otherwise
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in) :: inc
            integer (ip), intent (in) :: jump
            integer (ip), intent (in) :: n
            integer (ip), intent (in) :: lot
            logical                   :: return_value
            !--------------------------------------------------------------
        end function fft_consistent

        module subroutine fft_error_handler(subroutine_name, info)
            !
            !  Purpose:
            !
            !  An error handler for fftpack routines.
            !
            !  It is called by an routine if an input parameter has an
            !  invalid value. A message is printed to standard output
            !  and execution stops.
            !
            !  Installers may consider modifying the stop statement in order to
            !  call system-specific exception-handling facilities.
            !
            !  Parameters:
            !
            !  INPUT
            !
            !  character (len=*) subroutine_name, the name of the calling routine.
            !
            !  integer info, an error code.  When a single invalid
            !  parameter in the parameter list of the calling routine has been detected,
            !  info is the position of that parameter.
            !
            !  In the case when an illegal
            !  combination of lot, jump, n, and inc has been detected, the calling
            !  subprogram calls this routine with info = -1.
            !
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip),      intent (in) :: info
            character (len=*), intent (in) :: subroutine_name
            !--------------------------------------------------------------
        end subroutine fft_error_handler

    end interface

    !---------------------------------------------------------------------------------
    ! Variables confined to the module
    !---------------------------------------------------------------------------------
    real (wp), parameter :: TWO = 2.0_wp
    !---------------------------------------------------------------------------------

contains

    pure function get_real_1d_saved_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lensav => return_value )

            lensav = n+int(log(real(n, kind=wp))/log(TWO), kind=ip)+4

        end associate

    end function get_real_1d_saved_workspace_length



    pure function get_real_1d_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = n

        end associate

    end function get_real_1d_workspace_length



    pure function get_complex_1d_saved_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lensav => return_value )

            lensav = 2*n+int(log(real(n, kind=wp))/log(TWO), kind=ip)+4

        end associate

    end function get_complex_1d_saved_workspace_length



    pure function get_complex_1d_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*n

        end associate

    end function get_complex_1d_workspace_length



    pure function get_cost_1d_saved_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n ! length of each sequence
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lensav => return_value )

            lensav = 2*n+int(log(real(n, kind=wp))/log(TWO))+4

        end associate

    end function get_cost_1d_saved_workspace_length



    pure function get_cost_1d_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = n-1

        end associate

    end function get_cost_1d_workspace_length



    pure function get_real_2d_saved_workspace_length(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lensav => return_value )

            lensav = l+3*m &
                +int(log(real(l, kind=wp))/log(TWO)) &
                +2*int(log(real(m, kind=wp))/log(TWO)) + 12

        end associate

    end function get_real_2d_saved_workspace_length



    pure function get_real_2d_workspace_length(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = m*(l+1)

        end associate

    end function get_real_2d_workspace_length

    pure function get_complex_2d_saved_workspace_length(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lensav => return_value )

            lensav = 2*(l+m) &
                + int(log(real(l, kind=wp))/log(TWO)) &
                +int(log(real(m, kind=wp))/log(TWO)) + 8

        end associate

    end function get_complex_2d_saved_workspace_length



    pure function get_complex_2d_workspace_length(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*l*m

        end associate

    end function get_complex_2d_workspace_length



    pure function get_real_nd_saved_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n ! length of each sequence
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lensav => return_value )

            lensav = n+int(log(real(n, kind=wp))/log(TWO))+4

        end associate

    end function get_real_nd_saved_workspace_length



    pure function get_real_nd_workspace_length(n, lot) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n   ! length of each sequence
        integer (ip), intent (in) :: lot ! number of sequences
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = n*lot

        end associate

    end function get_real_nd_workspace_length



    pure function get_complex_nd_saved_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n ! length of each sequence
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lensav => return_value )

            lensav = 2*n+int(log(real(n, kind=wp))/log(TWO))+4

        end associate

    end function get_complex_nd_saved_workspace_length



    pure function get_complex_nd_workspace_length(n, lot) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n   ! length of each sequence
        integer (ip), intent (in) :: lot ! number of sequences
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*n*lot

        end associate

    end function get_complex_nd_workspace_length

    pure function get_1d_saved_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------
        real (wp) :: temp
        !------------------------------------------------------------------

        associate( lensav => return_value )

            temp = log(real(n, kind=wp))/log(TWO)
            lensav = 2*n + int(temp, kind=ip) + 4

        end associate

    end function get_1d_saved_workspace_length


    pure function get_1d_saved_workspace(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip),           intent (in) :: n
        real (wp), allocatable              :: return_value(:)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip) :: lensav
        !------------------------------------------------------------------

        lensav = get_1d_saved_workspace_length(n)

        !
        !==> Allocate memory
        !
        allocate( return_value(lensav) )

    end function get_1d_saved_workspace



    pure function get_1d_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*n

        end associate

    end function get_1d_workspace_length



    pure function get_1d_workspace(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        real (wp), allocatable    :: return_value(:)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip) :: lenwrk
        !------------------------------------------------------------------

        lenwrk = get_1d_workspace_length(n)

        !
        !==> Allocate memory
        !
        allocate( return_value(lenwrk) )

    end function get_1d_workspace



    pure function get_2d_saved_workspace_length(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lensav => return_value )
            associate( &
                log_l => log(real(l,kind=wp)), &
                log_m => log(real(m,kind=wp)) &
                )

                lensav = 2*(l+m)+int(log_l, kind=ip)+int(log_m, kind=ip)+8

            end associate
        end associate

    end function get_2d_saved_workspace_length


    pure function get_2d_saved_workspace(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        real (wp), allocatable    :: return_value(:)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip) :: lensav
        !------------------------------------------------------------------

        lensav = get_2d_saved_workspace_length(l, m)

        !
        !==> Allocate memory
        !
        allocate( return_value(lensav) )

    end function get_2d_saved_workspace



    pure function get_2d_workspace_length(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*l*m

        end associate

    end function get_2d_workspace_length



    pure function get_2d_workspace(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        real (wp), allocatable    :: return_value(:)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip) :: lenwrk
        !------------------------------------------------------------------

        lenwrk = get_2d_workspace_length(l, m)

        !
        !==> Allocate memory
        !
        allocate( return_value(lenwrk) )

    end function get_2d_workspace



    pure function get_1d_sin_workspace_length(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*n + 2

        end associate

    end function get_1d_sin_workspace_length


    pure function get_1d_sin_workspace(n) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip),           intent (in) :: n
        real (wp), allocatable              :: return_value(:)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip) :: lenwrk
        !------------------------------------------------------------------

        lenwrk = get_1d_sin_workspace_length(n)

        !
        !==> Allocate memory
        !
        allocate( return_value(lenwrk) )

    end function get_1d_sin_workspace


    subroutine w2r(ldr, ldw, l, m, r, w)
        !
        ! Purpose:
        !
        ! Copies a 2D array, allowing for different leading dimensions.
        !
        !----------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------
        integer (ip), intent (in)     :: ldr
        integer (ip), intent (in)     :: ldw
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: l
        real (wp),    intent (in out) :: r(ldr,m)
        real (wp),    intent (in)     :: w(ldw,m)
                    !--------------------------------------------------------------

        r(1:l,:) = w(1:l,:)

    end subroutine w2r

    subroutine r2w(ldr, ldw, l, m, r, w)
        !
        ! Purpose:
        !
        ! copies a 2D array, allowing for different leading dimensions.
        !
        integer (ip), intent (in)     :: ldr
        integer (ip), intent (in)     :: ldw
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: l
        real (wp),    intent (in)     :: r(ldr,m)
        real (wp),    intent (in out) :: w(ldw,m)

        w(1:l,:) = r(1:l,:)

    end subroutine r2w

end module auxiliary_routines
