module type_FFTpack

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stderr => ERROR_UNIT

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: FFTpack
    public :: cfft1i, cfft1b, cfft1f, cfft2i, cfft2b, cfft2f, cfftmi, cfftmb, cfftmf
    public :: rfft1i, rfft1b, rfft1f, rfft2i, rfft2b, rfft2f, rfftmi, rfftmb, rfftmf
    public :: cost1i, cost1b, cost1f, costmi, costmb, costmf
    public :: sint1i, sint1b, sint1f, sintmi, sintmb, sintmf
    public :: cosq1i, cosq1b, cosq1f, cosqmi, cosqmb, cosqmf
    public :: sinq1i, sinq1b, sinq1f, sinqmi, sinqmb, sinqmf

    ! Declare derived data type
    type, public :: FFTpack
        !----------------------------------------------------------------------
        ! Class variables
        !----------------------------------------------------------------------
        real (wp), allocatable :: saved_workspace(:)
        real (wp), allocatable :: workspace(:)
        !----------------------------------------------------------------------
    contains
        !----------------------------------------------------------------------
        ! Class methods
        !----------------------------------------------------------------------
        procedure, nopass, public :: get_1d_saved_workspace_size
        procedure, nopass, public :: get_1d_saved_workspace
        procedure, nopass, public :: get_1d_workspace_size
        procedure, nopass, public :: get_1d_workspace
        procedure, nopass, public :: get_1d_sin_workspace_size
        procedure, nopass, public :: get_1d_sin_workspace
        procedure, nopass, public :: get_2d_saved_workspace_size
        procedure, nopass, public :: get_2d_saved_workspace
        procedure, nopass, public :: get_2d_workspace_size
        procedure, nopass, public :: get_2d_workspace
        procedure,         public :: destroy => destroy_fftpack
        !----------------------------------------------------------------------
        ! Complex transform routines
        !----------------------------------------------------------------------
        procedure, nopass, public :: cfft1i ! 1-d complex initialization
        procedure, nopass, public :: cfft1b ! 1-d complex backward
        procedure, nopass, public :: cfft1f ! 1-d complex forward
        procedure, nopass, public :: cfft2i ! 2-d complex initialization
        procedure, nopass, public :: cfft2b ! 2-d complex backward
        procedure, nopass, public :: cfft2f ! 2-d complex forward
        procedure, nopass, public :: cfftmi ! multiple complex initialization
        procedure, nopass, public :: cfftmb ! multiple complex backward
        procedure, nopass, public :: cfftmf ! multiple complex forward
        !----------------------------------------------------------------------
        ! Real transform routines
        !----------------------------------------------------------------------
        procedure, nopass, public :: rfft1i ! 1-d real initialization
        procedure, nopass, public :: rfft1b ! 1-d real backward
        procedure, nopass, public :: rfft1f ! 1-d real forward
        procedure, nopass, public :: rfft2i ! 2-d real initialization
        procedure, nopass, public :: rfft2b ! 2-d real backward
        procedure, nopass, public :: rfft2f ! 2-d real forward
        procedure, nopass, public :: rfftmi ! multiple real initialization
        procedure, nopass, public :: rfftmb ! multiple real backward
        procedure, nopass, public :: rfftmf ! multiple real forward
        !----------------------------------------------------------------------
        ! Real cosine transform routines
        !----------------------------------------------------------------------
        procedure, nopass, public :: cost1i ! 1-d real cosine initialization
        procedure, nopass, public :: cost1b ! 1-d real cosine backward
        procedure, nopass, public :: cost1f ! 1-d real cosine forward
        procedure, nopass, public :: costmi ! multiple real cosine initialization
        procedure, nopass, public :: costmb ! multiple real cosine backward
        procedure, nopass, public :: costmf ! multiple real cosine forward
        !----------------------------------------------------------------------
        ! Real sine transform routines
        !----------------------------------------------------------------------
        procedure, nopass, public :: sint1i ! 1-d real sine initialization
        procedure, nopass, public :: sint1b ! 1-d real sine backward
        procedure, nopass, public :: sint1f ! 1-d real sine forward
        procedure, nopass, public :: sintmi ! multiple real sine initialization
        procedure, nopass, public :: sintmb ! multiple real sine backward
        procedure, nopass, public :: sintmf ! multiple real sine forward
        !----------------------------------------------------------------------
        ! Real quarter-cosine transform routines
        !----------------------------------------------------------------------
        procedure, nopass, public :: cosq1i ! 1-d real quarter-cosine initialization
        procedure, nopass, public :: cosq1b ! 1-d real quarter-cosine backward
        procedure, nopass, public :: cosq1f ! 1-d real quarter-cosine forward
        procedure, nopass, public :: cosqmi ! multiple real quarter-cosine initialization
        procedure, nopass, public :: cosqmb ! multiple real quarter-cosine backward
        procedure, nopass, public :: cosqmf ! multiple real quarter-cosine forward
        !----------------------------------------------------------------------
        ! Real quarter-sine transform routines
        !----------------------------------------------------------------------
        procedure, nopass, public :: sinq1i ! 1-d real quarter-sine initialization
        procedure, nopass, public :: sinq1b ! 1-d real quarter-sine backward
        procedure, nopass, public :: sinq1f ! 1-d real quarter-sine forward
        procedure, nopass, public :: sinqmi ! multiple real quarter-sine initialization
        procedure, nopass, public :: sinqmb ! multiple real quarter-sine backward
        procedure, nopass, public :: sinqmf ! multiple real quarter-sine forward
        !----------------------------------------------------------------------
    end type FFTpack

    ! Declare constructor
    interface FFTpack
        module procedure fftpack_2d_constructor
        module procedure fftpack_1d_constructor
    end interface

contains

    pure function fftpack_2d_constructor(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        type (FFTpack)            :: return_value
        !------------------------------------------------------------------

        !
        !==> Allocate memory
        !
        return_value %saved_workspace = get_2d_saved_workspace(l,m)
        return_value%workspace = get_2d_workspace(l,m)

    end function fftpack_2d_constructor


    pure function fftpack_1d_constructor(n) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        type (FFTpack)            :: return_value
        !------------------------------------------------------------------

        !
        !==> Allocate memory
        !
        return_value %saved_workspace = get_1d_saved_workspace(n)
        return_value%workspace = get_1d_workspace(n)

    end function fftpack_1d_constructor

    subroutine destroy_fftpack(this)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        class (FFTpack) , intent (in out) :: this
        !------------------------------------------------------------------

        !
        !==> Release memory
        !
        if (allocated(this%saved_workspace)) then
            deallocate( this%saved_workspace )
        end if

        if (allocated(this%workspace)) then
            deallocate( this%workspace )
        end if

    end subroutine destroy_fftpack


    pure function get_1d_saved_workspace_size(n) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------
        real (wp) :: temp
        !------------------------------------------------------------------

        associate( lensav => return_value )

            temp = log(real(n, kind=wp))/log(2.0_wp)
            lensav = 2*n + int(temp, kind=ip) + 4

        end associate

    end function get_1d_saved_workspace_size


    pure function get_1d_saved_workspace(n) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip),           intent (in) :: n
        real (wp), allocatable              :: return_value(:)
        !------------------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------------------
        integer (ip) :: lensav
        !------------------------------------------------------------------

        lensav = get_1d_saved_workspace_size(n)

        !
        !==> Allocate memory
        !
        allocate( return_value(lensav) )

    end function get_1d_saved_workspace



    pure function get_1d_workspace_size(n) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*n

        end associate

    end function get_1d_workspace_size



    pure function get_1d_workspace(n) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        real (wp), allocatable    :: return_value(:)
        !------------------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------------------
        integer (ip) :: lenwrk
        !------------------------------------------------------------------

        lenwrk = get_1d_workspace_size(n)

        !
        !==> Allocate memory
        !
        allocate( return_value(lenwrk) )

    end function get_1d_workspace



    pure function get_2d_saved_workspace_size(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
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

    end function get_2d_saved_workspace_size


    pure function get_2d_saved_workspace(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        real (wp), allocatable    :: return_value(:)
        !------------------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------------------
        integer (ip) :: lensav
        !------------------------------------------------------------------

        lensav = get_2d_saved_workspace_size(l, m)

        !
        !==> Allocate memory
        !
        allocate( return_value(lensav) )

    end function get_2d_saved_workspace



    pure function get_2d_workspace_size(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*l*m

        end associate

    end function get_2d_workspace_size



    pure function get_2d_workspace(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: l
        integer (ip), intent (in) :: m
        real (wp), allocatable    :: return_value(:)
        !------------------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------------------
        integer (ip) :: lenwrk
        !------------------------------------------------------------------

        lenwrk = get_2d_workspace_size(l, m)

        !
        !==> Allocate memory
        !
        allocate( return_value(lenwrk) )

    end function get_2d_workspace



    pure function get_1d_sin_workspace_size(n) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in) :: n
        integer (ip)              :: return_value
        !------------------------------------------------------------------

        associate( lenwrk => return_value )

            lenwrk = 2*n + 2

        end associate

    end function get_1d_sin_workspace_size



    pure function get_1d_sin_workspace(n) result (return_value)
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip),           intent (in) :: n
        real (wp), allocatable              :: return_value(:)
        !------------------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------------------
        integer (ip) :: lenwrk
        !------------------------------------------------------------------

        lenwrk = get_1d_sin_workspace_size(n)

        !
        !==> Allocate memory
        !
        allocate( return_value(lenwrk) )

    end function get_1d_sin_workspace



    subroutine cfft1i(n, wsave, lensav, ier)
        !
        ! cfft1i: initialization for cfft1b and cfft1f.
        !
        !  Purpose:
        !
        !  cfft1i initializes array wsave for use in its companion routines
        !  cfft1b and cfft1f. Routine cfft1i must be called before the first
        !  call to cfft1b or cfft1f, and after whenever the value of integer
        !  n changes.
        !
        !  Parameters:
        !
        !  input,
        !  n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n is a product
        !  of small primes.
        !
        !  input
        !  lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  output,
        !  wsave(lensav), containing the prime factors
        !  of n and  also containing certain trigonometric values which will be used
        !  in routines cfft1b or cfft1f.
        !
        !  output
        !  ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough.
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: wsave(lensav)
        integer (ip), intent (in)  :: lensav
        integer (ip), intent (out) :: ier
        !--------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if ( size(wsave) < get_1d_saved_workspace_size(n) ) then
            ier = 2
            call fft_error_handler('cfftmi ', 3)
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            associate( iw1 => 2*n+1 )

                call mcfti1(n,wsave,wsave(iw1),wsave(iw1+1))

            end associate
        end if

    end subroutine cfft1i


    subroutine cfft1b(n, inc, complex_data, lenc, wsave, lensav, work, lenwrk, ier)
        !
        !  INPUT
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array c, of two consecutive elements within the sequence to be transformed.
        !
        !  integer lenc, the dimension of the complex_data array.
        !  lenc must be at least inc*(n-1) + 1.
        !
        !  real wsave(lensav). wsave's contents must be initialized with a call
        !  to cfft1i before the first call to routine cfft1f
        !  or cfft1b for a given transform length n.  wsave's contents may be
        !  re-used for subsequent calls to cfft1f and cfft1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !
        !  input lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*n.
        !
        !  INPUT/OUTPUT
        !  complex complex_data(lenc) containing the sequence to be
        !  transformed.
        !
        !  real workspace work(lenwrk).
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenc not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  20, input error returned by lower level routine.
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: inc
        complex (wp), intent (in out) :: complex_data(lenc)
        integer (ip), intent (in)     :: lenc
        integer (ip), intent (in)     :: lensav
        integer (ip), intent (in)     :: lenwrk
        integer (ip), intent (out)    :: ier
        real (wp),    intent (in out) :: work(lenwrk)
        real (wp),    intent (in out) :: wsave(lensav)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip)           :: iw1
        real (wp), allocatable :: real_copy(:,:)
        !--------------------------------------------------------------

        !
        !==> Check validity of calling arguments
        !
        if (lenc < inc * (n - 1) + 1) then
            ier = 1
            call fft_error_handler('cfft1b ', 4)
        else if ( size(wsave) < get_1d_saved_workspace_size(n) ) then
            ier = 2
            call fft_error_handler('cfft1b ', 6)
        else if ( size(work) < 2 * n) then
            ier = 3
            call fft_error_handler('cfft1b ', 8)
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            !
            !==> Allocate memory
            !
            allocate( real_copy(2,size(complex_data)) )

            !
            !==> Make copy: complex to real
            !
            real_copy(1,:) = real(complex_data)
            real_copy(2,:) = aimag(complex_data)

            ! Set workspace index pointer
            iw1 = 2 * n + 1

            call c1fm1b(n, inc, real_copy, work, wsave, wsave(iw1), wsave(iw1+1:))

            !
            !==> Make copy: real to complex
            !
            complex_data =  cmplx(real_copy(1,:), real_copy(2,:), kind=wp)

            !
            !==> Release memory
            !
            deallocate( real_copy )
        end if

    contains

        subroutine c1fm1b(n, inc, c, ch, wa, fnf, fac)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: inc
            real (wp),    intent (in out) :: c(:, :)
            real (wp),    intent (out)    :: ch(:)
            real (wp),    intent (in out) :: wa(*)
            real (wp),    intent (in out) :: fnf
            real (wp),    intent (in out) :: fac(:)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip) :: ido, inc2, iip, iw
            integer (ip) :: k1, l1, l2, lid
            integer (ip) :: na, nbr, nf
            !----------------------------------------------------------------------

            inc2 = 2*inc
            nf = int(fnf, kind=ip)
            na = 0
            l1 = 1
            iw = 1

            do k1=1, nf
                iip = int(fac(k1), kind=ip)
                l2 = iip*l1
                ido = n/l2
                lid = l1*ido
                nbr = 1+na+2*min(iip-2, 4)
                select case (nbr)
                    case (1)
                        call c1f2kb(ido, l1, na, c, inc2, ch, 2, wa(iw))
                    case (2)
                        call c1f2kb(ido, l1, na, ch, 2, c, inc2, wa(iw))
                    case (3)
                        call c1f3kb(ido, l1, na, c, inc2, ch, 2, wa(iw))
                    case (4)
                        call c1f3kb(ido, l1, na, ch, 2, c, inc2, wa(iw))
                    case (5)
                        call c1f4kb(ido, l1, na, c, inc2, ch, 2, wa(iw))
                    case (6)
                        call c1f4kb(ido, l1, na, ch, 2, c, inc2, wa(iw))
                    case (7)
                        call c1f5kb(ido, l1, na, c, inc2, ch, 2, wa(iw))
                    case (8)
                        call c1f5kb(ido, l1, na, ch, 2, c, inc2, wa(iw))
                    case (9)
                        call c1fgkb(ido, iip, l1, lid, na, c, c, inc2, ch, ch, 2, wa(iw))
                    case (10)
                        call c1fgkb(ido, iip, l1, lid, na, ch, ch, 2, c, c, inc2, wa(iw))
                end select

                l1 = l2
                iw = iw+(iip-1)*(2*ido)

                if (iip <= 5) then
                    na = 1-na
                end if
            end do

        end subroutine c1fm1b

        subroutine c1f2kb(ido, l1, na, cc, in1, ch, in2, wa)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: na
            real (wp),    intent (in out) :: cc(in1,l1,ido,2)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,l1,2,ido)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido,1,2)
            !----------------------------------------------------------------------
            ! Dictionary: local variables
            !----------------------------------------------------------------------
            integer (ip)           :: i !! Counter
            real (wp), allocatable :: chold1(:), chold2(:)
            real (wp), allocatable :: ti2(:), tr2(:)
            !----------------------------------------------------------------------

            if (ido <= 1 .and. na /= 1) then
                !
                !==> Allocate memory
                !
                allocate( chold1(l1) )
                allocate( chold2(l1) )

                chold1 = cc(1,:,1,1)+cc(1,:,1,2)
                cc(1,:,1,2) = cc(1,:,1,1)-cc(1,:,1,2)
                cc(1,:,1,1) = chold1
                chold2 = cc(2,:,1,1)+cc(2,:,1,2)
                cc(2,:,1,2) = cc(2,:,1,1)-cc(2,:,1,2)
                cc(2,:,1,1) = chold2
                !
                !==> Release memory
                !
                deallocate( chold1 )
                deallocate( chold2 )
            else
                ch(1,:,1,1) = cc(1,:,1,1)+cc(1,:,1,2)
                ch(1,:,2,1) = cc(1,:,1,1)-cc(1,:,1,2)
                ch(2,:,1,1) = cc(2,:,1,1)+cc(2,:,1,2)
                ch(2,:,2,1) = cc(2,:,1,1)-cc(2,:,1,2)
                !
                !==> Allocate memory
                !
                allocate( tr2(l1) )
                allocate( ti2(l1) )

                do i=2,ido
                    ch(1,:,1,i) = cc(1,:,i,1)+cc(1,:,i,2)
                    tr2 = cc(1,:,i,1)-cc(1,:,i,2)
                    ch(2,:,1,i) = cc(2,:,i,1)+cc(2,:,i,2)
                    ti2 = cc(2,:,i,1)-cc(2,:,i,2)
                    ch(2,:,2,i) = wa(i,1,1)*ti2+wa(i,1,2)*tr2
                    ch(1,:,2,i) = wa(i,1,1)*tr2-wa(i,1,2)*ti2
                end do
                !
                !==> Release memory
                !
                deallocate( tr2 )
                deallocate( ti2 )
            end if

        end subroutine c1f2kb



        subroutine c1f3kb(ido, l1, na, cc, in1, ch, in2, wa)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: na
            real (wp),    intent (in out) :: cc(in1,l1,ido,3)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,l1,3,ido)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido,2,2)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip)           :: i !! Counter
            real (wp), allocatable :: di2(:), di3(:)
            real (wp), allocatable :: dr2(:), dr3(:)
            real (wp), allocatable :: ci2(:), ci3(:)
            real (wp), allocatable :: cr2(:), cr3(:)
            real (wp), allocatable :: ti2(:), tr2(:)
            real (wp), parameter   :: TAUI = sqrt(3.0_wp)/2!0.866025403784439_wp
            real (wp), parameter   :: TAUR = -0.5_wp
            !----------------------------------------------------------------------

            !
            !==> Allocate memory
            !
            allocate( ci2(l1), ci3(l1) )
            allocate( cr2(l1), cr3(l1) )
            allocate( ti2(l1), tr2(l1) )

            if (1 >= ido .and. na /= 1) then
                tr2 = cc(1,:,1,2)+cc(1,:,1,3)
                cr2 = cc(1,:,1,1)+TAUR*tr2
                cc(1,:,1,1) = cc(1,:,1,1)+tr2
                ti2 = cc(2,:,1,2)+cc(2,:,1,3)
                ci2 = cc(2,:,1,1)+TAUR*ti2
                cc(2,:,1,1) = cc(2,:,1,1)+ti2
                cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
                ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
                cc(1,:,1,2) = cr2-ci3
                cc(1,:,1,3) = cr2+ci3
                cc(2,:,1,2) = ci2+cr3
                cc(2,:,1,3) = ci2-cr3
            else
                tr2 = cc(1,:,1,2)+cc(1,:,1,3)
                cr2 = cc(1,:,1,1)+TAUR*tr2
                ch(1,:,1,1) = cc(1,:,1,1)+tr2
                ti2 = cc(2,:,1,2)+cc(2,:,1,3)
                ci2 = cc(2,:,1,1)+TAUR*ti2
                ch(2,:,1,1) = cc(2,:,1,1)+ti2
                cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
                ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
                ch(1,:,2,1) = cr2-ci3
                ch(1,:,3,1) = cr2+ci3
                ch(2,:,2,1) = ci2+cr3
                ch(2,:,3,1) = ci2-cr3

                !
                !==> Allocate memory
                !
                allocate( dr2(l1), dr3(l1) )
                allocate( di2(l1), di3(l1) )

                do i=2,ido
                    tr2 = cc(1,:,i,2)+cc(1,:,i,3)
                    cr2 = cc(1,:,i,1)+TAUR*tr2
                    ch(1,:,1,i) = cc(1,:,i,1)+tr2
                    ti2 = cc(2,:,i,2)+cc(2,:,i,3)
                    ci2 = cc(2,:,i,1)+TAUR*ti2
                    ch(2,:,1,i) = cc(2,:,i,1)+ti2
                    cr3 = TAUI*(cc(1,:,i,2)-cc(1,:,i,3))
                    ci3 = TAUI*(cc(2,:,i,2)-cc(2,:,i,3))

                    dr2 = cr2-ci3
                    dr3 = cr2+ci3
                    di2 = ci2+cr3
                    di3 = ci2-cr3

                    ch(2,:,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
                    ch(1,:,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
                    ch(2,:,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
                    ch(1,:,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3

                end do
                !
                !==> Release memory
                !
                deallocate( dr2, dr3 )
                deallocate( di2, di3 )
            end if

            !
            !==> Release memory
            !
            deallocate( ci2, ci3 )
            deallocate( cr2, cr3 )
            deallocate( ti2, tr2 )

        end subroutine c1f3kb


        subroutine c1f4kb(ido, l1, na, cc, in1, ch, in2, wa)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: na
            real (wp),    intent (in out) :: cc(in1,l1,ido,4)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,l1,4,ido)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido,3,2)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip)           :: i !! Counter
            real (wp), allocatable :: ci2(:), ci3(:), ci4(:)
            real (wp), allocatable :: cr2(:), cr3(:), cr4(:)
            real (wp), allocatable :: ti1(:), ti2(:), ti3(:), ti4(:)
            real (wp), allocatable :: tr1(:), tr2(:), tr3(:), tr4(:)
            !----------------------------------------------------------------------

            !
            !==> Allocate memory
            !
            allocate( ti1(l1), ti2(l1), ti3(l1), ti4(l1) )
            allocate( tr1(l1), tr2(l1), tr3(l1), tr4(l1) )

            if (1 >= ido .and. na /= 1) then
                ti1 = cc(2,:,1,1)-cc(2,:,1,3)
                ti2 = cc(2,:,1,1)+cc(2,:,1,3)
                tr4 = cc(2,:,1,4)-cc(2,:,1,2)
                ti3 = cc(2,:,1,2)+cc(2,:,1,4)
                tr1 = cc(1,:,1,1)-cc(1,:,1,3)
                tr2 = cc(1,:,1,1)+cc(1,:,1,3)
                ti4 = cc(1,:,1,2)-cc(1,:,1,4)
                tr3 = cc(1,:,1,2)+cc(1,:,1,4)
                cc(1,:,1,1) = tr2+tr3
                cc(1,:,1,3) = tr2-tr3
                cc(2,:,1,1) = ti2+ti3
                cc(2,:,1,3) = ti2-ti3
                cc(1,:,1,2) = tr1+tr4
                cc(1,:,1,4) = tr1-tr4
                cc(2,:,1,2) = ti1+ti4
                cc(2,:,1,4) = ti1-ti4
            else
                ti1 = cc(2,:,1,1)-cc(2,:,1,3)
                ti2 = cc(2,:,1,1)+cc(2,:,1,3)
                tr4 = cc(2,:,1,4)-cc(2,:,1,2)
                ti3 = cc(2,:,1,2)+cc(2,:,1,4)
                tr1 = cc(1,:,1,1)-cc(1,:,1,3)
                tr2 = cc(1,:,1,1)+cc(1,:,1,3)
                ti4 = cc(1,:,1,2)-cc(1,:,1,4)
                tr3 = cc(1,:,1,2)+cc(1,:,1,4)
                ch(1,:,1,1) = tr2+tr3
                ch(1,:,3,1) = tr2-tr3
                ch(2,:,1,1) = ti2+ti3
                ch(2,:,3,1) = ti2-ti3
                ch(1,:,2,1) = tr1+tr4
                ch(1,:,4,1) = tr1-tr4
                ch(2,:,2,1) = ti1+ti4
                ch(2,:,4,1) = ti1-ti4

                !
                !==> Allocate memory
                !
                allocate( ci2(l1), ci3(l1), ci4(l1) )
                allocate( cr2(l1), cr3(l1), cr4(l1) )

                do i=2,ido
                    ti1 = cc(2,:,i,1)-cc(2,:,i,3)
                    ti2 = cc(2,:,i,1)+cc(2,:,i,3)
                    ti3 = cc(2,:,i,2)+cc(2,:,i,4)
                    tr4 = cc(2,:,i,4)-cc(2,:,i,2)
                    tr1 = cc(1,:,i,1)-cc(1,:,i,3)
                    tr2 = cc(1,:,i,1)+cc(1,:,i,3)
                    ti4 = cc(1,:,i,2)-cc(1,:,i,4)
                    tr3 = cc(1,:,i,2)+cc(1,:,i,4)
                    ch(1,:,1,i) = tr2+tr3
                    cr3 = tr2-tr3
                    ch(2,:,1,i) = ti2+ti3
                    ci3 = ti2-ti3
                    cr2 = tr1+tr4
                    cr4 = tr1-tr4
                    ci2 = ti1+ti4
                    ci4 = ti1-ti4
                    ch(1,:,2,i) = wa(i,1,1)*cr2-wa(i,1,2)*ci2
                    ch(2,:,2,i) = wa(i,1,1)*ci2+wa(i,1,2)*cr2
                    ch(1,:,3,i) = wa(i,2,1)*cr3-wa(i,2,2)*ci3
                    ch(2,:,3,i) = wa(i,2,1)*ci3+wa(i,2,2)*cr3
                    ch(1,:,4,i) = wa(i,3,1)*cr4-wa(i,3,2)*ci4
                    ch(2,:,4,i) = wa(i,3,1)*ci4+wa(i,3,2)*cr4
                end do
                !
                !==> Release memory
                !
                deallocate( ci2, ci3, ci4 )
                deallocate( cr2, cr3, cr4 )
            end if

            !
            !==> Release memory
            !
            deallocate( ti1, ti2, ti3, ti4 )
            deallocate( tr1, tr2, tr3, tr4 )

        end subroutine c1f4kb

        subroutine c1f5kb(ido, l1, na, cc, in1, ch, in2, wa)
            !------------------------------------------------------------------
            ! Dictionary: calling arguments
            !------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: na
            real (wp),    intent (in out) :: cc(in1,l1,ido,5)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,l1,5,ido)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido,4,2)
            !------------------------------------------------------------------
            ! Dictionary: local variables
            !------------------------------------------------------------------
            integer (ip)           :: i !! Counter
            real (wp), allocatable :: chold1(:), chold2(:)
            real (wp), allocatable :: ci2(:), ci3(:), ci4(:), ci5(:)
            real (wp), allocatable :: cr2(:), cr3(:), cr4(:), cr5(:)
            real (wp), allocatable :: di2(:), di3(:), di4(:), di5(:)
            real (wp), allocatable :: dr2(:), dr3(:), dr4(:), dr5(:)
            real (wp), allocatable :: ti2(:), ti3(:), ti4(:), ti5(:)
            real (wp), allocatable :: tr2(:), tr3(:), tr4(:), tr5(:)
            real (wp), parameter   :: SQRT5 = sqrt(5.0_wp)
            real (wp), parameter   :: SQRT5_PLUS_5 = SQRT5 + 5.0_wp
            real (wp), parameter   :: TI11 = sqrt(SQRT5_PLUS_5/2)/2             ! 0.9510565162951536_wp
            real (wp), parameter   :: TI12 = sqrt(5.0_wp/(2.0_wp*SQRT5_PLUS_5)) ! 0.5877852522924731_wp
            real (wp), parameter   :: TR11 =  (SQRT5 - 1.0_wp)/4                 ! 0.3090169943749474_wp
            real (wp), parameter   :: TR12 = -(1.0_wp + SQRT5)/4                 !-0.8090169943749474_wp
            !------------------------------------------------------------------

            !
            !==> Allocate memory
            !
            allocate( chold1(l1), chold2(l1) )
            allocate( ti2(l1), ti3(l1), ti4(l1), ti5(l1) )
            allocate( tr2(l1),  tr3(l1), tr4(l1), tr5(l1) )
            allocate( cr2(l1), cr3(l1), cr4(l1), cr5(l1) )
            allocate( ci2(l1), ci3(l1), ci4(l1), ci5(l1) )

            if (1 >= ido .and. na /= 1) then
                ti5 = cc(2,:,1,2)-cc(2,:,1,5)
                ti2 = cc(2,:,1,2)+cc(2,:,1,5)
                ti4 = cc(2,:,1,3)-cc(2,:,1,4)
                ti3 = cc(2,:,1,3)+cc(2,:,1,4)
                tr5 = cc(1,:,1,2)-cc(1,:,1,5)
                tr2 = cc(1,:,1,2)+cc(1,:,1,5)
                tr4 = cc(1,:,1,3)-cc(1,:,1,4)
                tr3 = cc(1,:,1,3)+cc(1,:,1,4)
                chold1 = cc(1,:,1,1)+tr2+tr3
                chold2 = cc(2,:,1,1)+ti2+ti3
                cr2 = cc(1,:,1,1)+TR11*tr2+TR12*tr3
                ci2 = cc(2,:,1,1)+TR11*ti2+TR12*ti3
                cr3 = cc(1,:,1,1)+TR12*tr2+TR11*tr3
                ci3 = cc(2,:,1,1)+TR12*ti2+TR11*ti3
                cc(1,:,1,1) = chold1
                cc(2,:,1,1) = chold2
                cr5 = TI11*tr5+TI12*tr4
                ci5 = TI11*ti5+TI12*ti4
                cr4 = TI12*tr5-TI11*tr4
                ci4 = TI12*ti5-TI11*ti4
                cc(1,:,1,2) = cr2-ci5
                cc(1,:,1,5) = cr2+ci5
                cc(2,:,1,2) = ci2+cr5
                cc(2,:,1,3) = ci3+cr4
                cc(1,:,1,3) = cr3-ci4
                cc(1,:,1,4) = cr3+ci4
                cc(2,:,1,4) = ci3-cr4
                cc(2,:,1,5) = ci2-cr5
            else
                ti5 = cc(2,:,1,2)-cc(2,:,1,5)
                ti2 = cc(2,:,1,2)+cc(2,:,1,5)
                ti4 = cc(2,:,1,3)-cc(2,:,1,4)
                ti3 = cc(2,:,1,3)+cc(2,:,1,4)
                tr5 = cc(1,:,1,2)-cc(1,:,1,5)
                tr2 = cc(1,:,1,2)+cc(1,:,1,5)
                tr4 = cc(1,:,1,3)-cc(1,:,1,4)
                tr3 = cc(1,:,1,3)+cc(1,:,1,4)
                ch(1,:,1,1) = cc(1,:,1,1)+tr2+tr3
                ch(2,:,1,1) = cc(2,:,1,1)+ti2+ti3
                cr2 = cc(1,:,1,1)+TR11*tr2+TR12*tr3
                ci2 = cc(2,:,1,1)+TR11*ti2+TR12*ti3
                cr3 = cc(1,:,1,1)+TR12*tr2+TR11*tr3
                ci3 = cc(2,:,1,1)+TR12*ti2+TR11*ti3
                cr5 = TI11*tr5+TI12*tr4
                ci5 = TI11*ti5+TI12*ti4
                cr4 = TI12*tr5-TI11*tr4
                ci4 = TI12*ti5-TI11*ti4
                ch(1,:,2,1) = cr2-ci5
                ch(1,:,5,1) = cr2+ci5
                ch(2,:,2,1) = ci2+cr5
                ch(2,:,3,1) = ci3+cr4
                ch(1,:,3,1) = cr3-ci4
                ch(1,:,4,1) = cr3+ci4
                ch(2,:,4,1) = ci3-cr4
                ch(2,:,5,1) = ci2-cr5

                !
                !==> Allocate memory
                !
                allocate( di2(l1), di3(l1), di4(l1), di5(l1) )
                allocate( dr2(l1), dr3(l1), dr4(l1), dr5(l1) )

                do i=2,ido
                    ti5 = cc(2,:,i,2)-cc(2,:,i,5)
                    ti2 = cc(2,:,i,2)+cc(2,:,i,5)
                    ti4 = cc(2,:,i,3)-cc(2,:,i,4)
                    ti3 = cc(2,:,i,3)+cc(2,:,i,4)
                    tr5 = cc(1,:,i,2)-cc(1,:,i,5)
                    tr2 = cc(1,:,i,2)+cc(1,:,i,5)
                    tr4 = cc(1,:,i,3)-cc(1,:,i,4)
                    tr3 = cc(1,:,i,3)+cc(1,:,i,4)
                    ch(1,:,1,i) = cc(1,:,i,1)+tr2+tr3
                    ch(2,:,1,i) = cc(2,:,i,1)+ti2+ti3
                    cr2 = cc(1,:,i,1)+TR11*tr2+TR12*tr3
                    ci2 = cc(2,:,i,1)+TR11*ti2+TR12*ti3
                    cr3 = cc(1,:,i,1)+TR12*tr2+TR11*tr3
                    ci3 = cc(2,:,i,1)+TR12*ti2+TR11*ti3
                    cr5 = TI11*tr5+TI12*tr4
                    ci5 = TI11*ti5+TI12*ti4
                    cr4 = TI12*tr5-TI11*tr4
                    ci4 = TI12*ti5-TI11*ti4
                    dr3 = cr3-ci4
                    dr4 = cr3+ci4
                    di3 = ci3+cr4
                    di4 = ci3-cr4
                    dr5 = cr2+ci5
                    dr2 = cr2-ci5
                    di5 = ci2-cr5
                    di2 = ci2+cr5
                    ch(1,:,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
                    ch(2,:,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
                    ch(1,:,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
                    ch(2,:,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
                    ch(1,:,4,i) = wa(i,3,1)*dr4-wa(i,3,2)*di4
                    ch(2,:,4,i) = wa(i,3,1)*di4+wa(i,3,2)*dr4
                    ch(1,:,5,i) = wa(i,4,1)*dr5-wa(i,4,2)*di5
                    ch(2,:,5,i) = wa(i,4,1)*di5+wa(i,4,2)*dr5
                end do
                !
                !==> Release memory
                !
                deallocate( di2, di3, di4, di5 )
                deallocate( dr2, dr3, dr4, dr5 )
            end if

            !
            !==> Release memory
            !
            deallocate( ti2, ti3, ti4, ti5 )
            deallocate( tr2,  tr3, tr4, tr5 )
            deallocate( chold1, chold2 )
            deallocate( cr2, cr3, cr4, cr5 )
            deallocate( ci2, ci3, ci4, ci5 )

        end subroutine c1f5kb


        subroutine c1fgkb(ido, iip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa)
            !------------------------------------------------------------------
            ! Dictionary: calling arguments
            !------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: iip
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: lid
            integer (ip), intent (in)     :: na
            real (wp),    intent (in out) :: cc(in1,l1,iip,ido)
            real (wp),    intent (in out) :: cc1(in1,lid,iip)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,l1,ido,iip)
            real (wp),    intent (in out) :: ch1(in2,lid,iip)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido,iip-1,2)
            !------------------------------------------------------------------
            ! Dictionary: local variables
            !------------------------------------------------------------------
            integer (ip)           :: i, idlj, iipp2, iipph
            integer (ip)           :: j, jc, l, lc
            real (wp)              :: wai, war
            real (wp), allocatable :: chold1(:), chold2(:)
            !------------------------------------------------------------------

            iipp2 = iip+2
            iipph = (iip+1)/2

            ch1(1,:,1) = cc1(1,:,1)
            ch1(2,:,1) = cc1(2,:,1)

            do j=2,iipph
                jc = iipp2-j
                ch1(1,:,j) =  cc1(1,:,j)+cc1(1,:,jc)
                ch1(1,:,jc) = cc1(1,:,j)-cc1(1,:,jc)
                ch1(2,:,j) =  cc1(2,:,j)+cc1(2,:,jc)
                ch1(2,:,jc) = cc1(2,:,j)-cc1(2,:,jc)
            end do

            do j=2,iipph
                cc1(1,:,1) = cc1(1,:,1)+ch1(1,:,j)
                cc1(2,:,1) = cc1(2,:,1)+ch1(2,:,j)
            end do

            do l=2,iipph

                lc = iipp2-l
                cc1(1,:,l) = ch1(1,:,1)+wa(1,l-1,1)*ch1(1,:,2)
                cc1(1,:,lc) = wa(1,l-1,2)*ch1(1,:,iip)
                cc1(2,:,l) = ch1(2,:,1)+wa(1,l-1,1)*ch1(2,:,2)
                cc1(2,:,lc) = wa(1,l-1,2)*ch1(2,:,iip)

                do j=3,iipph
                    jc = iipp2-j
                    idlj = mod((l-1)*(j-1),iip)
                    war = wa(1,idlj,1)
                    wai = wa(1,idlj,2)
                    cc1(1,:,l) = cc1(1,:,l)+war*ch1(1,:,j)
                    cc1(1,:,lc) = cc1(1,:,lc)+wai*ch1(1,:,jc)
                    cc1(2,:,l) = cc1(2,:,l)+war*ch1(2,:,j)
                    cc1(2,:,lc) = cc1(2,:,lc)+wai*ch1(2,:,jc)
                end do
            end do

            if (1 >= ido .and. na /= 1) then
                !
                !==> Allocate memory
                !
                allocate( chold1(lid), chold2(lid) )

                do j=2,iipph
                    jc = iipp2-j
                    chold1 = cc1(1,:,j)-cc1(2,:,jc)
                    chold2 = cc1(1,:,j)+cc1(2,:,jc)
                    cc1(1,:,j) = chold1
                    cc1(2,:,jc) = cc1(2,:,j)-cc1(1,:,jc)
                    cc1(2,:,j) = cc1(2,:,j)+cc1(1,:,jc)
                    cc1(1,:,jc) = chold2
                end do
                !
                !==> Release memory
                !
                deallocate( chold1, chold2 )
            else
                ch1(1,:,1) = cc1(1,:,1)
                ch1(2,:,1) = cc1(2,:,1)

                do j=2,iipph
                    jc = iipp2-j
                    ch1(1,:,j) = cc1(1,:,j)-cc1(2,:,jc)
                    ch1(1,:,jc) = cc1(1,:,j)+cc1(2,:,jc)
                    ch1(2,:,jc) = cc1(2,:,j)-cc1(1,:,jc)
                    ch1(2,:,j) = cc1(2,:,j)+cc1(1,:,jc)
                end do

                if (ido /= 1) then
                    cc(1,:,1,:) = ch(1,:,:,1)
                    cc(2,:,1,:) = ch(2,:,:,1)

                    cc(1,:,2:iip,1) = ch(1,:,1,2:iip)
                    cc(2,:,2:iip,1) = ch(2,:,1,2:iip)

                    do j=2,iip
                        do i=2,ido
                            cc(1,:,j,i) = wa(i,j-1,1)*ch(1,:,i,j) &
                                -wa(i,j-1,2)*ch(2,:,i,j)
                            cc(2,:,j,i) = wa(i,j-1,1)*ch(2,:,i,j) &
                                +wa(i,j-1,2)*ch(1,:,i,j)
                        end do
                    end do
                end if
            end if

        end subroutine c1fgkb

    end subroutine cfft1b



    subroutine cfft1f(n, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
        !
        ! cfft1f: complex 64-bit precision forward fast fourier transform, 1d.
        !
        !  Purpose:
        !
        !  cfft1f computes the one-dimensional fourier transform of a single
        !  periodic sequence within a complex array. This transform is referred
        !  to as the forward transform or fourier analysis, transforming the
        !  sequence from physical to spectral space.
        !
        !  This transform is normalized since a call to cfft1f followed
        !  by a call to cfft1b (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  INPUT
        !
        !  integer n, the length of the sequence to be
        !  transformed. the transform is most efficient when
        !  n is a product of small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array c, of two consecutive elements within the sequence to be transformed.
        !
        !  real wsave(lensav).  wsave's contents must be
        !  initialized with a call to cfft1i before the first call to routine cfft1f
        !  or cfft1b for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to cfft1f and cfft1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*n.
        !
        !  INPUT/OUTPUT
        !  complex c(lenc) containing the sequence to be transformed.
        !
        !  real work(lenwrk), workspace array
        !  integer lenc, the dimension of the c array.
        !  lenc must be at least inc*(n-1) + 1.
        !
        !  OUTPUT
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenc not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  20, input error returned by lower level routine.
        !
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: inc
        complex (wp), intent (in out) :: c(lenc)
        integer (ip), intent (in)     :: lenc
        real (wp),    intent (in)     :: wsave(lensav)
        integer (ip), intent (in)     :: lensav
        real (wp),    intent (out)    :: work(lenwrk)
        integer (ip), intent (in)     :: lenwrk
        integer (ip), intent (out)    :: ier
        !------------------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------------------
        real (wp), allocatable :: real_copy(:,:)
        integer (ip)           :: iw1
        !------------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (lenc < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cfft1f ', 4)
        else if (lensav < get_1d_saved_workspace_size(n) ) then
            ier = 2
            call fft_error_handler('cfft1f ', 6)
        else if (lenwrk < 2*n) then
            ier = 3
            call fft_error_handler('cfft1f ', 8)
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then

            !
            !==> Allocate memory
            !
            allocate( real_copy(2,size(c)) )

            ! Make copy
            real_copy(1,:) = real(c)
            real_copy(2,:) = aimag(c)

            ! Set workspace index pointer
            iw1 = 2 * n+1

            call c1fm1f(n,inc,real_copy,work,wsave,wsave(iw1),wsave(iw1+1:))

            ! Make copy
            c =  cmplx(real_copy(1,:), real_copy(2,:), kind=wp)

            !
            !==> Release memory
            !
            deallocate( real_copy )
        end if

    contains

        subroutine c1fm1f(n, inc, c, ch, wa, fnf, fac)
            !----------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: inc
            real (wp),    intent (in out) :: c(:,:)
            real (wp),    intent (out)    :: ch(:)
            real (wp),    intent (in)     :: wa(*)
            real (wp),    intent (in)     :: fnf
            real (wp),    intent (in)     :: fac(:)
            !----------------------------------------------------------
            ! Dictionary: local variables
            !----------------------------------------------------------
            integer (ip) :: ido, inc2, iip, iw
            integer (ip) :: k1, l1, l2, lid
            integer (ip) :: na, nbr, nf
            !----------------------------------------------------------

            inc2 = 2*inc
            nf = int(fnf, kind=ip)
            na = 0
            l1 = 1
            iw = 1

            do k1=1,nf
                iip = int(fac(k1), kind=ip)
                l2 = iip*l1
                ido = n/l2
                lid = l1*ido
                nbr = 1+na+2*min(iip-2,4)
                select case (nbr)
                    case (1)
                        call c1f2kf(ido,l1,na,c,inc2,ch,2,wa(iw))
                    case (2)
                        call c1f2kf(ido,l1,na,ch,2,c,inc2,wa(iw))
                    case (3)
                        call c1f3kf(ido,l1,na,c,inc2,ch,2,wa(iw))
                    case (4)
                        call c1f3kf(ido,l1,na,ch,2,c,inc2,wa(iw))
                    case (5)
                        call c1f4kf(ido,l1,na,c,inc2,ch,2,wa(iw))
                    case (6)
                        call c1f4kf(ido,l1,na,ch,2,c,inc2,wa(iw))
                    case (7)
                        call c1f5kf(ido,l1,na,c,inc2,ch,2,wa(iw))
                    case (8)
                        call c1f5kf(ido,l1,na,ch,2,c,inc2,wa(iw))
                    case (9)
                        call c1fgkf(ido,iip,l1,lid,na,c,c,inc2,ch,ch,2,wa(iw))
                    case (10)
                        call c1fgkf(ido,iip,l1,lid,na,ch,ch,2,c,c,inc2,wa(iw))
                end select

                l1 = l2
                iw = iw+(iip-1)*(2*ido)
                if (iip <= 5) then
                    na = 1-na
                end if
            end do

        end subroutine c1fm1f


        subroutine c1f2kf(ido, l1, na, cc, in1, ch, in2, wa)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: na
            real (wp),    intent (in out) :: cc(in1,l1,ido,2)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,l1,2,ido)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido,1,2)
            !----------------------------------------------------------------------
            ! Dictionary: local variables
            !----------------------------------------------------------------------
            integer (ip)           :: i !! counter
            real (wp)              :: sn
            real (wp), allocatable :: temp1(:), temp2(:)
            real (wp), allocatable :: ti2(:), tr2(:)
            !----------------------------------------------------------------------

            if (1 >= ido) then
                sn = 1.0_wp/(2 * l1)
                if (na /= 1) then
                    !
                    !==> Allocate memory
                    !
                    allocate( temp1(l1) )
                    allocate( temp2(l1) )

                    temp1 = sn*(cc(1,:,1,1)+cc(1,:,1,2))
                    cc(1,:,1,2) = sn*(cc(1,:,1,1)-cc(1,:,1,2))
                    cc(1,:,1,1) = temp1
                    temp2 = sn*(cc(2,:,1,1)+cc(2,:,1,2))
                    cc(2,:,1,2) = sn*(cc(2,:,1,1)-cc(2,:,1,2))
                    cc(2,:,1,1) = temp2
                    !
                    !==> Release memory
                    !
                    deallocate( temp1 )
                    deallocate( temp2 )
                else
                    ch(1,:,1,1) = sn*(cc(1,:,1,1)+cc(1,:,1,2))
                    ch(1,:,2,1) = sn*(cc(1,:,1,1)-cc(1,:,1,2))
                    ch(2,:,1,1) = sn*(cc(2,:,1,1)+cc(2,:,1,2))
                    ch(2,:,2,1) = sn*(cc(2,:,1,1)-cc(2,:,1,2))
                end if
            else
                ch(1,:,1,1) = cc(1,:,1,1)+cc(1,:,1,2)
                ch(1,:,2,1) = cc(1,:,1,1)-cc(1,:,1,2)
                ch(2,:,1,1) = cc(2,:,1,1)+cc(2,:,1,2)
                ch(2,:,2,1) = cc(2,:,1,1)-cc(2,:,1,2)
                !
                !==> Allocate memory
                !
                allocate( tr2(l1) )
                allocate( ti2(l1) )

                do i=2,ido
                    ch(1,:,1,i) = cc(1,:,i,1)+cc(1,:,i,2)
                    tr2 = cc(1,:,i,1)-cc(1,:,i,2)
                    ch(2,:,1,i) = cc(2,:,i,1)+cc(2,:,i,2)
                    ti2 = cc(2,:,i,1)-cc(2,:,i,2)
                    ch(2,:,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
                    ch(1,:,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
                end do
                !
                !==> Release memory
                !
                deallocate( tr2 )
                deallocate( ti2 )
            end if

        end subroutine c1f2kf


        subroutine c1f3kf(ido, l1, na, cc, in1, ch, in2, wa)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: na
            real (wp),    intent (in out) :: cc(in1,l1,ido,3)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,l1,3,ido)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido,2,2)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip)           :: i !! Counter
            real (wp), allocatable :: ci2(:), ci3(:)
            real (wp), allocatable :: cr2(:), cr3(:)
            real (wp), allocatable :: di2(:), di3(:)
            real (wp), allocatable :: dr2(:), dr3(:)
            real (wp), allocatable :: ti2(:), tr2(:)
            real (wp), parameter   :: TAUI = -sqrt(3.0)/2!-0.866025403784439_wp
            real (wp), parameter   :: TAUR = -0.5_wp
            real (wp)              :: sn
            !----------------------------------------------------------------------

            !
            !==> Allocate memory
            !
            allocate( ci2(l1), ci3(l1) )
            allocate( cr2(l1), cr3(l1) )
            allocate( ti2(l1), tr2(l1) )

            if (1 >= ido) then
                sn = 1.0_wp/(3 * l1)
                if (na /= 1) then
                    tr2 = cc(1,:,1,2)+cc(1,:,1,3)
                    cr2 = cc(1,:,1,1)+TAUR*tr2
                    cc(1,:,1,1) = sn*(cc(1,:,1,1)+tr2)
                    ti2 = cc(2,:,1,2)+cc(2,:,1,3)
                    ci2 = cc(2,:,1,1)+TAUR*ti2
                    cc(2,:,1,1) = sn*(cc(2,:,1,1)+ti2)
                    cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
                    ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
                    cc(1,:,1,2) = sn*(cr2-ci3)
                    cc(1,:,1,3) = sn*(cr2+ci3)
                    cc(2,:,1,2) = sn*(ci2+cr3)
                    cc(2,:,1,3) = sn*(ci2-cr3)
                else
                    tr2 = cc(1,:,1,2)+cc(1,:,1,3)
                    cr2 = cc(1,:,1,1)+TAUR*tr2
                    ch(1,:,1,1) = sn*(cc(1,:,1,1)+tr2)
                    ti2 = cc(2,:,1,2)+cc(2,:,1,3)
                    ci2 = cc(2,:,1,1)+TAUR*ti2
                    ch(2,:,1,1) = sn*(cc(2,:,1,1)+ti2)
                    cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
                    ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
                    ch(1,:,2,1) = sn*(cr2-ci3)
                    ch(1,:,3,1) = sn*(cr2+ci3)
                    ch(2,:,2,1) = sn*(ci2+cr3)
                    ch(2,:,3,1) = sn*(ci2-cr3)
                end if
            else
                tr2 = cc(1,:,1,2)+cc(1,:,1,3)
                cr2 = cc(1,:,1,1)+TAUR*tr2
                ch(1,:,1,1) = cc(1,:,1,1)+tr2
                ti2 = cc(2,:,1,2)+cc(2,:,1,3)
                ci2 = cc(2,:,1,1)+TAUR*ti2
                ch(2,:,1,1) = cc(2,:,1,1)+ti2
                cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
                ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
                ch(1,:,2,1) = cr2-ci3
                ch(1,:,3,1) = cr2+ci3
                ch(2,:,2,1) = ci2+cr3
                ch(2,:,3,1) = ci2-cr3

                !
                !==> Allocate memory
                !
                allocate( di2(l1), di3(l1) )
                allocate( dr2(l1), dr3(l1) )

                do i=2,ido
                    tr2 = cc(1,:,i,2)+cc(1,:,i,3)
                    cr2 = cc(1,:,i,1)+TAUR*tr2
                    ch(1,:,1,i) = cc(1,:,i,1)+tr2
                    ti2 = cc(2,:,i,2)+cc(2,:,i,3)
                    ci2 = cc(2,:,i,1)+TAUR*ti2
                    ch(2,:,1,i) = cc(2,:,i,1)+ti2
                    cr3 = TAUI*(cc(1,:,i,2)-cc(1,:,i,3))
                    ci3 = TAUI*(cc(2,:,i,2)-cc(2,:,i,3))
                    dr2 = cr2-ci3
                    dr3 = cr2+ci3
                    di2 = ci2+cr3
                    di3 = ci2-cr3
                    ch(2,:,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
                    ch(1,:,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
                    ch(2,:,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
                    ch(1,:,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
                end do
                !
                !==> Release memory
                !
                deallocate( di2, di3 )
                deallocate( dr2, dr3 )
            end if

            !
            !==> Release memory
            !
            deallocate( ci2, ci3 )
            deallocate( cr2, cr3 )
            deallocate( ti2, tr2 )

        end subroutine c1f3kf


        subroutine c1f4kf(ido, l1, na, cc, in1, ch, in2, wa)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: na
            real (wp),    intent (in out) :: cc(in1,l1,ido,4)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,l1,4,ido)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido,3,2)
            !----------------------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------------------
            integer (ip) :: i
            real (wp)    :: sn
            !----------------------------------------------------------------------

            if (1 >= ido) then
                sn = 1.0_wp/(4 * l1)
                if (na /= 1) then
                    associate(&
                        ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
                        ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
                        tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
                        ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
                        tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
                        tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
                        ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
                        tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
                        )
                        cc(1,:,1,1) = sn*(tr2+tr3)
                        cc(1,:,1,3) = sn*(tr2-tr3)
                        cc(2,:,1,1) = sn*(ti2+ti3)
                        cc(2,:,1,3) = sn*(ti2-ti3)
                        cc(1,:,1,2) = sn*(tr1+tr4)
                        cc(1,:,1,4) = sn*(tr1-tr4)
                        cc(2,:,1,2) = sn*(ti1+ti4)
                        cc(2,:,1,4) = sn*(ti1-ti4)
                    end associate
                else
                    associate( &
                        ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
                        ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
                        tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
                        ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
                        tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
                        tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
                        ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
                        tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
                        )
                        ch(1,:,1,1) = sn*(tr2+tr3)
                        ch(1,:,3,1) = sn*(tr2-tr3)
                        ch(2,:,1,1) = sn*(ti2+ti3)
                        ch(2,:,3,1) = sn*(ti2-ti3)
                        ch(1,:,2,1) = sn*(tr1+tr4)
                        ch(1,:,4,1) = sn*(tr1-tr4)
                        ch(2,:,2,1) = sn*(ti1+ti4)
                        ch(2,:,4,1) = sn*(ti1-ti4)
                    end associate
                end if
            else
                associate( &
                    ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
                    ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
                    tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
                    ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
                    tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
                    tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
                    ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
                    tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
                    )
                    ch(1,:,1,1) = tr2+tr3
                    ch(1,:,3,1) = tr2-tr3
                    ch(2,:,1,1) = ti2+ti3
                    ch(2,:,3,1) = ti2-ti3
                    ch(1,:,2,1) = tr1+tr4
                    ch(1,:,4,1) = tr1-tr4
                    ch(2,:,2,1) = ti1+ti4
                    ch(2,:,4,1) = ti1-ti4
                end associate
                do i=2,ido
                    associate( &
                        ti1 => cc(2,:,i,1)-cc(2,:,i,3), &
                        ti2 => cc(2,:,i,1)+cc(2,:,i,3), &
                        ti3 => cc(2,:,i,2)+cc(2,:,i,4), &
                        tr4 => cc(2,:,i,2)-cc(2,:,i,4), &
                        tr1 => cc(1,:,i,1)-cc(1,:,i,3), &
                        tr2 => cc(1,:,i,1)+cc(1,:,i,3), &
                        ti4 => cc(1,:,i,4)-cc(1,:,i,2), &
                        tr3 => cc(1,:,i,2)+cc(1,:,i,4) &
                        )
                        ch(1,:,1,i) = tr2+tr3
                        associate( cr3 => tr2-tr3 )
                            ch(2,:,1,i) = ti2+ti3
                            associate( &
                                ci3 => ti2-ti3, &
                                cr2 => tr1+tr4, &
                                cr4 => tr1-tr4, &
                                ci2 => ti1+ti4, &
                                ci4 => ti1-ti4 &
                                )
                                ch(1,:,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
                                ch(2,:,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
                                ch(1,:,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
                                ch(2,:,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
                                ch(1,:,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
                                ch(2,:,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
                            end associate
                        end associate
                    end associate
                end do
            end if

        end subroutine c1f4kf


        subroutine c1f5kf(ido, l1, na, cc, in1, ch, in2, wa)
            !----------------------------------------------------------
            ! Dictionary: calling arguments
            !----------------------------------------------------------
            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1
            integer (ip) na
            real (wp) cc(in1,l1,ido,5)
            real (wp) ch(in2,l1,5,ido)
            real (wp) wa(ido,4,2)
            !----------------------------------------------------------
            ! Dictionary: local variables
            !----------------------------------------------------------
            integer (ip)         :: i, k
            real (wp)            :: sn
            real (wp)            :: chold1, chold2
            real (wp)            :: ci2, ci3, ci4, ci5
            real (wp)            :: cr2, cr3, cr4, cr5
            real (wp)            :: di2, di3, di4, di5
            real (wp)            :: dr2, dr3, dr4, dr5
            real (wp)            :: ti2, ti3, ti4, ti5
            real (wp)            :: tr2, tr3, tr4, tr5
            real (wp), parameter :: SQRT5 = sqrt(5.0_wp)
            real (wp), parameter :: SQRT5_PLUS_5 = SQRT5 + 5.0_wp
            real (wp), parameter :: TI11 = -sqrt(SQRT5_PLUS_5/2)/2             !-0.9510565162951536_wp
            real (wp), parameter :: TI12 = -sqrt(5.0_wp/(2.0_wp*SQRT5_PLUS_5)) !-0.5877852522924731_wp
            real (wp), parameter :: TR11 =  (SQRT5 - 1.0_wp)/4                 ! 0.3090169943749474_wp
            real (wp), parameter :: TR12 = -(1.0_wp + SQRT5)/4                 !-0.8090169943749474_wp
            !----------------------------------------------------------

            if (1 >= ido) then
                sn = 1.0_wp/(5 * l1)
                if (na /= 1) then
                    do k=1,l1
                        ti5 = cc(2,k,1,2)-cc(2,k,1,5)
                        ti2 = cc(2,k,1,2)+cc(2,k,1,5)
                        ti4 = cc(2,k,1,3)-cc(2,k,1,4)
                        ti3 = cc(2,k,1,3)+cc(2,k,1,4)
                        tr5 = cc(1,k,1,2)-cc(1,k,1,5)
                        tr2 = cc(1,k,1,2)+cc(1,k,1,5)
                        tr4 = cc(1,k,1,3)-cc(1,k,1,4)
                        tr3 = cc(1,k,1,3)+cc(1,k,1,4)
                        chold1 = sn*(cc(1,k,1,1)+tr2+tr3)
                        chold2 = sn*(cc(2,k,1,1)+ti2+ti3)
                        cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
                        ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
                        cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
                        ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
                        cc(1,k,1,1) = chold1
                        cc(2,k,1,1) = chold2
                        cr5 = TI11*tr5+TI12*tr4
                        ci5 = TI11*ti5+TI12*ti4
                        cr4 = TI12*tr5-TI11*tr4
                        ci4 = TI12*ti5-TI11*ti4
                        cc(1,k,1,2) = sn*(cr2-ci5)
                        cc(1,k,1,5) = sn*(cr2+ci5)
                        cc(2,k,1,2) = sn*(ci2+cr5)
                        cc(2,k,1,3) = sn*(ci3+cr4)
                        cc(1,k,1,3) = sn*(cr3-ci4)
                        cc(1,k,1,4) = sn*(cr3+ci4)
                        cc(2,k,1,4) = sn*(ci3-cr4)
                        cc(2,k,1,5) = sn*(ci2-cr5)
                    end do
                else
                    do k=1,l1
                        ti5 = cc(2,k,1,2)-cc(2,k,1,5)
                        ti2 = cc(2,k,1,2)+cc(2,k,1,5)
                        ti4 = cc(2,k,1,3)-cc(2,k,1,4)
                        ti3 = cc(2,k,1,3)+cc(2,k,1,4)
                        tr5 = cc(1,k,1,2)-cc(1,k,1,5)
                        tr2 = cc(1,k,1,2)+cc(1,k,1,5)
                        tr4 = cc(1,k,1,3)-cc(1,k,1,4)
                        tr3 = cc(1,k,1,3)+cc(1,k,1,4)
                        ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2+tr3)
                        ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2+ti3)
                        cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
                        ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
                        cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
                        ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
                        cr5 = TI11*tr5+TI12*tr4
                        ci5 = TI11*ti5+TI12*ti4
                        cr4 = TI12*tr5-TI11*tr4
                        ci4 = TI12*ti5-TI11*ti4
                        ch(1,k,2,1) = sn*(cr2-ci5)
                        ch(1,k,5,1) = sn*(cr2+ci5)
                        ch(2,k,2,1) = sn*(ci2+cr5)
                        ch(2,k,3,1) = sn*(ci3+cr4)
                        ch(1,k,3,1) = sn*(cr3-ci4)
                        ch(1,k,4,1) = sn*(cr3+ci4)
                        ch(2,k,4,1) = sn*(ci3-cr4)
                        ch(2,k,5,1) = sn*(ci2-cr5)
                    end do
                end if
            else
                do k=1,l1
                    ti5 = cc(2,k,1,2)-cc(2,k,1,5)
                    ti2 = cc(2,k,1,2)+cc(2,k,1,5)
                    ti4 = cc(2,k,1,3)-cc(2,k,1,4)
                    ti3 = cc(2,k,1,3)+cc(2,k,1,4)
                    tr5 = cc(1,k,1,2)-cc(1,k,1,5)
                    tr2 = cc(1,k,1,2)+cc(1,k,1,5)
                    tr4 = cc(1,k,1,3)-cc(1,k,1,4)
                    tr3 = cc(1,k,1,3)+cc(1,k,1,4)
                    ch(1,k,1,1) = cc(1,k,1,1)+tr2+tr3
                    ch(2,k,1,1) = cc(2,k,1,1)+ti2+ti3
                    cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
                    ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
                    cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
                    ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
                    cr5 = TI11*tr5+TI12*tr4
                    ci5 = TI11*ti5+TI12*ti4
                    cr4 = TI12*tr5-TI11*tr4
                    ci4 = TI12*ti5-TI11*ti4
                    ch(1,k,2,1) = cr2-ci5
                    ch(1,k,5,1) = cr2+ci5
                    ch(2,k,2,1) = ci2+cr5
                    ch(2,k,3,1) = ci3+cr4
                    ch(1,k,3,1) = cr3-ci4
                    ch(1,k,4,1) = cr3+ci4
                    ch(2,k,4,1) = ci3-cr4
                    ch(2,k,5,1) = ci2-cr5
                end do
                do i=2,ido
                    do k=1,l1
                        ti5 = cc(2,k,i,2)-cc(2,k,i,5)
                        ti2 = cc(2,k,i,2)+cc(2,k,i,5)
                        ti4 = cc(2,k,i,3)-cc(2,k,i,4)
                        ti3 = cc(2,k,i,3)+cc(2,k,i,4)
                        tr5 = cc(1,k,i,2)-cc(1,k,i,5)
                        tr2 = cc(1,k,i,2)+cc(1,k,i,5)
                        tr4 = cc(1,k,i,3)-cc(1,k,i,4)
                        tr3 = cc(1,k,i,3)+cc(1,k,i,4)
                        ch(1,k,1,i) = cc(1,k,i,1)+tr2+tr3
                        ch(2,k,1,i) = cc(2,k,i,1)+ti2+ti3
                        cr2 = cc(1,k,i,1)+TR11*tr2+TR12*tr3
                        ci2 = cc(2,k,i,1)+TR11*ti2+TR12*ti3
                        cr3 = cc(1,k,i,1)+TR12*tr2+TR11*tr3
                        ci3 = cc(2,k,i,1)+TR12*ti2+TR11*ti3
                        cr5 = TI11*tr5+TI12*tr4
                        ci5 = TI11*ti5+TI12*ti4
                        cr4 = TI12*tr5-TI11*tr4
                        ci4 = TI12*ti5-TI11*ti4
                        dr3 = cr3-ci4
                        dr4 = cr3+ci4
                        di3 = ci3+cr4
                        di4 = ci3-cr4
                        dr5 = cr2+ci5
                        dr2 = cr2-ci5
                        di5 = ci2-cr5
                        di2 = ci2+cr5
                        ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
                        ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
                        ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
                        ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
                        ch(1,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
                        ch(2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
                        ch(1,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
                        ch(2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
                    end do
                end do
            end if

        end subroutine c1f5kf


        subroutine c1fgkf(ido, iip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) iip
            integer (ip) l1
            integer (ip) lid

            real (wp) cc(in1,l1,iip,ido)
            real (wp) cc1(in1,lid,iip)
            real (wp) ch(in2,l1,ido,iip)
            real (wp) ch1(in2,lid,iip)
            real (wp) chold1
            real (wp) chold2
            integer (ip) i
            integer (ip) idlj
            integer (ip) iipp2
            integer (ip) iipph
            integer (ip) j
            integer (ip) jc
            integer (ip) ki
            integer (ip) l
            integer (ip) lc
            integer (ip) na
            real (wp) sn
            real (wp) wa(ido,iip-1,2)
            real (wp) wai
            real (wp) war

            iipp2 = iip+2
            iipph = (iip+1)/2
            do ki=1,lid
                ch1(1,ki,1) = cc1(1,ki,1)
                ch1(2,ki,1) = cc1(2,ki,1)
            end do

            do j=2,iipph
                jc = iipp2-j
                do ki=1,lid
                    ch1(1,ki,j) =  cc1(1,ki,j)+cc1(1,ki,jc)
                    ch1(1,ki,jc) = cc1(1,ki,j)-cc1(1,ki,jc)
                    ch1(2,ki,j) =  cc1(2,ki,j)+cc1(2,ki,jc)
                    ch1(2,ki,jc) = cc1(2,ki,j)-cc1(2,ki,jc)
                end do
            end do

            do j=2,iipph
                do ki=1,lid
                    cc1(1,ki,1) = cc1(1,ki,1)+ch1(1,ki,j)
                    cc1(2,ki,1) = cc1(2,ki,1)+ch1(2,ki,j)
                end do
            end do

            do l=2,iipph
                lc = iipp2-l
                do ki=1,lid
                    cc1(1,ki,l) = ch1(1,ki,1)+wa(1,l-1,1)*ch1(1,ki,2)
                    cc1(1,ki,lc) = -wa(1,l-1,2)*ch1(1,ki,iip)
                    cc1(2,ki,l) = ch1(2,ki,1)+wa(1,l-1,1)*ch1(2,ki,2)
                    cc1(2,ki,lc) = -wa(1,l-1,2)*ch1(2,ki,iip)
                end do
                do j=3,iipph
                    jc = iipp2-j
                    idlj = mod((l-1)*(j-1),iip)
                    war = wa(1,idlj,1)
                    wai = -wa(1,idlj,2)
                    do ki=1,lid
                        cc1(1,ki,l) = cc1(1,ki,l)+war*ch1(1,ki,j)
                        cc1(1,ki,lc) = cc1(1,ki,lc)+wai*ch1(1,ki,jc)
                        cc1(2,ki,l) = cc1(2,ki,l)+war*ch1(2,ki,j)
                        cc1(2,ki,lc) = cc1(2,ki,lc)+wai*ch1(2,ki,jc)
                    end do
                end do
            end do

            if (1 >= ido)then
                sn = 1.0_wp/(iip * l1)
                if (na /= 1) then
                    do ki=1,lid
                        cc1(1,ki,1) = sn*cc1(1,ki,1)
                        cc1(2,ki,1) = sn*cc1(2,ki,1)
                    end do
                    do j=2,iipph
                        jc = iipp2-j
                        do ki=1,lid
                            chold1 = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
                            chold2 = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
                            cc1(1,ki,j) = chold1
                            cc1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
                            cc1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
                            cc1(1,ki,jc) = chold2
                        end do
                    end do
                else
                    do ki=1,lid
                        ch1(1,ki,1) = sn*cc1(1,ki,1)
                        ch1(2,ki,1) = sn*cc1(2,ki,1)
                    end do
                    do j=2,iipph
                        jc = iipp2-j
                        do ki=1,lid
                            ch1(1,ki,j) = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
                            ch1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
                            ch1(1,ki,jc) = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
                            ch1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
                        end do
                    end do
                end if
            else
                do ki=1,lid
                    ch1(1,ki,1) = cc1(1,ki,1)
                    ch1(2,ki,1) = cc1(2,ki,1)
                end do
                do j=2,iipph
                    jc = iipp2-j
                    do ki=1,lid
                        ch1(1,ki,j) = cc1(1,ki,j)-cc1(2,ki,jc)
                        ch1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
                        ch1(1,ki,jc) = cc1(1,ki,j)+cc1(2,ki,jc)
                        ch1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
                    end do
                end do
                do i=1,ido
                    cc(1,:,1,i) = ch(1,:,i,1)
                    cc(2,:,1,i) = ch(2,:,i,1)
                end do
                do j=2,iip
                    cc(1,:,j,1) = ch(1,:,1,j)
                    cc(2,:,j,1) = ch(2,:,1,j)
                end do
                do j=2,iip
                    do i=2,ido
                        cc(1,:,j,i) = wa(i,j-1,1)*ch(1,:,i,j)+wa(i,j-1,2)*ch(2,:,i,j)
                        cc(2,:,j,i) = wa(i,j-1,1)*ch(2,:,i,j)-wa(i,j-1,2)*ch(1,:,i,j)
                    end do
                end do
            end if

        end subroutine c1fgkf


    end subroutine cfft1f


    subroutine cfft2b(ldim, l, m, c, wsave, lensav, work, lenwrk, ier)
        !
        ! cfft2b: complex 64-bit precision backward fast fourier transform, 2d.
        !
        !  Purpose:
        !
        !  cfft2b computes the two-dimensional discrete fourier transform of a
        !  complex periodic array.  this transform is known as the backward
        !  transform or fourier synthesis, transforming from spectral to
        !  physical space.  routine cfft2b is normalized, in that a call to
        !  cfft2b followed by a call to cfft2f (or vice-versa) reproduces the
        !  original array within roundoff error.
        !
        !  BUG FIX
        !  On 10 May 2010, this code was modified by changing the value
        !  of an index into the wsave array.
        !
        !  parameters:
        !
        !  INPUT
        !  integer ldim, the first dimension of c.
        !
        !  integer l, the number of elements to be transformed
        !  in the first dimension of the two-dimensional complex array c.  the value
        !  of l must be less than or equal to that of ldim.  the transform is
        !  most efficient when l is a product of small primes.
        !
        !  integer m, the number of elements to be transformed in
        !  the second dimension of the two-dimensional complex array c.  the transform
        !  is most efficient when m is a product of small primes.
        !
        !  real wsave(lensav). wsave's contents must be initialized with
        !  a call to cfft2i before the first call to routine cfft2f
        !  or cfft2b with transform lengths l and m.  wsave's contents may be
        !  re-used for subsequent calls to cfft2f and cfft2b with the same
        !  transform lengths l and m.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least
        !
        !  2*(l+m) + int(log(real(l))) + int(log(real(m))) + 8.
        !
        !  real work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*l*m.
        !
        !  INPUT/OUTPUT
        !  complex c(ldim,m), on intput, the array of two dimensions
        !  containing the (l,m) subarray to be transformed.  on
        !  output, the transformed data.
        !
        !
        !  OUTPUT
        !
        !  integer ier, the error flag.
        !
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  5, input parameter ldim < l;
        !  20, input error returned by lower level routine.
        !
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: ldim
        integer (ip), intent (in)     :: l
        integer (ip), intent (in)     :: m
        complex (wp), intent (in out) :: c(ldim,m)
        real (wp),    intent (in)     :: wsave(lensav)
        integer (ip), intent (in)     :: lensav
        real (wp),    intent (in out) :: work(lenwrk)
        integer (ip), intent (in)     :: lenwrk
        integer (ip), intent (out)    :: ier
        !------------------------------------------------------------------
        ! Dictionary: local variables
        !------------------------------------------------------------------
        integer (ip)           :: local_error_flag
        real (wp), allocatable :: real_copy(:,:,:)
        !------------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if ( ldim < l ) then
            ier = 5
            call fft_error_handler('cfft2b', -2)
            return
        else if (lensav < get_2d_saved_workspace_size(l, m)) then
            ier = 2
            call fft_error_handler('cfft2b', 6)
            return
        else if (lenwrk < get_2d_workspace_size(l, m)) then
            ier = 3
            call fft_error_handler('cfft2b', 8)
            return
        else
            ier = 0
        end if


        !
        !==> Allocate memory
        !
        allocate( real_copy(2,ldim,m) )

        ! Make a copy: complex to real
        real_copy(1,:,:) = real(c)
        real_copy(2,:,:) = aimag(c)

        !
        !==> transform x lines of real_copy
        !
        associate( &
            iw =>  get_1d_saved_workspace_size(m)-1, &
            iw1 => (l-1) + ldim*(m-1)+1, &
            iw2 => get_1d_saved_workspace_size(m), &
            iw3 => get_2d_workspace_size(l, m), &
            local_error_flag => local_error_flag &
            )

            ! Perform transform
            call cfftmb(l, 1, m, ldim, real_copy, iw1 , wsave(iw), iw2, work, iw3, local_error_flag)

        end associate

        ! Check error_flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cfft2b',-5)
            return
        end if
        !
        !==>  transform y lines of real_copy
        !
        associate( &
            iw => 1, &
            iw1 => (m-1)*ldim + l, &
            iw2 => get_1d_saved_workspace_size(l), &
            iw3 => get_2d_workspace_size(l, m), &
            local_error_flag => local_error_flag &
            )

            ! Perform transform
            call cfftmb(m, ldim, l, 1, real_copy, iw1, wsave(iw), iw2, work, iw3, local_error_flag)

        end associate

        ! Check error flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cfft2b',-5)
        end if

        ! Make copy: real to complex
        c =  cmplx(real_copy(1,:,:), real_copy(2,:,:), kind=wp)

        !
        !==> Release memory
        !
        deallocate( real_copy )

    end subroutine cfft2b


    subroutine cfft2f(ldim, l, m, c, wsave, lensav, work, lenwrk, ier)
        !
        ! cfft2f: complex 64-bit precision forward fast fourier transform, 2d.
        !
        !  Purpose:
        !
        !  cfft2f computes the two-dimensional discrete fourier transform of
        !  a complex periodic array. This transform is known as the forward
        !  transform or fourier analysis, transforming from physical to
        !  spectral space. routine cfft2f is normalized, in that a call to
        !  cfft2f followed by a call to cfft2b (or vice-versa) reproduces the
        !  original array within roundoff error.
        !
        !  BUG FIX
        !  On 10 May 2010, this code was modified by changing the value
        !  of an index into the wsave array.
        !
        !  INPUT
        !
        !  integer ldim, the first dimension of the array c.
        !
        !  integer l, the number of elements to be transformed
        !  in the first dimension of the two-dimensional complex array c.  the value
        !  of l must be less than or equal to that of ldim.  the transform is most
        !  efficient when l is a product of small primes.
        !
        !  integer m, the number of elements to be transformed
        !  in the second dimension of the two-dimensional complex array c.  the
        !  transform is most efficient when m is a product of small primes.
        !
        !  real wsave(lensav). wsave's contents must be
        !  initialized with a call to cfft2i before the first call to routine cfft2f
        !  or cfft2b with transform lengths l and m.  wsave's contents may be re-used
        !  for subsequent calls to cfft2f and cfft2b having those same
        !  transform lengths.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least
        !
        !  2*(l+m) + int(log(real(l))) + int(log(real(m))) + 8.
        !
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*l*m.
        !
        !  INPUT/OUTPUT
        !  complex c(ldim,m), on input, the array of two
        !  dimensions containing the (l,m) subarray to be transformed. On output, the
        !  transformed data.
        !
        !  real work(lenwrk), workspace array
        !
        !  OUTPUT
        !  integer ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  5, input parameter ldim < l;
        !  20, input error returned by lower level routine.
        !
        !------------------------------------------------------------------
        ! Dictionary: calling arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: ldim
        integer (ip), intent (in)     :: l
        integer (ip), intent (in)     :: m
        complex (wp), intent (in out) :: c(ldim,m)
        real (wp),    intent (in)     :: wsave(lensav)
        integer (ip), intent (in)     :: lensav
        real (wp),    intent (in out) :: work(lenwrk)
        integer (ip), intent (in)     :: lenwrk
        integer (ip), intent (out)    :: ier
        !------------------------------------------------------------------
        integer (ip)           :: local_error_flag
        real (wp), allocatable :: real_copy(:,:,:)
        !------------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if ( ldim < l ) then
            ier = 5
            call fft_error_handler('cfft2f', -2)
            return
        else if (lensav < get_2d_saved_workspace_size(l, m)) then
            ier = 2
            call fft_error_handler('cfft2f', 6)
            return
        else if (lenwrk < get_2d_workspace_size(l, m)) then
            ier = 3
            call fft_error_handler('cfft2f', 8)
            return
        else
            ier = 0
        end if

        !
        !==> Allocate memory
        !
        allocate( real_copy(2,ldim,m) )

        ! Make a copy: complex to real array
        real_copy(1,:,:)=real(c)
        real_copy(2,:,:)=aimag(c)

        !
        !==> Transform x lines of real_copy
        !
        associate( &
            iw =>  2*l+int(log(real(l, kind=wp) )/log(2.0_wp)) + 3, &
            iw1 => (l-1) + ldim*(m-1) +1, &
            iw2 => get_1d_saved_workspace_size(m), &
            iw3 => get_2d_workspace_size(l, m), &
            local_error_flag => local_error_flag &
            )

            call cfftmf(l, 1, m, ldim, real_copy, iw1, wsave(iw), iw2 , work, iw3, local_error_flag)

        end associate

        ! Check error flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cfft2f',-5)
            return
        end if
        !
        !==> Transform y lines of real_copy
        !
        associate( &
            iw => 1, &
            iw1 => (m-1)*ldim + l, &
            iw2 => get_1d_saved_workspace_size(l), &
            iw3 => 2*m*l, &
            local_error_flag => local_error_flag &
            )

            call cfftmf(m, ldim, l, 1, real_copy, iw1, wsave(iw), iw2, work, iw3, local_error_flag)

        end associate

        ! Check error flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cfft2f',-5)
        end if

        ! Make copy: real to complex
        c =  cmplx(real_copy(1,:,:), real_copy(2,:,:), kind=wp)

        !
        !==> Release memory
        !
        deallocate( real_copy )

    end subroutine cfft2f


    subroutine cfft2i(l, m, wsave, lensav, ier)
        !
        ! cfft2i: initialization for cfft2b and cfft2f.
        !
        !  purpose:
        !
        !  cfft2i initializes real array wsave for use in its companion
        !  routines cfft2f and cfft2b for computing two-dimensional fast
        !  fourier transforms of complex data.  prime factorizations of l and m,
        !  together with tabulations of the trigonometric functions, are
        !  computed and stored in array wsave.
        !
        !  on 10 may 2010, this code was modified by changing the value
        !  of an index into the wsave array.
        !
        !  INPUT
        !  integer l, the number of elements to be transformed
        !  in the first dimension. The transform is most efficient when l is a
        !  product of small primes.
        !
        !  integer m, the number of elements to be transformed
        !  in the second dimension. The transform is most efficient when m is a
        !  product of small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least
        !
        !  2*(l+m) + int(log(real(l)))+ int(log(real(m))) + 8.
        !
        !  OUTPUT
        !  real wsave(lensav), contains the prime factors of l
        !  and m, and also certain trigonometric values which will be used in
        !  routines cfft2b or cfft2f.
        !
        !  integer  ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  20, input error returned by lower level routine.
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: l
        integer (ip), intent (in)  :: m
        real (wp),    intent (out) :: wsave(lensav)
        integer (ip), intent (in)  :: lensav
        integer (ip), intent (out) :: ier
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip) local_error_flag
        !--------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if ( lensav < get_2d_saved_workspace_size(l, m)) then
            ier = 2
            call fft_error_handler('cfft2i', 4)
            return
        else
            ier = 0
        end if

        associate( temp => get_1d_saved_workspace_size(l) )

            call cfftmi(l, wsave(1), temp, local_error_flag)

        end associate

        ! Check error flag
        if ( local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cfft2i',-5)
            return
        end if

        associate( &
            iw => get_1d_saved_workspace_size(l) - 1, &
            temp => get_1d_saved_workspace_size(m) &
            )

            call cfftmi(m, wsave(iw), temp, local_error_flag)

        end associate

        ! Check error flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cfft2i',-5)
        end if

    end subroutine cfft2i



    subroutine cfftmb(lot, jump, n, inc, c, lenc, wsave, lensav, work, &
        lenwrk, ier)
        !
        ! cfftmb: complex 64-bit precision backward fft, 1d, multiple vectors.
        !
        !  purpose:
        !
        !  cfftmb computes the one-dimensional fourier transform of multiple
        !  periodic sequences within a complex array.  this transform is referred
        !  to as the backward transform or fourier synthesis, transforming the
        !  sequences from spectral to physical space.  this transform is
        !  normalized since a call to cfftmf followed by a call to cfftmb (or
        !  vice-versa) reproduces the original array within roundoff error.
        !
        !  the parameters inc, jump, n and lot are consistent if equality
        !  i1*inc + j1*jump = i2*inc + j2*jump for i1,i2 < n and j1,j2 < lot
        !  implies i1=i2 and j1=j2.  for multiple ffts to execute correctly,
        !  input variables inc, jump, n and lot must be consistent, otherwise
        !  at least one array element mistakenly is transformed more than once.
        !
        !  parameters:
        !
        !  integer lot, the number of sequences to be transformed
        !  within array c.
        !
        !  integer jump, the increment between the locations, in
        !  array c, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array c, of two consecutive elements within the same sequence to be
        !  transformed.
        !
        !  input/output, complex (wp) c(lenc), an array containing lot
        !  sequences, each having length n, to be transformed.  c can have any
        !  number of dimensions, but the total number of locations must be at least
        !  lenc.  on output, c contains the transformed sequences.
        !
        !  integer lenc, the dimension of the c array.
        !  lenc must be at least (lot-1)*jump + inc*(n-1) + 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to cfftmi before the first call to routine cfftmf
        !  or cfftmb for a given transform length n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*lot*n.
        !
        !  integer ier, error flag.
        !  0, successful exit
        !  1, input parameter lenc not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  4, input parameters inc, jump, n, lot are not consistent.
        !
        integer (ip) lenc
        integer (ip) lensav
        integer (ip) lenwrk
        real (wp) c(2,lenc)
        integer (ip) ier
        integer (ip) inc
        integer (ip) iw1
        integer (ip) jump
        integer (ip) lot
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lenc < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cfftmb ', 6)
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('cfftmb ', 8)
        else if (lenwrk < 2*lot*n) then
            ier = 3
            call fft_error_handler('cfftmb ', 10)
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('cfftmb ', -1)
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then

            ! Set workspace pointer
            iw1 = 2*n+1

            call cmfm1b(lot,jump,n,inc,c,work,wsave,wsave(iw1),wsave(iw1+1))

        end if

    end subroutine cfftmb



    subroutine cfftmf(lot, jump, n, inc, c, lenc, wsave, lensav, work, &
        lenwrk, ier)
        !
        ! cfftmf: complex 64-bit precision forward fft, 1d, multiple vectors.
        !
        !  Purpose:
        !
        !  cfftmf computes the one-dimensional fourier transform of multiple
        !  periodic sequences within a complex array. this transform is referred
        !  to as the forward transform or fourier analysis, transforming the
        !  sequences from physical to spectral space. this transform is
        !  normalized since a call to cfftmf followed by a call to cfftmb
        !  (or vice-versa) reproduces the original array within roundoff error.
        !
        !  the parameters integers inc, jump, n and lot are consistent if equality
        !  i1*inc + j1*jump = i2*inc + j2*jump for i1,i2 < n and j1,j2 < lot
        !  implies i1=i2 and j1=j2. for multiple ffts to execute correctly,
        !  input variables inc, jump, n and lot must be consistent, otherwise
        !  at least one array element mistakenly is transformed more than once.
        !
        !  parameters:
        !
        !  integer lot, the number of sequences to be
        !  transformed within array c.
        !
        !  integer jump, the increment between the locations,
        !  in array c, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array c, of two consecutive elements within the same sequence to be
        !  transformed.
        !
        !  input/output, complex (wp) c(lenc), array containing lot sequences,
        !  each having length n, to be transformed.  c can have any number of
        !  dimensions, but the total number of locations must be at least lenc.
        !
        !  integer lenc, the dimension of the c array.
        !  lenc must be at least (lot-1)*jump + inc*(n-1) + 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to cfftmi before the first call to routine cfftmf
        !  or cfftmb for a given transform length n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*lot*n.
        !
        !  integer ier, error flag.
        !  0 successful exit;
        !  1 input parameter lenc not big enough;
        !  2 input parameter lensav not big enough;
        !  3 input parameter lenwrk not big enough;
        !  4 input parameters inc, jump, n, lot are not consistent.
        !


        integer (ip) lenc
        integer (ip) lensav
        integer (ip) lenwrk

        real(wp) c(2,lenc)
        integer (ip) ier
        integer (ip) inc
        integer (ip) iw1
        integer (ip) jump
        integer (ip) lot
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lenc < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cfftmf ', 6)
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('cfftmf ', 8)
        else if (lenwrk < 2*lot*n) then
            ier = 3
            call fft_error_handler('cfftmf ', 10)
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('cfftmf ', -1)
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then

            ! Set workspace index pointer
            iw1 = 2*n+1
            call cmfm1f(lot,jump,n,inc,c,work,wsave,wsave(iw1),wsave(iw1+1))
        end if

    end subroutine cfftmf


    subroutine cfftmi(n, wsave, lensav, ier)


        !
        !! cfftmi: initialization for cfftmb and cfftmf.
        !
        !  purpose:
        !
        !  cfftmi initializes array wsave for use in its companion routines
        !  cfftmb and cfftmf.  cfftmi must be called before the first call
        !  to cfftmb or cfftmf, and after whenever the value of integer n changes.
        !
        !  parameters:
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  real wsave(lensav), containing the prime factors
        !  of n and also containing certain trigonometric values which will be used in
        !  routines cfftmb or cfftmf.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough.
        !


        integer (ip) lensav

        integer (ip) ier
        integer (ip) iw1
        integer (ip) n
        real (wp) wsave(lensav)

        ier = 0

        if (lensav < get_1d_saved_workspace_size(n) ) then
            ier = 2
            call fft_error_handler('cfftmi ', 3)
        end if

        if (n /= 1) then
            iw1 = n+n+1
            call mcfti1(n,wsave,wsave(iw1),wsave(iw1+1))
        end if

    end subroutine cfftmi


    subroutine cmfm1b(lot, jump, n, inc, c, ch, wa, fnf, fac)

        real (wp) c(2,*)
        real (wp) ch(*)
        real (wp) fac(*)
        real (wp) fnf
        integer (ip) ido
        integer (ip) inc
        integer (ip) iip
        integer (ip) iw
        integer (ip) jump
        integer (ip) k1
        integer (ip) l1
        integer (ip) l2
        integer (ip) lid
        integer (ip) lot
        integer (ip) n
        integer (ip) na
        integer (ip) nbr
        integer (ip) nf
        real (wp) wa(*)

        nf = int(fnf, kind=ip)
        na = 0
        l1 = 1
        iw = 1
        do k1=1,nf
            iip = int(fac(k1), kind=ip)
            l2 = iip*l1
            ido = n/l2
            lid = l1*ido
            nbr = 1+na+2*min(iip-2,4)

            select case (nbr)
                case (1)
                    call cmf2kb(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (2)
                    call cmf2kb(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (3)
                    call cmf3kb(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (4)
                    call cmf3kb(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (5)
                    call cmf4kb(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (6)
                    call cmf4kb(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (7)
                    call cmf5kb(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (8)
                    call cmf5kb(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (9)
                    call cmfgkb(lot,ido,iip,l1,lid,na,c,c,jump,inc,ch,ch,1,lot,wa(iw))
                case (10)
                    call cmfgkb(lot,ido,iip,l1,lid,na,ch,ch,1,lot,c,c, &
                        jump,inc,wa(iw))
            end select

            l1 = l2
            iw = iw+(iip-1)*(2*ido)

            if (iip <= 5) then
                na = 1-na
            end if

        end do


    contains


        subroutine cmf2kb(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(2,in1,l1,ido,2)
            real (wp) ch(2,in2,l1,2,ido)
            real (wp) chold1
            real (wp) chold2
            integer (ip) i
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp) ti2
            real (wp) tr2
            real (wp) wa(ido,1,2)

            m1d = (lot-1)*im1+1
            m2s = 1-im2

            if (1 >= ido .and. na /= 1) then
                do k=1,l1
                    do m1=1,m1d,im1
                        chold1 = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
                        cc(1,m1,k,1,2) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
                        cc(1,m1,k,1,1) = chold1
                        chold2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
                        cc(2,m1,k,1,2) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
                        cc(2,m1,k,1,1) = chold2
                    end do
                end do
            else
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
                        ch(1,m2,k,2,1) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
                        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
                        ch(2,m2,k,2,1) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
                    end do
                end do

                do i=2,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+cc(1,m1,k,i,2)
                            tr2 = cc(1,m1,k,i,1)-cc(1,m1,k,i,2)
                            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+cc(2,m1,k,i,2)
                            ti2 = cc(2,m1,k,i,1)-cc(2,m1,k,i,2)
                            ch(2,m2,k,2,i) = wa(i,1,1)*ti2+wa(i,1,2)*tr2
                            ch(1,m2,k,2,i) = wa(i,1,1)*tr2-wa(i,1,2)*ti2
                        end do
                    end do
                end do
            end if

        end subroutine cmf2kb

        subroutine cmf3kb(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(2,in1,l1,ido,3)
            real (wp) ch(2,in2,l1,3,ido)
            real (wp) ci2
            real (wp) ci3
            real (wp) cr2
            real (wp) cr3
            real (wp) di2
            real (wp) di3
            real (wp) dr2
            real (wp) dr3
            integer (ip) i
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp), parameter :: TAUI =  0.866025403784439_wp
            real (wp), parameter :: TAUR = -0.5_wp
            real (wp) ti2
            real (wp) tr2
            real (wp) wa(ido,2,2)

            m1d = (lot-1)*im1+1
            m2s = 1-im2

            if (1 >= ido .and. na /= 1) then
                do k=1,l1
                    do m1=1,m1d,im1
                        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
                        cr2 = cc(1,m1,k,1,1)+TAUR*tr2
                        cc(1,m1,k,1,1) = cc(1,m1,k,1,1)+tr2
                        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
                        ci2 = cc(2,m1,k,1,1)+TAUR*ti2
                        cc(2,m1,k,1,1) = cc(2,m1,k,1,1)+ti2
                        cr3 = TAUI*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
                        ci3 = TAUI*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
                        cc(1,m1,k,1,2) = cr2-ci3
                        cc(1,m1,k,1,3) = cr2+ci3
                        cc(2,m1,k,1,2) = ci2+cr3
                        cc(2,m1,k,1,3) = ci2-cr3
                    end do
                end do
            else
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
                        cr2 = cc(1,m1,k,1,1)+TAUR*tr2
                        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2
                        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
                        ci2 = cc(2,m1,k,1,1)+TAUR*ti2
                        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2
                        cr3 = TAUI*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
                        ci3 = TAUI*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
                        ch(1,m2,k,2,1) = cr2-ci3
                        ch(1,m2,k,3,1) = cr2+ci3
                        ch(2,m2,k,2,1) = ci2+cr3
                        ch(2,m2,k,3,1) = ci2-cr3
                    end do
                end do

                do i=2,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,3)
                            cr2 = cc(1,m1,k,i,1)+TAUR*tr2
                            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2
                            ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,3)
                            ci2 = cc(2,m1,k,i,1)+TAUR*ti2
                            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2
                            cr3 = TAUI*(cc(1,m1,k,i,2)-cc(1,m1,k,i,3))
                            ci3 = TAUI*(cc(2,m1,k,i,2)-cc(2,m1,k,i,3))
                            dr2 = cr2-ci3
                            dr3 = cr2+ci3
                            di2 = ci2+cr3
                            di3 = ci2-cr3
                            ch(2,m2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
                            ch(1,m2,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
                            ch(2,m2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
                            ch(1,m2,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
                        end do
                    end do
                end do
            end if
        end subroutine cmf3kb

        subroutine cmf4kb(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(2,in1,l1,ido,4)
            real (wp) ch(2,in2,l1,4,ido)
            real (wp) ci2
            real (wp) ci3
            real (wp) ci4
            real (wp) cr2
            real (wp) cr3
            real (wp) cr4
            integer (ip) i
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp) ti1
            real (wp) ti2
            real (wp) ti3
            real (wp) ti4
            real (wp) tr1
            real (wp) tr2
            real (wp) tr3
            real (wp) tr4
            real (wp) wa(ido,3,2)

            m1d = (lot-1)*im1+1
            m2s = 1-im2

            if (1 >= ido .and. na /= 1) then
                do k=1,l1
                    do m1=1,m1d,im1
                        ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
                        ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
                        tr4 = cc(2,m1,k,1,4)-cc(2,m1,k,1,2)
                        ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
                        tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
                        tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
                        ti4 = cc(1,m1,k,1,2)-cc(1,m1,k,1,4)
                        tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
                        cc(1,m1,k,1,1) = tr2+tr3
                        cc(1,m1,k,1,3) = tr2-tr3
                        cc(2,m1,k,1,1) = ti2+ti3
                        cc(2,m1,k,1,3) = ti2-ti3
                        cc(1,m1,k,1,2) = tr1+tr4
                        cc(1,m1,k,1,4) = tr1-tr4
                        cc(2,m1,k,1,2) = ti1+ti4
                        cc(2,m1,k,1,4) = ti1-ti4
                    end do
                end do
            else
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
                        ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
                        tr4 = cc(2,m1,k,1,4)-cc(2,m1,k,1,2)
                        ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
                        tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
                        tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
                        ti4 = cc(1,m1,k,1,2)-cc(1,m1,k,1,4)
                        tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
                        ch(1,m2,k,1,1) = tr2+tr3
                        ch(1,m2,k,3,1) = tr2-tr3
                        ch(2,m2,k,1,1) = ti2+ti3
                        ch(2,m2,k,3,1) = ti2-ti3
                        ch(1,m2,k,2,1) = tr1+tr4
                        ch(1,m2,k,4,1) = tr1-tr4
                        ch(2,m2,k,2,1) = ti1+ti4
                        ch(2,m2,k,4,1) = ti1-ti4
                    end do
                end do

                do i=2,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ti1 = cc(2,m1,k,i,1)-cc(2,m1,k,i,3)
                            ti2 = cc(2,m1,k,i,1)+cc(2,m1,k,i,3)
                            ti3 = cc(2,m1,k,i,2)+cc(2,m1,k,i,4)
                            tr4 = cc(2,m1,k,i,4)-cc(2,m1,k,i,2)
                            tr1 = cc(1,m1,k,i,1)-cc(1,m1,k,i,3)
                            tr2 = cc(1,m1,k,i,1)+cc(1,m1,k,i,3)
                            ti4 = cc(1,m1,k,i,2)-cc(1,m1,k,i,4)
                            tr3 = cc(1,m1,k,i,2)+cc(1,m1,k,i,4)
                            ch(1,m2,k,1,i) = tr2+tr3
                            cr3 = tr2-tr3
                            ch(2,m2,k,1,i) = ti2+ti3
                            ci3 = ti2-ti3
                            cr2 = tr1+tr4
                            cr4 = tr1-tr4
                            ci2 = ti1+ti4
                            ci4 = ti1-ti4
                            ch(1,m2,k,2,i) = wa(i,1,1)*cr2-wa(i,1,2)*ci2
                            ch(2,m2,k,2,i) = wa(i,1,1)*ci2+wa(i,1,2)*cr2
                            ch(1,m2,k,3,i) = wa(i,2,1)*cr3-wa(i,2,2)*ci3
                            ch(2,m2,k,3,i) = wa(i,2,1)*ci3+wa(i,2,2)*cr3
                            ch(1,m2,k,4,i) = wa(i,3,1)*cr4-wa(i,3,2)*ci4
                            ch(2,m2,k,4,i) = wa(i,3,1)*ci4+wa(i,3,2)*cr4
                        end do
                    end do
                end do
            end if

        end subroutine cmf4kb

        subroutine cmf5kb(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(2,in1,l1,ido,5)
            real (wp) ch(2,in2,l1,5,ido)
            real (wp) chold1
            real (wp) chold2
            real (wp) ci2
            real (wp) ci3
            real (wp) ci4
            real (wp) ci5
            real (wp) cr2
            real (wp) cr3
            real (wp) cr4
            real (wp) cr5
            real (wp) di2
            real (wp) di3
            real (wp) di4
            real (wp) di5
            real (wp) dr2
            real (wp) dr3
            real (wp) dr4
            real (wp) dr5
            integer (ip) i
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp) ti2
            real (wp) ti3
            real (wp) ti4
            real (wp) ti5
            real (wp) tr2
            real (wp) tr3
            real (wp) tr4
            real (wp) tr5
            real (wp), parameter :: SQRT5 = sqrt(5.0_wp)
            real (wp), parameter :: SQRT5_PLUS_5 = SQRT5 + 5.0_wp
            real (wp), parameter :: TI11 = sqrt(SQRT5_PLUS_5/2)/2             ! 0.9510565162951536_wp
            real (wp), parameter :: TI12 = sqrt(5.0_wp/(2.0_wp*SQRT5_PLUS_5)) ! 0.5877852522924731_wp
            real (wp), parameter :: TR11 =  (SQRT5 - 1.0_wp)/4                 ! 0.3090169943749474_wp
            real (wp), parameter :: TR12 = -(1.0_wp + SQRT5)/4                 !-0.8090169943749474_wp

            real (wp) wa(ido,4,2)

            m1d = (lot-1)*im1+1
            m2s = 1-im2

            if (1 >= ido .and. na /= 1) then
                do k=1,l1
                    do m1=1,m1d,im1
                        ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
                        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
                        ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
                        ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
                        tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
                        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
                        tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
                        tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
                        chold1 = cc(1,m1,k,1,1)+tr2+tr3
                        chold2 = cc(2,m1,k,1,1)+ti2+ti3
                        cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
                        ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
                        cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
                        ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
                        cc(1,m1,k,1,1) = chold1
                        cc(2,m1,k,1,1) = chold2
                        cr5 = ti11*tr5+ti12*tr4
                        ci5 = ti11*ti5+ti12*ti4
                        cr4 = ti12*tr5-ti11*tr4
                        ci4 = ti12*ti5-ti11*ti4
                        cc(1,m1,k,1,2) = cr2-ci5
                        cc(1,m1,k,1,5) = cr2+ci5
                        cc(2,m1,k,1,2) = ci2+cr5
                        cc(2,m1,k,1,3) = ci3+cr4
                        cc(1,m1,k,1,3) = cr3-ci4
                        cc(1,m1,k,1,4) = cr3+ci4
                        cc(2,m1,k,1,4) = ci3-cr4
                        cc(2,m1,k,1,5) = ci2-cr5
                    end do
                end do
            else
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
                        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
                        ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
                        ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
                        tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
                        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
                        tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
                        tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
                        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2+tr3
                        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2+ti3
                        cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
                        ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
                        cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
                        ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
                        cr5 = ti11*tr5+ti12*tr4
                        ci5 = ti11*ti5+ti12*ti4
                        cr4 = ti12*tr5-ti11*tr4
                        ci4 = ti12*ti5-ti11*ti4
                        ch(1,m2,k,2,1) = cr2-ci5
                        ch(1,m2,k,5,1) = cr2+ci5
                        ch(2,m2,k,2,1) = ci2+cr5
                        ch(2,m2,k,3,1) = ci3+cr4
                        ch(1,m2,k,3,1) = cr3-ci4
                        ch(1,m2,k,4,1) = cr3+ci4
                        ch(2,m2,k,4,1) = ci3-cr4
                        ch(2,m2,k,5,1) = ci2-cr5
                    end do
                end do

                do i=2,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ti5 = cc(2,m1,k,i,2)-cc(2,m1,k,i,5)
                            ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,5)
                            ti4 = cc(2,m1,k,i,3)-cc(2,m1,k,i,4)
                            ti3 = cc(2,m1,k,i,3)+cc(2,m1,k,i,4)
                            tr5 = cc(1,m1,k,i,2)-cc(1,m1,k,i,5)
                            tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,5)
                            tr4 = cc(1,m1,k,i,3)-cc(1,m1,k,i,4)
                            tr3 = cc(1,m1,k,i,3)+cc(1,m1,k,i,4)
                            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2+tr3
                            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2+ti3
                            cr2 = cc(1,m1,k,i,1)+tr11*tr2+tr12*tr3
                            ci2 = cc(2,m1,k,i,1)+tr11*ti2+tr12*ti3
                            cr3 = cc(1,m1,k,i,1)+tr12*tr2+tr11*tr3
                            ci3 = cc(2,m1,k,i,1)+tr12*ti2+tr11*ti3
                            cr5 = ti11*tr5+ti12*tr4
                            ci5 = ti11*ti5+ti12*ti4
                            cr4 = ti12*tr5-ti11*tr4
                            ci4 = ti12*ti5-ti11*ti4
                            dr3 = cr3-ci4
                            dr4 = cr3+ci4
                            di3 = ci3+cr4
                            di4 = ci3-cr4
                            dr5 = cr2+ci5
                            dr2 = cr2-ci5
                            di5 = ci2-cr5
                            di2 = ci2+cr5
                            ch(1,m2,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
                            ch(2,m2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
                            ch(1,m2,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
                            ch(2,m2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
                            ch(1,m2,k,4,i) = wa(i,3,1)*dr4-wa(i,3,2)*di4
                            ch(2,m2,k,4,i) = wa(i,3,1)*di4+wa(i,3,2)*dr4
                            ch(1,m2,k,5,i) = wa(i,4,1)*dr5-wa(i,4,2)*di5
                            ch(2,m2,k,5,i) = wa(i,4,1)*di5+wa(i,4,2)*dr5
                        end do
                    end do
                end do
            end if

        end subroutine cmf5kb

        subroutine cmfgkb(lot, ido, iip, l1, lid, na, cc, cc1, im1, in1, &
            ch, ch1, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) iip
            integer (ip) l1
            integer (ip) lid

            real (wp) cc(2,in1,l1,iip,ido)
            real (wp) cc1(2,in1,lid,iip)
            real (wp) ch(2,in2,l1,ido,iip)
            real (wp) ch1(2,in2,lid,iip)
            real (wp) chold1
            real (wp) chold2
            integer (ip) i
            integer (ip) idlj
            integer (ip) im1
            integer (ip) im2
            integer (ip) iipp2
            integer (ip) iipph
            integer (ip) j
            integer (ip) jc
            integer (ip) k
            integer (ip) ki
            integer (ip) l
            integer (ip) lc
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp) wa(ido,iip-1,2)
            real (wp) wai
            real (wp) war

            m1d = (lot-1)*im1+1
            m2s = 1-im2
            iipp2 = iip+2
            iipph = (iip+1)/2

            do ki=1,lid
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
                    ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
                end do
            end do

            do j=2,iipph
                jc = iipp2-j
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch1(1,m2,ki,j) =  cc1(1,m1,ki,j)+cc1(1,m1,ki,jc)
                        ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)-cc1(1,m1,ki,jc)
                        ch1(2,m2,ki,j) =  cc1(2,m1,ki,j)+cc1(2,m1,ki,jc)
                        ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(2,m1,ki,jc)
                    end do
                end do
            end do

            do j=2,iipph
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        cc1(1,m1,ki,1) = cc1(1,m1,ki,1)+ch1(1,m2,ki,j)
                        cc1(2,m1,ki,1) = cc1(2,m1,ki,1)+ch1(2,m2,ki,j)
                    end do
                end do
            end do

            do l=2,iipph
                lc = iipp2-l
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        cc1(1,m1,ki,l) = ch1(1,m2,ki,1)+wa(1,l-1,1)*ch1(1,m2,ki,2)
                        cc1(1,m1,ki,lc) = wa(1,l-1,2)*ch1(1,m2,ki,iip)
                        cc1(2,m1,ki,l) = ch1(2,m2,ki,1)+wa(1,l-1,1)*ch1(2,m2,ki,2)
                        cc1(2,m1,ki,lc) = wa(1,l-1,2)*ch1(2,m2,ki,iip)
                    end do
                end do
                do j=3,iipph
                    jc = iipp2-j
                    idlj = mod((l-1)*(j-1),iip)
                    war = wa(1,idlj,1)
                    wai = wa(1,idlj,2)
                    do ki=1,lid
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc1(1,m1,ki,l) = cc1(1,m1,ki,l)+war*ch1(1,m2,ki,j)
                            cc1(1,m1,ki,lc) = cc1(1,m1,ki,lc)+wai*ch1(1,m2,ki,jc)
                            cc1(2,m1,ki,l) = cc1(2,m1,ki,l)+war*ch1(2,m2,ki,j)
                            cc1(2,m1,ki,lc) = cc1(2,m1,ki,lc)+wai*ch1(2,m2,ki,jc)
                        end do
                    end do
                end do
            end do

            if (1 >= ido .and. na /= 1) then
                do j=2,iipph
                    jc = iipp2-j
                    do ki=1,lid
                        do m1=1,m1d,im1
                            chold1 = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
                            chold2 = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
                            cc1(1,m1,ki,j) = chold1
                            cc1(2,m1,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
                            cc1(2,m1,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
                            cc1(1,m1,ki,jc) = chold2
                        end do
                    end do
                end do
            else
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
                        ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
                    end do
                end do

                do j=2,iipph
                    jc = iipp2-j
                    do ki=1,lid
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch1(1,m2,ki,j) = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
                            ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
                            ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
                            ch1(2,m2,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
                        end do
                    end do
                end do

                if (ido /= 1) then
                    do i=1,ido
                        do k=1,l1
                            m2 = m2s
                            do m1=1,m1d,im1
                                m2 = m2+im2
                                cc(1,m1,k,1,i) = ch(1,m2,k,i,1)
                                cc(2,m1,k,1,i) = ch(2,m2,k,i,1)
                            end do
                        end do
                    end do

                    do j=2,iip
                        do k=1,l1
                            m2 = m2s
                            do m1=1,m1d,im1
                                m2 = m2+im2
                                cc(1,m1,k,j,1) = ch(1,m2,k,1,j)
                                cc(2,m1,k,j,1) = ch(2,m2,k,1,j)
                            end do
                        end do
                    end do

                    do j=2,iip
                        do i=2,ido
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    cc(1,m1,k,j,i) = wa(i,j-1,1)*ch(1,m2,k,i,j) &
                                        -wa(i,j-1,2)*ch(2,m2,k,i,j)
                                    cc(2,m1,k,j,i) = wa(i,j-1,1)*ch(2,m2,k,i,j) &
                                        +wa(i,j-1,2)*ch(1,m2,k,i,j)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

        end subroutine cmfgkb

    end subroutine cmfm1b

    subroutine cmfm1f(lot, jump, n, inc, c, ch, wa, fnf, fac)

        real (wp) c(2,*)
        real (wp) ch(*)
        real (wp) fac(*)
        real (wp) fnf
        integer (ip) ido
        integer (ip) inc
        integer (ip) iip
        integer (ip) iw
        integer (ip) jump
        integer (ip) k1
        integer (ip) l1
        integer (ip) l2
        integer (ip) lid
        integer (ip) lot
        integer (ip) n
        integer (ip) na
        integer (ip) nbr
        integer (ip) nf
        real (wp) wa(*)

        nf = int(fnf, kind=ip)
        na = 0
        l1 = 1
        iw = 1

        do k1=1,nf
            iip = int(fac(k1), kind=ip)
            l2 = iip*l1
            ido = n/l2
            lid = l1*ido
            nbr = 1+na+2*min(iip-2,4)
            select case (nbr)
                case (1)
                    call cmf2kf(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (2)
                    call cmf2kf(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (3)
                    call cmf3kf(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (4)
                    call cmf3kf(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (5)
                    call cmf4kf(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (6)
                    call cmf4kf(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (7)
                    call cmf5kf(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (8)
                    call cmf5kf(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (9)
                    call cmfgkf(lot,ido,iip,l1,lid,na,c,c,jump,inc,ch,ch,1,lot,wa(iw))
                case (10)
                    call cmfgkf(lot,ido,iip,l1,lid,na,ch,ch,1,lot,c,c, &
                        jump,inc,wa(iw))
            end select

            l1 = l2
            iw = iw+(iip-1)*(2*ido)

            if (iip <= 5) then
                na = 1-na
            end if
        end do

    contains

        subroutine cmf2kf(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)
            !--------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------
            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1
            real (wp) cc(2,in1,l1,ido,2)

            real (wp) ch(2,in2,l1,2,ido)
            real (wp) chold1
            real (wp) chold2

            integer (ip) i
            integer (ip) im1
            integer (ip) im2
            integer (ip) k

            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s

            integer (ip) na
            real (wp) sn
            real (wp) ti2
            real (wp) tr2
            real (wp) wa(ido,1,2)

            m1d = (lot-1)*im1+1
            m2s = 1-im2

            if (1 >= ido) then
                sn = 1.0_wp/(2 * l1)
                if (na /= 1) then
                    do k=1,l1
                        do m1=1,m1d,im1
                            chold1 = sn*(cc(1,m1,k,1,1)+cc(1,m1,k,1,2))
                            cc(1,m1,k,1,2) = sn*(cc(1,m1,k,1,1)-cc(1,m1,k,1,2))
                            cc(1,m1,k,1,1) = chold1
                            chold2 = sn*(cc(2,m1,k,1,1)+cc(2,m1,k,1,2))
                            cc(2,m1,k,1,2) = sn*(cc(2,m1,k,1,1)-cc(2,m1,k,1,2))
                            cc(2,m1,k,1,1) = chold2
                        end do
                    end do
                else
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(1,m2,k,1,1) = sn*(cc(1,m1,k,1,1)+cc(1,m1,k,1,2))
                            ch(1,m2,k,2,1) = sn*(cc(1,m1,k,1,1)-cc(1,m1,k,1,2))
                            ch(2,m2,k,1,1) = sn*(cc(2,m1,k,1,1)+cc(2,m1,k,1,2))
                            ch(2,m2,k,2,1) = sn*(cc(2,m1,k,1,1)-cc(2,m1,k,1,2))
                        end do
                    end do
                end if
            else
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
                        ch(1,m2,k,2,1) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
                        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
                        ch(2,m2,k,2,1) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
                    end do
                end do

                do i=2,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+cc(1,m1,k,i,2)
                            tr2 = cc(1,m1,k,i,1)-cc(1,m1,k,i,2)
                            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+cc(2,m1,k,i,2)
                            ti2 = cc(2,m1,k,i,1)-cc(2,m1,k,i,2)
                            ch(2,m2,k,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
                            ch(1,m2,k,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
                        end do
                    end do
                end do
            end if

        end subroutine cmf2kf



        subroutine cmf3kf(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(2,in1,l1,ido,3)
            real (wp) ch(2,in2,l1,3,ido)
            real (wp) ci2
            real (wp) ci3
            real (wp) cr2
            real (wp) cr3
            real (wp) di2
            real (wp) di3
            real (wp) dr2
            real (wp) dr3
            integer (ip) i
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp) sn
            real (wp), parameter :: taui = -0.866025403784439_wp
            real (wp), parameter :: taur = -0.5_wp
            real (wp) ti2
            real (wp) tr2
            real (wp) wa(ido,2,2)

            m1d = (lot-1)*im1+1
            m2s = 1-im2

            if (1 >= ido) then
                sn = 1.0_wp/(3 * l1)
                if (na /= 1) then
                    do k=1,l1
                        do m1=1,m1d,im1
                            tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
                            cr2 = cc(1,m1,k,1,1)+taur*tr2
                            cc(1,m1,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2)
                            ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
                            ci2 = cc(2,m1,k,1,1)+taur*ti2
                            cc(2,m1,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2)
                            cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
                            ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
                            cc(1,m1,k,1,2) = sn*(cr2-ci3)
                            cc(1,m1,k,1,3) = sn*(cr2+ci3)
                            cc(2,m1,k,1,2) = sn*(ci2+cr3)
                            cc(2,m1,k,1,3) = sn*(ci2-cr3)
                        end do
                    end do
                else
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
                            cr2 = cc(1,m1,k,1,1)+taur*tr2
                            ch(1,m2,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2)
                            ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
                            ci2 = cc(2,m1,k,1,1)+taur*ti2
                            ch(2,m2,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2)
                            cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
                            ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
                            ch(1,m2,k,2,1) = sn*(cr2-ci3)
                            ch(1,m2,k,3,1) = sn*(cr2+ci3)
                            ch(2,m2,k,2,1) = sn*(ci2+cr3)
                            ch(2,m2,k,3,1) = sn*(ci2-cr3)
                        end do
                    end do
                end if
            else
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
                        cr2 = cc(1,m1,k,1,1)+taur*tr2
                        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2
                        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
                        ci2 = cc(2,m1,k,1,1)+taur*ti2
                        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2
                        cr3 = taui*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
                        ci3 = taui*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
                        ch(1,m2,k,2,1) = cr2-ci3
                        ch(1,m2,k,3,1) = cr2+ci3
                        ch(2,m2,k,2,1) = ci2+cr3
                        ch(2,m2,k,3,1) = ci2-cr3
                    end do
                end do
                do i=2,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,3)
                            cr2 = cc(1,m1,k,i,1)+taur*tr2
                            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2
                            ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,3)
                            ci2 = cc(2,m1,k,i,1)+taur*ti2
                            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2
                            cr3 = taui*(cc(1,m1,k,i,2)-cc(1,m1,k,i,3))
                            ci3 = taui*(cc(2,m1,k,i,2)-cc(2,m1,k,i,3))
                            dr2 = cr2-ci3
                            dr3 = cr2+ci3
                            di2 = ci2+cr3
                            di3 = ci2-cr3
                            ch(2,m2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
                            ch(1,m2,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
                            ch(2,m2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
                            ch(1,m2,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
                        end do
                    end do
                end do
            end if

        end subroutine cmf3kf



        subroutine cmf4kf(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(2,in1,l1,ido,4)
            real (wp) ch(2,in2,l1,4,ido)
            real (wp) ci2
            real (wp) ci3
            real (wp) ci4
            real (wp) cr2
            real (wp) cr3
            real (wp) cr4
            integer (ip) i
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp) sn
            real (wp) ti1
            real (wp) ti2
            real (wp) ti3
            real (wp) ti4
            real (wp) tr1
            real (wp) tr2
            real (wp) tr3
            real (wp) tr4
            real (wp) wa(ido,3,2)

            m1d = (lot-1)*im1+1
            m2s = 1-im2

            if (1 >= ido) then
                sn = 1.0_wp /(4 * l1)
                if (na /= 1) then
                    do k=1,l1
                        do m1=1,m1d,im1
                            ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
                            ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
                            tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
                            ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
                            tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
                            tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
                            ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
                            tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
                            cc(1,m1,k,1,1) = sn*(tr2+tr3)
                            cc(1,m1,k,1,3) = sn*(tr2-tr3)
                            cc(2,m1,k,1,1) = sn*(ti2+ti3)
                            cc(2,m1,k,1,3) = sn*(ti2-ti3)
                            cc(1,m1,k,1,2) = sn*(tr1+tr4)
                            cc(1,m1,k,1,4) = sn*(tr1-tr4)
                            cc(2,m1,k,1,2) = sn*(ti1+ti4)
                            cc(2,m1,k,1,4) = sn*(ti1-ti4)
                        end do
                    end do
                else
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
                            ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
                            tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
                            ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
                            tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
                            tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
                            ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
                            tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
                            ch(1,m2,k,1,1) = sn*(tr2+tr3)
                            ch(1,m2,k,3,1) = sn*(tr2-tr3)
                            ch(2,m2,k,1,1) = sn*(ti2+ti3)
                            ch(2,m2,k,3,1) = sn*(ti2-ti3)
                            ch(1,m2,k,2,1) = sn*(tr1+tr4)
                            ch(1,m2,k,4,1) = sn*(tr1-tr4)
                            ch(2,m2,k,2,1) = sn*(ti1+ti4)
                            ch(2,m2,k,4,1) = sn*(ti1-ti4)
                        end do
                    end do
                end if
            else
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
                        ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
                        tr4 = cc(2,m1,k,1,2)-cc(2,m1,k,1,4)
                        ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
                        tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
                        tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
                        ti4 = cc(1,m1,k,1,4)-cc(1,m1,k,1,2)
                        tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
                        ch(1,m2,k,1,1) = tr2+tr3
                        ch(1,m2,k,3,1) = tr2-tr3
                        ch(2,m2,k,1,1) = ti2+ti3
                        ch(2,m2,k,3,1) = ti2-ti3
                        ch(1,m2,k,2,1) = tr1+tr4
                        ch(1,m2,k,4,1) = tr1-tr4
                        ch(2,m2,k,2,1) = ti1+ti4
                        ch(2,m2,k,4,1) = ti1-ti4
                    end do
                end do
                do i=2,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ti1 = cc(2,m1,k,i,1)-cc(2,m1,k,i,3)
                            ti2 = cc(2,m1,k,i,1)+cc(2,m1,k,i,3)
                            ti3 = cc(2,m1,k,i,2)+cc(2,m1,k,i,4)
                            tr4 = cc(2,m1,k,i,2)-cc(2,m1,k,i,4)
                            tr1 = cc(1,m1,k,i,1)-cc(1,m1,k,i,3)
                            tr2 = cc(1,m1,k,i,1)+cc(1,m1,k,i,3)
                            ti4 = cc(1,m1,k,i,4)-cc(1,m1,k,i,2)
                            tr3 = cc(1,m1,k,i,2)+cc(1,m1,k,i,4)
                            ch(1,m2,k,1,i) = tr2+tr3
                            cr3 = tr2-tr3
                            ch(2,m2,k,1,i) = ti2+ti3
                            ci3 = ti2-ti3
                            cr2 = tr1+tr4
                            cr4 = tr1-tr4
                            ci2 = ti1+ti4
                            ci4 = ti1-ti4
                            ch(1,m2,k,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
                            ch(2,m2,k,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
                            ch(1,m2,k,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
                            ch(2,m2,k,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
                            ch(1,m2,k,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
                            ch(2,m2,k,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
                        end do
                    end do
                end do
            end if

        end subroutine cmf4kf


        subroutine cmf5kf(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(2,in1,l1,ido,5)
            real (wp) ch(2,in2,l1,5,ido)
            real (wp) chold1
            real (wp) chold2
            real (wp) ci2
            real (wp) ci3
            real (wp) ci4
            real (wp) ci5
            real (wp) cr2
            real (wp) cr3
            real (wp) cr4
            real (wp) cr5
            real (wp) di2
            real (wp) di3
            real (wp) di4
            real (wp) di5
            real (wp) dr2
            real (wp) dr3
            real (wp) dr4
            real (wp) dr5
            integer (ip) i
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp) sn
            real (wp) ti2
            real (wp) ti3
            real (wp) ti4
            real (wp) ti5
            real (wp) tr2
            real (wp) tr3
            real (wp) tr4
            real (wp) tr5

            real (wp), parameter :: SQRT5 = sqrt(5.0_wp)
            real (wp), parameter :: SQRT5_PLUS_5 = SQRT5 + 5.0_wp
            real (wp), parameter :: TI11 = -sqrt(SQRT5_PLUS_5/2)/2             !-0.9510565162951536_wp
            real (wp), parameter :: TI12 = -sqrt(5.0_wp/(2.0_wp*SQRT5_PLUS_5)) !-0.5877852522924731_wp
            real (wp), parameter :: TR11 =  (SQRT5 - 1.0_wp)/4                 ! 0.3090169943749474_wp
            real (wp), parameter :: TR12 = -(1.0_wp + SQRT5)/4                 !-0.8090169943749474_wp

            real (wp) wa(ido,4,2)

            m1d = (lot-1)*im1+1
            m2s = 1-im2

            if (1 >= ido) then
                sn = 1.0_wp/(5 * l1)
                if (na /= 1) then
                    do k=1,l1
                        do m1=1,m1d,im1
                            ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
                            ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
                            ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
                            ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
                            tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
                            tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
                            tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
                            tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
                            chold1 = sn*(cc(1,m1,k,1,1)+tr2+tr3)
                            chold2 = sn*(cc(2,m1,k,1,1)+ti2+ti3)
                            cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
                            ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
                            cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
                            ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
                            cc(1,m1,k,1,1) = chold1
                            cc(2,m1,k,1,1) = chold2
                            cr5 = ti11*tr5+ti12*tr4
                            ci5 = ti11*ti5+ti12*ti4
                            cr4 = ti12*tr5-ti11*tr4
                            ci4 = ti12*ti5-ti11*ti4
                            cc(1,m1,k,1,2) = sn*(cr2-ci5)
                            cc(1,m1,k,1,5) = sn*(cr2+ci5)
                            cc(2,m1,k,1,2) = sn*(ci2+cr5)
                            cc(2,m1,k,1,3) = sn*(ci3+cr4)
                            cc(1,m1,k,1,3) = sn*(cr3-ci4)
                            cc(1,m1,k,1,4) = sn*(cr3+ci4)
                            cc(2,m1,k,1,4) = sn*(ci3-cr4)
                            cc(2,m1,k,1,5) = sn*(ci2-cr5)
                        end do
                    end do
                else
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
                            ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
                            ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
                            ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
                            tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
                            tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
                            tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
                            tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
                            ch(1,m2,k,1,1) = sn*(cc(1,m1,k,1,1)+tr2+tr3)
                            ch(2,m2,k,1,1) = sn*(cc(2,m1,k,1,1)+ti2+ti3)
                            cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
                            ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
                            cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
                            ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
                            cr5 = ti11*tr5+ti12*tr4
                            ci5 = ti11*ti5+ti12*ti4
                            cr4 = ti12*tr5-ti11*tr4
                            ci4 = ti12*ti5-ti11*ti4
                            ch(1,m2,k,2,1) = sn*(cr2-ci5)
                            ch(1,m2,k,5,1) = sn*(cr2+ci5)
                            ch(2,m2,k,2,1) = sn*(ci2+cr5)
                            ch(2,m2,k,3,1) = sn*(ci3+cr4)
                            ch(1,m2,k,3,1) = sn*(cr3-ci4)
                            ch(1,m2,k,4,1) = sn*(cr3+ci4)
                            ch(2,m2,k,4,1) = sn*(ci3-cr4)
                            ch(2,m2,k,5,1) = sn*(ci2-cr5)
                        end do
                    end do
                end if
            else
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
                        ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
                        ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
                        ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
                        tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
                        tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
                        tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
                        tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
                        ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2+tr3
                        ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2+ti3
                        cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
                        ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
                        cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
                        ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
                        cr5 = ti11*tr5+ti12*tr4
                        ci5 = ti11*ti5+ti12*ti4
                        cr4 = ti12*tr5-ti11*tr4
                        ci4 = ti12*ti5-ti11*ti4
                        ch(1,m2,k,2,1) = cr2-ci5
                        ch(1,m2,k,5,1) = cr2+ci5
                        ch(2,m2,k,2,1) = ci2+cr5
                        ch(2,m2,k,3,1) = ci3+cr4
                        ch(1,m2,k,3,1) = cr3-ci4
                        ch(1,m2,k,4,1) = cr3+ci4
                        ch(2,m2,k,4,1) = ci3-cr4
                        ch(2,m2,k,5,1) = ci2-cr5
                    end do
                end do
                do i=2,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ti5 = cc(2,m1,k,i,2)-cc(2,m1,k,i,5)
                            ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,5)
                            ti4 = cc(2,m1,k,i,3)-cc(2,m1,k,i,4)
                            ti3 = cc(2,m1,k,i,3)+cc(2,m1,k,i,4)
                            tr5 = cc(1,m1,k,i,2)-cc(1,m1,k,i,5)
                            tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,5)
                            tr4 = cc(1,m1,k,i,3)-cc(1,m1,k,i,4)
                            tr3 = cc(1,m1,k,i,3)+cc(1,m1,k,i,4)
                            ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2+tr3
                            ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2+ti3
                            cr2 = cc(1,m1,k,i,1)+tr11*tr2+tr12*tr3
                            ci2 = cc(2,m1,k,i,1)+tr11*ti2+tr12*ti3
                            cr3 = cc(1,m1,k,i,1)+tr12*tr2+tr11*tr3
                            ci3 = cc(2,m1,k,i,1)+tr12*ti2+tr11*ti3
                            cr5 = ti11*tr5+ti12*tr4
                            ci5 = ti11*ti5+ti12*ti4
                            cr4 = ti12*tr5-ti11*tr4
                            ci4 = ti12*ti5-ti11*ti4
                            dr3 = cr3-ci4
                            dr4 = cr3+ci4
                            di3 = ci3+cr4
                            di4 = ci3-cr4
                            dr5 = cr2+ci5
                            dr2 = cr2-ci5
                            di5 = ci2-cr5
                            di2 = ci2+cr5
                            ch(1,m2,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
                            ch(2,m2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
                            ch(1,m2,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
                            ch(2,m2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
                            ch(1,m2,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
                            ch(2,m2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
                            ch(1,m2,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
                            ch(2,m2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
                        end do
                    end do
                end do
            end if

        end subroutine cmf5kf


        subroutine cmfgkf(lot, ido, iip, l1, lid, na, cc, cc1, im1, in1, &
            ch, ch1, im2, in2, wa)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) iip
            integer (ip) l1
            integer (ip) lid

            real (wp) cc(2,in1,l1,iip,ido)
            real (wp) cc1(2,in1,lid,iip)
            real (wp) ch(2,in2,l1,ido,iip)
            real (wp) ch1(2,in2,lid,iip)
            real (wp) chold1
            real (wp) chold2
            integer (ip) i
            integer (ip) idlj
            integer (ip) im1
            integer (ip) im2
            integer (ip) iipp2
            integer (ip) iipph
            integer (ip) j
            integer (ip) jc
            integer (ip) k
            integer (ip) ki
            integer (ip) l
            integer (ip) lc
            integer (ip) lot
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) na
            real (wp) sn
            real (wp) wa(ido,iip-1,2)
            real (wp) wai
            real (wp) war

            m1d = (lot-1)*im1+1
            m2s = 1-im2
            iipp2 = iip+2
            iipph = (iip+1)/2

            do  ki=1,lid
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
                    ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
                end do
            end do

            do j=2,iipph
                jc = iipp2-j
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch1(1,m2,ki,j) =  cc1(1,m1,ki,j)+cc1(1,m1,ki,jc)
                        ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)-cc1(1,m1,ki,jc)
                        ch1(2,m2,ki,j) =  cc1(2,m1,ki,j)+cc1(2,m1,ki,jc)
                        ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(2,m1,ki,jc)
                    end do
                end do
            end do

            do j=2,iipph
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        cc1(1,m1,ki,1) = cc1(1,m1,ki,1)+ch1(1,m2,ki,j)
                        cc1(2,m1,ki,1) = cc1(2,m1,ki,1)+ch1(2,m2,ki,j)
                    end do
                end do
            end do

            do l=2,iipph
                lc = iipp2-l
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        cc1(1,m1,ki,l) = ch1(1,m2,ki,1)+wa(1,l-1,1)*ch1(1,m2,ki,2)
                        cc1(1,m1,ki,lc) = -wa(1,l-1,2)*ch1(1,m2,ki,iip)
                        cc1(2,m1,ki,l) = ch1(2,m2,ki,1)+wa(1,l-1,1)*ch1(2,m2,ki,2)
                        cc1(2,m1,ki,lc) = -wa(1,l-1,2)*ch1(2,m2,ki,iip)
                    end do
                end do
                do j=3,iipph
                    jc = iipp2-j
                    idlj = mod((l-1)*(j-1),iip)
                    war = wa(1,idlj,1)
                    wai = -wa(1,idlj,2)
                    do ki=1,lid
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc1(1,m1,ki,l) = cc1(1,m1,ki,l)+war*ch1(1,m2,ki,j)
                            cc1(1,m1,ki,lc) = cc1(1,m1,ki,lc)+wai*ch1(1,m2,ki,jc)
                            cc1(2,m1,ki,l) = cc1(2,m1,ki,l)+war*ch1(2,m2,ki,j)
                            cc1(2,m1,ki,lc) = cc1(2,m1,ki,lc)+wai*ch1(2,m2,ki,jc)
                        end do
                    end do
                end do
            end do

            if (1 >= ido) then
                sn = 1.0_wp /(iip * l1)
                if (na /= 1) then
                    do ki=1,lid
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc1(1,m1,ki,1) = sn*cc1(1,m1,ki,1)
                            cc1(2,m1,ki,1) = sn*cc1(2,m1,ki,1)
                        end do
                    end do
                    do j=2,iipph
                        jc = iipp2-j
                        do ki=1,lid
                            do m1=1,m1d,im1
                                chold1 = sn*(cc1(1,m1,ki,j)-cc1(2,m1,ki,jc))
                                chold2 = sn*(cc1(1,m1,ki,j)+cc1(2,m1,ki,jc))
                                cc1(1,m1,ki,j) = chold1
                                cc1(2,m1,ki,jc) = sn*(cc1(2,m1,ki,j)-cc1(1,m1,ki,jc))
                                cc1(2,m1,ki,j) = sn*(cc1(2,m1,ki,j)+cc1(1,m1,ki,jc))
                                cc1(1,m1,ki,jc) = chold2
                            end do
                        end do
                    end do
                else
                    do ki=1,lid
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch1(1,m2,ki,1) = sn*cc1(1,m1,ki,1)
                            ch1(2,m2,ki,1) = sn*cc1(2,m1,ki,1)
                        end do
                    end do
                    do j=2,iipph
                        jc = iipp2-j
                        do ki=1,lid
                            m2 = m2s
                            do m1=1,m1d,im1
                                m2 = m2+im2
                                ch1(1,m2,ki,j) = sn*(cc1(1,m1,ki,j)-cc1(2,m1,ki,jc))
                                ch1(2,m2,ki,j) = sn*(cc1(2,m1,ki,j)+cc1(1,m1,ki,jc))
                                ch1(1,m2,ki,jc) = sn*(cc1(1,m1,ki,j)+cc1(2,m1,ki,jc))
                                ch1(2,m2,ki,jc) = sn*(cc1(2,m1,ki,j)-cc1(1,m1,ki,jc))
                            end do
                        end do
                    end do
                end if
            else
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
                        ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
                    end do
                end do
                do j=2,iipph
                    jc = iipp2-j
                    do ki=1,lid
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch1(1,m2,ki,j) = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
                            ch1(2,m2,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
                            ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
                            ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
                        end do
                    end do
                end do
                do i=1,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc(1,m1,k,1,i) = ch(1,m2,k,i,1)
                            cc(2,m1,k,1,i) = ch(2,m2,k,i,1)
                        end do
                    end do
                end do
                do j=2,iip
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc(1,m1,k,j,1) = ch(1,m2,k,1,j)
                            cc(2,m1,k,j,1) = ch(2,m2,k,1,j)
                        end do
                    end do
                end do
                do j=2,iip
                    do i=2,ido
                        do k=1,l1
                            m2 = m2s
                            do m1=1,m1d,im1
                                m2 = m2+im2
                                cc(1,m1,k,j,i) = wa(i,j-1,1)*ch(1,m2,k,i,j) &
                                    +wa(i,j-1,2)*ch(2,m2,k,i,j)
                                cc(2,m1,k,j,i) = wa(i,j-1,1)*ch(2,m2,k,i,j) &
                                    -wa(i,j-1,2)*ch(1,m2,k,i,j)
                            end do
                        end do
                    end do
                end do
            end if

        end subroutine cmfgkf

    end subroutine cmfm1f

    subroutine cosq1b(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! cosq1b: 64-bit float precision backward cosine quarter wave transform, 1d.
        !
        !  Purpose:
        !
        !  cosq1b computes the one-dimensional fourier transform of a sequence
        !  which is a cosine series with odd wave numbers.  this transform is
        !  referred to as the backward transform or fourier synthesis, transforming
        !  the sequence from spectral to physical space.
        !
        !  this transform is normalized since a call to cosq1b followed
        !  by a call to cosq1f (or vice-versa) reproduces the original
        !  array  within roundoff error.
        !
        !  integer n, the number of elements to be transformed
        !  in the sequence.  the transform is most efficient when n is a
        !  product of small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr); on input, containing the sequence
        !  to be transformed, and on output, containing the transformed sequence.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to cosq1i before the first call to routine cosq1f
        !  or cosq1b for a given transform length n.  wsave's contents may be
        !  re-used for subsequent calls to cosq1f and cosq1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least n.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  20, input error returned by lower level routine.
        !
        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) lenx
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)
        real (wp) x1

        !
        !==> Check validity of input arguments
        !
        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cosq1b', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n) ) then
            ier = 2
            call fft_error_handler('cosq1b', 8)
            return
        else if (lenwrk < n) then
            ier = 3
            call fft_error_handler('cosq1b', 10)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n < 2) then
            return
        else if (n == 2) then
            x1 = x(1,1)+x(1,2)
            x(1,2) = (x(1,1)-x(1,2))/sqrt(2.0_wp)
            x(1,1) = x1
        else
            call cosqb1(n,inc,x,wsave,work,local_error_flag)

            ! check error flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('cosq1b',-5)
            end if
        end if

    end subroutine cosq1b


    subroutine cosq1f(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! cosq1f: 64-bit float precision forward cosine quarter wave transform, 1d.
        !
        !  Purpose:
        !
        !  cosq1f computes the one-dimensional fourier transform of a sequence
        !  which is a cosine series with odd wave numbers.  this transform is
        !  referred to as the forward transform or fourier analysis, transforming
        !  the sequence from physical to spectral space.
        !
        !  this transform is normalized since a call to cosq1f followed
        !  by a call to cosq1b (or vice-versa) reproduces the original
        !  array  within roundoff error.
        !
        !  parameters:
        !
        !  integer n, the number of elements to be transformed
        !  in the sequence.  the transform is most efficient when n is a
        !  product of small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr); on input, containing the sequence
        !  to be transformed, and on output, containing the transformed sequence.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to cosq1i before the first call to routine cosq1f
        !  or cosq1b for a given transform length n.  wsave's contents may be
        !  re-used for subsequent calls to cosq1f and cosq1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least n.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  20, input error returned by lower level routine.
        !
        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) n
        integer (ip) lenx
        real (wp) tsqx
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Check validity of calling arguments
        !
        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cosq1f', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n) ) then
            ier = 2
            call fft_error_handler('cosq1f', 8)
            return
        else if (lenwrk < n) then
            ier = 3
            call fft_error_handler('cosq1f', 10)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n < 2) then
            return
        else if (n == 2) then
            tsqx = x(1,2)/sqrt(2.0_wp)
            x(1,2) = 0.5_wp *x(1,1)-tsqx
            x(1,1) = 0.5_wp *x(1,1)+tsqx
        else
            ! Peform cosine transform
            call cosqf1(n,inc,x,wsave,work,local_error_flag)

            ! Check error flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('cosq1f',-5)
            end if
        end if


    end subroutine cosq1f



    subroutine cosq1i(n, wsave, lensav, ier)
        !
        ! cosq1i: initialization for cosq1b and cosq1f.
        !
        !  Purpose:
        !
        !  cosq1i initializes array wsave for use in its companion routines
        !  cosq1f and cosq1b.  the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  separate wsave arrays are required for different
        !  values of n.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n is a product
        !  of small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  real wsave(lensav), containing the prime factors of n
        !  and also containing certain trigonometric values which will be used
        !  in routines cosq1b or cosq1f.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  20, input error returned by lower level routine.
        !
        integer (ip) lensav
        real (wp) dt
        real (wp) fk
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) n
        real (wp), parameter :: HALF_PI = acos(-1.0_wp)/2
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lensav < get_1d_saved_workspace_size(n) ) then
            ier = 2
            call fft_error_handler('cosq1i', 3)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        dt = HALF_PI /n
        fk = 0.0_wp

        do k=1,n
            fk = fk + 1.0_wp
            wsave(k) = cos(fk*dt)
        end do

        associate( lnsv => n+int(log(real(n, kind=wp))/log(2.0_wp))+4 )

            call rfft1i(n, wsave(n+1), lnsv, local_error_flag)

        end associate

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cosq1i',-5)
        end if

    end subroutine cosq1i


    subroutine cosqb1(n, inc, x, wsave, work, ier)

        integer (ip) inc

        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) modn
        integer (ip) n
        integer (ip) np2
        integer (ip) ns2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) xim1

        ier = 0
        ns2 = (n+1)/2
        np2 = n+2

        do i=3,n,2
            xim1 = x(1,i-1)+x(1,i)
            x(1,i) = 0.5_wp * (x(1,i-1)-x(1,i))
            x(1,i-1) = 0.5_wp * xim1
        end do

        x(1,1) = 0.5_wp * x(1,1)
        modn = mod(n,2)

        if (modn == 0) then
            x(1,n) = 0.5_wp * x(1,n)
        end if

        associate( &
            lenx => inc*(n-1)  + 1, &
            lnsv => n + int(log(real(n, kind=wp) )/log(2.0_wp)) + 4, &
            lnwk => n &
            )

            call rfft1b(n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

        end associate

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cosqb1',-5)
            return
        end if

        do k=2,ns2
            kc = np2-k
            work(k) = wsave(k-1)*x(1,kc)+wsave(kc-1)*x(1,k)
            work(kc) = wsave(k-1)*x(1,k)-wsave(kc-1)*x(1,kc)
        end do

        if (modn == 0) then
            x(1,ns2+1) = 2.0_wp * wsave(ns2) * x(1,ns2+1)
        end if

        do k=2,ns2
            kc = np2-k
            x(1,k) = work(k)+work(kc)
            x(1,kc) = work(k)-work(kc)
        end do

        x(1,1) = 2.0_wp * x(1,1)

    end subroutine cosqb1


    subroutine cosqf1(n, inc, x, wsave, work, ier)

        integer (ip) inc

        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) modn
        integer (ip) n
        integer (ip) np2
        integer (ip) ns2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) xim1

        ier = 0
        ns2 = (n+1)/2
        np2 = n+2

        do k=2,ns2
            kc = np2-k
            work(k)  = x(1,k)+x(1,kc)
            work(kc) = x(1,k)-x(1,kc)
        end do

        modn = mod(n,2)

        if (modn == 0) then
            work(ns2+1) = 2.0_wp * x(1,ns2+1)
        end if

        do k=2,ns2
            kc = np2-k
            x(1,k)  = wsave(k-1)*work(kc)+wsave(kc-1)*work(k)
            x(1,kc) = wsave(k-1)*work(k) -wsave(kc-1)*work(kc)
        end do

        if (modn == 0) then
            x(1,ns2+1) = wsave(ns2)*work(ns2+1)
        end if

        associate( &
            lenx => inc*(n-1)  + 1, &
            lnsv => n + int(log(real(n, kind=wp) )/log(2.0_wp)) + 4, &
            lnwk => n &
            )

            call rfft1f(n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

        end associate

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cosqf1',-5)
            return
        end if

        do i=3,n,2
            xim1 = 0.5_wp * (x(1,i-1)+x(1,i))
            x(1,i) = 0.5_wp * (x(1,i-1)-x(1,i))
            x(1,i-1) = xim1
        end do

    end subroutine cosqf1



    subroutine cosqmb(lot, jump, n, inc, x, lenx, wsave, lensav, work, lenwrk, &
        ier)
        ! cosqmb: 64-bit float precision backward cosine quarter wave, multiple vectors.
        !
        !  Purpose:
        !
        !  cosqmb computes the one-dimensional fourier transform of multiple
        !  sequences, each of which is a cosine series with odd wave numbers.
        !  this transform is referred to as the backward transform or fourier
        !  synthesis, transforming the sequences from spectral to physical space.
        !
        !  this transform is normalized since a call to cosqmb followed
        !  by a call to cosqmf (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  integer lot, the number of sequences to be transformed
        !  within array r.
        !
        !  integer jump, the increment between the locations,
        !  in array r, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the same sequence.
        !
        !  input/real r(lenr), array containing lot sequences,
        !  each having length n.  r can have any number of dimensions, but the total
        !  number of locations must be at least lenr.  on input, r contains the data
        !  to be transformed, and on output, the transformed data.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least (lot-1)*jump + inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to cosqmi before the first call to routine cosqmf
        !  or cosqmb for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to cosqmf and cosqmb with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least lot*n.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  4, input parameters inc,jump,n,lot are not consistent;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) lenx
        integer (ip) lj
        integer (ip) lot
        integer (ip) m
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)
        real (wp) x1

        !
        !==> Check validity of calling arguments
        !
        if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cosqmb', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('cosqmb', 8)
            return
        else if (lenwrk < lot*n) then
            ier = 3
            call fft_error_handler('cosqmb', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('cosqmb', -1)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        lj = (lot-1)*jump+1
        if (n < 2) then
            do m=1,lj,jump
                x(m,1) = x(m,1)
            end do
        else if (n == 2) then
            do m=1,lj,jump
                x1 = x(m,1)+x(m,2)
                x(m,2) = (x(m,1)-x(m,2))/sqrt(2.0_wp)
                x(m,1) = x1
            end do
        else
            call mcsqb1(lot,jump,n,inc,x,wsave,work,local_error_flag)

            ! Check error flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('cosqmb',-5)
            end if
        end if

    end subroutine cosqmb



    subroutine cosqmf(lot, jump, n, inc, x, lenx, wsave, lensav, work, &
        lenwrk, ier)
        ! cosqmf: 64-bit float precision forward cosine quarter wave, multiple vectors.
        !
        !  purpose:
        !
        !  cosqmf computes the one-dimensional fourier transform of multiple
        !  sequences within a real array, where each of the sequences is a
        !  cosine series with odd wave numbers.  this transform is referred to
        !  as the forward transform or fourier synthesis, transforming the
        !  sequences from spectral to physical space.
        !
        !  this transform is normalized since a call to cosqmf followed
        !  by a call to cosqmb (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  parameters:
        !
        !  integer lot, the number of sequences to be transformed
        !  within array r.
        !
        !  integer jump, the increment between the locations, in
        !  array r, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the same sequence.
        !
        !  input/real r(lenr), array containing lot sequences,
        !  each having length n.  r can have any number of dimensions, but the total
        !  number of locations must be at least lenr.  on input, r contains the data
        !  to be transformed, and on output, the transformed data.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least (lot-1)*jump + inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to cosqmi before the first call to routine cosqmf
        !  or cosqmb for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to cosqmf and cosqmb with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least lot*n.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  4, input parameters inc,jump,n,lot are not consistent;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) lenx
        integer (ip) lj
        integer (ip) lot
        integer (ip) m
        integer (ip) n
        real (wp) tsqx
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Validity of input arguments
        !
        if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cosqmf', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('cosqmf', 8)
            return
        else if (lenwrk < lot*n) then
            ier = 3
            call fft_error_handler('cosqmf', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('cosqmf', -1)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        lj = (lot-1)*jump+1

        if (n < 2) then
            return
        else if (n == 2) then
            do m=1,lj,jump
                tsqx = x(m,2)/sqrt(2.0_wp)
                x(m,2) = 0.5_wp * x(m,1)-tsqx
                x(m,1) = 0.5_wp * x(m,1)+tsqx
            end do
        else
            call mcsqf1(lot,jump,n,inc,x,wsave,work,local_error_flag)

            ! Check error flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('cosqmf',-5)
            end if
        end if

    end subroutine cosqmf


    subroutine cosqmi(n, wsave, lensav, ier)
        !
        ! cosqmi: initialization for cosqmb and cosqmf.
        !
        !  purpose:
        !
        !  cosqmi initializes array wsave for use in its companion routines
        !  cosqmf and cosqmb.  the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  separate wsave arrays are required for different
        !  values of n.
        !
        !  parameters:
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  real wsave(lensav), containing the prime factors of
        !  n and also containing certain trigonometric values which will be used
        !  in routines cosqmb or cosqmf.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  20, input error returned by lower level routine.
        !
        integer (ip) lensav
        real (wp) dt
        real (wp) fk
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) lnsv
        integer (ip) n
        real (wp), parameter ::  HALF_PI = acos(-1.0_wp)/2
        real (wp) wsave(lensav)

        if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('cosqmi', 3)
            return
        else
            ier = 0
        end if

        dt = HALF_PI/n
        fk = 0.0_wp

        do k=1,n
            fk = fk + 1.0_wp
            wsave(k) = cos(fk*dt)
        end do

        lnsv = get_1d_saved_workspace_size(n) - n

        call rfftmi(n, wsave(n+1), lnsv, local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('cosqmi',-5)
        end if

    end subroutine cosqmi


    subroutine cost1b(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! cost1b: 64-bit float precision backward cosine transform, 1d.
        !
        !  purpose:
        !
        !  cost1b computes the one-dimensional fourier transform of an even
        !  sequence within a real array.  this transform is referred to as
        !  the backward transform or fourier synthesis, transforming the sequence
        !  from spectral to physical space.
        !
        !  this transform is normalized since a call to cost1b followed
        !  by a call to cost1f (or vice-versa) reproduces the original array
        !  within roundoff error.
        !
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n-1 is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr), containing the sequence to
        !   be transformed.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to cost1i before the first call to routine cost1f
        !  or cost1b for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to cost1f and cost1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least n-1.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) lenx
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Check validity of input arguments
        !
        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cost1b', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('cost1b', 8)
            return
        else if (lenwrk < n-1) then
            ier = 3
            call fft_error_handler('cost1b', 10)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then

            call costb1(n,inc,x,wsave,work,local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('cost1b',-5)
            end if
        end if

    end subroutine cost1b


    subroutine cost1f(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! cost1f: 64-bit float precision forward cosine transform, 1d.
        !
        !  purpose:
        !
        !  cost1f computes the one-dimensional fourier transform of an even
        !  sequence within a real array.  this transform is referred to as the
        !  forward transform or fourier analysis, transforming the sequence
        !  from  physical to spectral space.
        !
        !  this transform is normalized since a call to cost1f followed by a call
        !  to cost1b (or vice-versa) reproduces the original array within
        !  roundoff error.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n-1 is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr), containing the sequence to
        !  be transformed.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to cost1i before the first call to routine cost1f
        !  or cost1b for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to cost1f and cost1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least n-1.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) lenx
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Check validity of input arguments
        !
        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cost1f', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('cost1f', 8)
            return
        else if (lenwrk < n-1) then
            ier = 3
            call fft_error_handler('cost1f', 10)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then

            call costf1(n,inc,x,wsave,work,local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('cost1f',-5)
            end if
        end if

    end subroutine cost1f



    subroutine cost1i(n, wsave, lensav, ier)
        !
        ! cost1i: initialization for cost1b and cost1f.
        !
        !  purpose:
        !
        !  cost1i initializes array wsave for use in its companion routines
        !  cost1f and cost1b.  the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  separate wsave arrays are required for different
        !  values of n.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n-1 is a product
        !  of small primes.
        !
        !  integer lensav, dimension of wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  real wsave(lensav), containing the prime factors of
        !  n and also containing certain trigonometric values which will be used in
        !  routines cost1b or cost1f.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  20, input error returned by lower level routine.
        !


        integer (ip) lensav

        real (wp) dt
        real (wp) fk
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) lnsv
        integer (ip) n
        integer (ip) nm1
        integer (ip) np1
        integer (ip) ns2
        real (wp), parameter :: PI = acos(-1.0_wp)
        real (wp) wsave(lensav)


        !
        !==> Check validity of input arguments
        !
        if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('cost1i', 3)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n > 3) then

            nm1 = n-1
            np1 = n+1
            ns2 = n/2
            dt = pi/ nm1
            fk = 0.0_wp

            do k=2,ns2
                kc = np1-k
                fk = fk + 1.0_wp
                wsave(k) = 2.0_wp * sin(fk*dt)
                wsave(kc) = 2.0_wp * cos(fk*dt)
            end do

            lnsv = nm1 + int(log(real(nm1, kind=wp) )/log(2.0_wp)) +4

            call rfft1i(nm1, wsave(n+1), lnsv, local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('cost1i',-5)
            end if
        end if

    end subroutine cost1i



    subroutine costb1(n, inc, x, wsave, work, ier)

        integer (ip) inc
        real (wp) dsum
        real (wp) fnm1s2
        real (wp) fnm1s4
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) modn
        integer (ip) n
        integer (ip) nm1
        integer (ip) np1
        integer (ip) ns2
        real (wp) t1
        real (wp) t2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) x1h
        real (wp) x1p3
        real (wp) x2
        real (wp) xi

        ier = 0

        nm1 = n-1
        np1 = n+1
        ns2 = n/2

        if (n < 2) then
            return
        else if (n == 2) then
            x1h = x(1,1)+x(1,2)
            x(1,2) = x(1,1)-x(1,2)
            x(1,1) = x1h
        else
            if (3 >= n) then
                x1p3 = x(1,1)+x(1,3)
                x2 = x(1,2)
                x(1,2) = x(1,1)-x(1,3)
                x(1,1) = x1p3+x2
                x(1,3) = x1p3-x2
            else
                x(1,1) = x(1,1)+x(1,1)
                x(1,n) = x(1,n)+x(1,n)
                dsum = x(1,1)-x(1,n)
                x(1,1) = x(1,1)+x(1,n)

                do k=2,ns2
                    kc = np1-k
                    t1 = x(1,k)+x(1,kc)
                    t2 = x(1,k)-x(1,kc)
                    dsum = dsum+wsave(kc)*t2
                    t2 = wsave(k)*t2
                    x(1,k) = t1-t2
                    x(1,kc) = t1+t2
                end do

                modn = mod(n,2)

                if (modn /= 0) then
                    x(1,ns2+1) = 2.0_wp * x(1,ns2+1)
                end if

                lenx = inc*(nm1-1) + 1
                lnsv = nm1 + int(log(real(nm1, kind=wp))/log(2.0_wp), kind=ip) + 4
                lnwk = nm1

                call rfft1f(nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

                if (local_error_flag /= 0) then
                    ier = 20
                    call fft_error_handler('costb1',-5)
                    return
                end if

                fnm1s2 = real(nm1, kind=wp)/2
                dsum = 0.5_wp * dsum
                x(1,1) = fnm1s2*x(1,1)

                if (mod(nm1,2) == 0) then
                    x(1,nm1) = 2.0_wp * x(1,nm1)
                end if

                fnm1s4 = real(nm1, kind=wp)/4

                do i=3,n,2
                    xi = fnm1s4*x(1,i)
                    x(1,i) = fnm1s4*x(1,i-1)
                    x(1,i-1) = dsum
                    dsum = dsum+xi
                end do

                if (modn == 0) then
                    x(1,n) = dsum
                end if
            end if
        end if

    end subroutine costb1

    subroutine costf1(n, inc, x, wsave, work, ier)

        integer (ip) inc
        real (wp) dsum
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) modn
        integer (ip) n
        integer (ip) nm1
        integer (ip) np1
        integer (ip) ns2
        real (wp) snm1
        real (wp) t1
        real (wp) t2
        real (wp) tx2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) x1h
        real (wp) x1p3
        real (wp) xi

        ier = 0

        nm1 = n-1
        np1 = n+1
        ns2 = n/2

        if (n < 2) then
            return
        else if (n == 2) then
            x1h = x(1,1)+x(1,2)
            x(1,2) = 0.5_wp * (x(1,1)-x(1,2))
            x(1,1) = 0.5_wp * x1h
        else
            if (3 >= n) then
                x1p3 = x(1,1)+x(1,3)
                tx2 = x(1,2)+x(1,2)
                x(1,2) = 0.5_wp * (x(1,1)-x(1,3))
                x(1,1) = 0.25_wp *(x1p3+tx2)
                x(1,3) = 0.25_wp *(x1p3-tx2)
            else
                dsum = x(1,1)-x(1,n)
                x(1,1) = x(1,1)+x(1,n)
                do k=2,ns2
                    kc = np1-k
                    t1 = x(1,k)+x(1,kc)
                    t2 = x(1,k)-x(1,kc)
                    dsum = dsum+wsave(kc)*t2
                    t2 = wsave(k)*t2
                    x(1,k) = t1-t2
                    x(1,kc) = t1+t2
                end do

                modn = mod(n,2)

                if (modn /= 0) then
                    x(1,ns2+1) = x(1,ns2+1)+x(1,ns2+1)
                end if

                lenx = inc*(nm1-1)  + 1
                lnsv = nm1 + int(log(real(nm1, kind=wp) )/log(2.0_wp)) + 4
                lnwk = nm1

                call rfft1f(nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

                if (local_error_flag /= 0) then
                    ier = 20
                    call fft_error_handler('costf1',-5)
                    return
                end if

                snm1 = 1.0_wp /nm1
                dsum = snm1*dsum

                if (mod(nm1,2) == 0) then
                    x(1,nm1) = x(1,nm1)+x(1,nm1)
                end if

                do i=3,n,2
                    xi = 0.5_wp * x(1,i)
                    x(1,i) = 0.5_wp * x(1,i-1)
                    x(1,i-1) = dsum
                    dsum = dsum+xi
                end do

                if (modn == 0) then
                    x(1,n) = dsum
                end if

                x(1,1) = 0.5_wp * x(1,1)
                x(1,n) = 0.5_wp * x(1,n)
            end if
        end if

    end subroutine costf1

    subroutine costmb(lot, jump, n, inc, x, lenx, wsave, lensav, work, &
        lenwrk, ier)
        !
        ! costmb: 64-bit float precision backward cosine transform, multiple vectors.
        !
        !  purpose:
        !
        !  costmb computes the one-dimensional fourier transform of multiple
        !  even sequences within a real array.  this transform is referred to
        !  as the backward transform or fourier synthesis, transforming the
        !  sequences from spectral to physical space.
        !
        !  this transform is normalized since a call to costmb followed
        !  by a call to costmf (or vice-versa) reproduces the original
        !  array  within roundoff error.
        !
        !  parameters:
        !
        !  integer lot, the number of sequences to be transformed
        !  within array r.
        !
        !  integer jump, the increment between the locations, in
        !  array r, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n-1 is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array r, of two consecutive elements within the same sequence.
        !
        !  input/real r(lenr), array containing lot sequences,
        !  each having length n.  on input, the data to be transformed; on output,
        !  the transormed data.  r can have any number of dimensions, but the total
        !  number of locations must be at least lenr.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least (lot-1)*jump + inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to costmi before the first call to routine costmf
        !  or costmb for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to costmf and costmb with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least lot*(n+1).
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  4, input parameters inc,jump,n,lot are not consistent;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) iw1
        integer (ip) jump
        integer (ip) lenx
        integer (ip) lot
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Check validity of input arguments
        !
        if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('costmb', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('costmb', 8)
            return
        else if (lenwrk < lot*(n+1)) then
            ier = 3
            call fft_error_handler('costmb', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('costmb', -1)
            return
        else
            ier = 0
        end if

        iw1 = 2*lot+1
        call mcstb1(lot,jump,n,inc,x,wsave,work,work(iw1),local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('costmb',-5)
        end if

    end subroutine costmb



    subroutine costmf(lot, jump, n, inc, x, lenx, wsave, lensav, work, &
        lenwrk, ier)


        !
        !! COSTMF: 64-bit float precision forward cosine transform, multiple vectors.
        !
        !  Purpose:
        !
        !  COSTMF computes the one-dimensional Fourier transform of multiple
        !  even sequences within a real array.  This transform is referred to
        !  as the forward transform or Fourier analysis, transforming the
        !  sequences from physical to spectral space.
        !
        !  This transform is normalized since a call to COSTMF followed
        !  by a call to COSTMB (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  Parameters:
        !
        !  integer LOT, the number of sequences to be transformed
        !  within array R.
        !
        !  integer JUMP, the increment between the locations,
        !  in array R, of the first elements of two consecutive sequences to
        !  be transformed.
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N-1 is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations,
        !  in array R, of two consecutive elements within the same sequence.
        !
        !  Input/real R(LENR), array containing LOT sequences,
        !  each having length N.  On input, the data to be transformed; on output,
        !  the transormed data.  R can have any number of dimensions, but the total
        !  number of locations must be at least LENR.
        !
        !  integer LENR, the dimension of the  R array.
        !  LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to COSTMI before the first call to routine COSTMF
        !  or COSTMB for a given transform length N.  wsave's contents may be re-used
        !  for subsequent calls to COSTMF and COSTMB with the same N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least LOT*(N+1).
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  4, input parameters INC,JUMP,N,LOT are not consistent;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) iw1
        integer (ip) jump
        integer (ip) lenx
        integer (ip) lot
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Check validity of input arguments
        !
        if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('costmf', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('costmf', 8)
            return
        else if (lenwrk < lot*(n+1)) then
            ier = 3
            call fft_error_handler('costmf', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('costmf', -1)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        iw1 = 2*lot+1
        call mcstf1(lot,jump,n,inc,x,wsave,work,work(iw1),local_error_flag)

        ! Check error flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('costmf',-5)
        end if

    end subroutine costmf



    subroutine costmi(n, wsave, lensav, ier)
        !
        ! costmi: initialization for costmb and costmf.
        !
        !  purpose:
        !
        !  costmi initializes array wsave for use in its companion routines
        !  costmf and costmb.  the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  separate wsave arrays are required for different
        !  values of n.
        !
        !  parameters:
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4
        !
        !  real wsave(lensav), containing the prime factors of n
        !  and also containing certain trigonometric values which will be used
        !  in routines costmb or costmf.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  20, input error returned by lower level routine.
        !
        integer (ip) lensav
        real (wp) dt
        real (wp) fk
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) lnsv
        integer (ip) n
        integer (ip) nm1
        integer (ip) np1
        integer (ip) ns2
        real (wp), parameter :: PI = acos(-1.0_wp)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('costmi', 3)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n > 3) then

            nm1 = n-1
            np1 = n+1
            ns2 = n/2
            dt = PI/nm1
            fk = 0.0_wp

            do k=2,ns2
                kc = np1-k
                fk = fk + 1.0_wp
                wsave(k) = 2.0_wp * sin(fk*dt)
                wsave(kc) = 2.0_wp * cos(fk*dt)
            end do

            lnsv = nm1 + int(log(real(nm1, kind=wp) )/log(2.0_wp)) + 4

            call rfftmi(nm1, wsave(n+1), lnsv, local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('costmi',-5)
            end if
        end if

    end subroutine costmi


    subroutine mcsqb1(lot,jump,n,inc,x,wsave,work,ier)
        integer (ip) inc
        integer (ip) lot
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) lj
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) m
        integer (ip) m1
        integer (ip) modn
        integer (ip) n
        integer (ip) np2
        integer (ip) ns2
        real (wp) work(lot,*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) xim1

        ier = 0
        lj = (lot-1)*jump+1
        ns2 = (n+1)/2
        np2 = n+2

        do i=3,n,2
            do m=1,lj,jump
                xim1 = x(m,i-1)+x(m,i)
                x(m,i) = 0.5_wp * (x(m,i-1)-x(m,i))
                x(m,i-1) = 0.5_wp * xim1
            end do
        end do

        do m=1,lj,jump
            x(m,1) = 0.5_wp * x(m,1)
        end do

        modn = mod(n,2)
        if (modn == 0) then
            do m=1,lj,jump
                x(m,n) = 0.5_wp * x(m,n)
            end do
        end if

        lenx = (lot-1)*jump + inc*(n-1)  + 1
        lnsv = n + int(log(real(n, kind=wp) )/log(2.0_wp)) + 4
        lnwk = lot*n

        call rfftmb(lot,jump,n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('mcsqb1',-5)
            return
        end if

        do k=2,ns2
            kc = np2-k
            m1 = 0
            do m=1,lj,jump
                m1 = m1 + 1
                work(m1,k) = wsave(k-1)*x(m,kc)+wsave(kc-1)*x(m,k)
                work(m1,kc) = wsave(k-1)*x(m,k)-wsave(kc-1)*x(m,kc)
            end do
        end do

        if (modn == 0) then
            do m=1,lj,jump
                x(m,ns2+1) = wsave(ns2)*(x(m,ns2+1)+x(m,ns2+1))
            end do
        end if

        do k=2,ns2
            kc = np2-k
            m1 = 0
            do m=1,lj,jump
                m1 = m1 + 1
                x(m,k) = work(m1,k)+work(m1,kc)
                x(m,kc) = work(m1,k)-work(m1,kc)
            end do
        end do

        do m=1,lj,jump
            x(m,1) = x(m,1)+x(m,1)
        end do

    end subroutine mcsqb1



    subroutine mcsqf1(lot,jump,n,inc,x,wsave,work,ier)

        integer (ip) inc
        integer (ip) lot
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) lj
        integer (ip) m
        integer (ip) m1
        integer (ip) modn
        integer (ip) n
        integer (ip) np2
        integer (ip) ns2
        real (wp) work(lot,*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) xim1

        ier = 0

        lj = (lot-1)*jump+1
        ns2 = (n+1)/2
        np2 = n+2

        do k=2,ns2
            kc = np2-k
            m1 = 0
            do m=1,lj,jump
                m1 = m1 + 1
                work(m1,k) = x(m,k)+x(m,kc)
                work(m1,kc) = x(m,k)-x(m,kc)
            end do
        end do

        modn = mod(n,2)

        if (modn == 0) then
            m1 = 0
            do m=1,lj,jump
                m1 = m1 + 1
                work(m1,ns2+1) = x(m,ns2+1)+x(m,ns2+1)
            end do
        end if

        do k=2,ns2
            kc = np2-k
            m1 = 0
            do m=1,lj,jump
                m1 = m1 + 1
                x(m,k)  = wsave(k-1)*work(m1,kc)+wsave(kc-1)*work(m1,k)
                x(m,kc) = wsave(k-1)*work(m1,k) -wsave(kc-1)*work(m1,kc)
            end do
        end do

        if (modn == 0) then
            m1 = 0
            do m=1,lj,jump
                m1 = m1 + 1
                x(m,ns2+1) = wsave(ns2)*work(m1,ns2+1)
            end do
        end if

        lenx = (lot-1)*jump + inc*(n-1)  + 1
        lnsv = n + int(log(real(n, kind=wp) )/log(2.0_wp)) + 4
        lnwk = lot*n

        call rfftmf(lot,jump,n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('mcsqf1',-5)
            return
        end if

        do i=3,n,2
            do m=1,lj,jump
                xim1 = 0.5_wp * (x(m,i-1)+x(m,i))
                x(m,i) = 0.5_wp * (x(m,i-1)-x(m,i))
                x(m,i-1) = xim1
            end do
        end do

    end subroutine mcsqf1



    subroutine mcstb1(lot,jump,n,inc,x,wsave,dsum,work,ier)

        integer (ip) inc
        real (wp) dsum(*)
        real (wp) fnm1s2
        real (wp) fnm1s4
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) lj
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) lot
        integer (ip) m
        integer (ip) m1
        integer (ip) modn
        integer (ip) n
        integer (ip) nm1
        integer (ip) np1
        integer (ip) ns2
        real (wp) t1
        real (wp) t2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) x1h
        real (wp) x1p3
        real (wp) x2
        real (wp) xi

        ier = 0

        nm1 = n-1
        np1 = n+1
        ns2 = n/2
        lj = (lot-1)*jump+1

        if (n < 2) then
            return
        else if (n == 2) then
            do m=1,lj,jump
                x1h = x(m,1)+x(m,2)
                x(m,2) = x(m,1)-x(m,2)
                x(m,1) = x1h
            end do
        else
            if (3 >= n) then
                do m=1,lj,jump
                    x1p3 = x(m,1)+x(m,3)
                    x2 = x(m,2)
                    x(m,2) = x(m,1)-x(m,3)
                    x(m,1) = x1p3+x2
                    x(m,3) = x1p3-x2
                end do
            else
                do m=1,lj,jump
                    x(m,1) = x(m,1)+x(m,1)
                    x(m,n) = x(m,n)+x(m,n)
                end do

                m1 = 0

                do m=1,lj,jump
                    m1 = m1+1
                    dsum(m1) = x(m,1)-x(m,n)
                    x(m,1) = x(m,1)+x(m,n)
                end do

                do k=2,ns2
                    m1 = 0
                    do m=1,lj,jump
                        m1 = m1+1
                        kc = np1-k
                        t1 = x(m,k)+x(m,kc)
                        t2 = x(m,k)-x(m,kc)
                        dsum(m1) = dsum(m1)+wsave(kc)*t2
                        t2 = wsave(k)*t2
                        x(m,k) = t1-t2
                        x(m,kc) = t1+t2
                    end do
                end do

                modn = mod(n,2)

                if (modn /= 0) then
                    do m=1,lj,jump
                        x(m,ns2+1) = x(m,ns2+1)+x(m,ns2+1)
                    end do
                end if

                lenx = (lot-1)*jump + inc*(nm1-1)  + 1
                lnsv = nm1 + int(log(real(nm1, kind=wp))/log(2.0_wp)) + 4
                lnwk = lot*nm1

                call rfftmf(lot,jump,nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

                if (local_error_flag /= 0) then
                    ier = 20
                    call fft_error_handler('mcstb1',-5)
                    return
                end if

                fnm1s2 = real(nm1, kind=wp)/2
                m1 = 0

                do m=1,lj,jump
                    m1 = m1+1
                    dsum(m1) = 0.5_wp * dsum(m1)
                    x(m,1) = fnm1s2 * x(m,1)
                end do

                if (mod(nm1,2) == 0) then
                    do m=1,lj,jump
                        x(m,nm1) = x(m,nm1)+x(m,nm1)
                    end do
                end if

                fnm1s4 = real(nm1, kind=wp)/4

                do i=3,n,2
                    m1 = 0
                    do m=1,lj,jump
                        m1 = m1+1
                        xi = fnm1s4*x(m,i)
                        x(m,i) = fnm1s4*x(m,i-1)
                        x(m,i-1) = dsum(m1)
                        dsum(m1) = dsum(m1)+xi
                    end do
                end do
                if (modn == 0) then
                    m1 = 0
                    do m=1,lj,jump
                        m1 = m1+1
                        x(m,n) = dsum(m1)
                    end do
                end if
            end if
        end if

    end subroutine mcstb1

    subroutine mcstf1(lot,jump,n,inc,x,wsave,dsum,work,ier)

        integer (ip) inc
        real (wp) dsum(*)
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) lj
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) lot
        integer (ip) m
        integer (ip) m1
        integer (ip) modn
        integer (ip) n
        integer (ip) nm1
        integer (ip) np1
        integer (ip) ns2
        real (wp) snm1
        real (wp) t1
        real (wp) t2
        real (wp) tx2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) x1h
        real (wp) x1p3
        real (wp) xi

        ier = 0

        nm1 = n-1
        np1 = n+1
        ns2 = n/2
        lj = (lot-1)*jump+1

        if (n < 2) then
            return
        else if (n == 2) then
            do m=1,lj,jump
                x1h = x(m,1)+x(m,2)
                x(m,2) = 0.5_wp * (x(m,1)-x(m,2))
                x(m,1) = 0.5_wp * x1h
            end do
        else
            if (3 >= n) then
                do m=1,lj,jump
                    x1p3 = x(m,1)+x(m,3)
                    tx2 = x(m,2)+x(m,2)
                    x(m,2) = 0.5_wp * (x(m,1)-x(m,3))
                    x(m,1) = 0.25_wp * (x1p3+tx2)
                    x(m,3) = 0.25_wp * (x1p3-tx2)
                end do
            else
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    dsum(m1) = x(m,1)-x(m,n)
                    x(m,1) = x(m,1)+x(m,n)
                end do
                do k=2,ns2
                    m1 = 0
                    do m=1,lj,jump
                        m1 = m1+1
                        kc = np1-k
                        t1 = x(m,k)+x(m,kc)
                        t2 = x(m,k)-x(m,kc)
                        dsum(m1) = dsum(m1)+wsave(kc)*t2
                        t2 = wsave(k)*t2
                        x(m,k) = t1-t2
                        x(m,kc) = t1+t2
                    end do
                end do

                modn = mod(n,2)

                if (modn /= 0) then
                    do m=1,lj,jump
                        x(m,ns2+1) = x(m,ns2+1)+x(m,ns2+1)
                    end do
                end if

                lenx = (lot-1)*jump + inc*(nm1-1)  + 1
                lnsv = nm1 + int(log(real(nm1, kind=wp))/log(2.0_wp)) + 4
                lnwk = lot*nm1

                call rfftmf(lot,jump,nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

                if (local_error_flag /= 0) then
                    ier = 20
                    call fft_error_handler('mcstf1',-5)
                    return
                end if

                snm1 = 1.0_wp/nm1

                do m=1,lot
                    dsum(m) = snm1*dsum(m)
                end do

                if (mod(nm1,2) == 0) then
                    do m=1,lj,jump
                        x(m,nm1) = x(m,nm1)+x(m,nm1)
                    end do
                end if

                do i=3,n,2
                    m1 = 0
                    do m=1,lj,jump
                        m1 = m1+1
                        xi = 0.5_wp * x(m,i)
                        x(m,i) = 0.5_wp * x(m,i-1)
                        x(m,i-1) = dsum(m1)
                        dsum(m1) = dsum(m1)+xi
                    end do
                end do

                if (modn == 0) then
                    m1 = 0
                    do m=1,lj,jump
                        m1 = m1+1
                        x(m,n) = dsum(m1)
                    end do
                end if

                do m=1,lj,jump
                    x(m,1) = 0.5_wp * x(m,1)
                    x(m,n) = 0.5_wp * x(m,n)
                end do
            end if
        end if

    end subroutine mcstf1


    subroutine mrftb1(m,im,n,in,c,ch,wa,fac)

        integer (ip) in
        integer (ip) m
        integer (ip) n

        real (wp) c(in,*)
        real (wp) ch(m,*)
        real (wp) fac(15)
        real (wp) half
        real (wp) halfm
        integer (ip) i
        integer (ip) idl1
        integer (ip) ido
        integer (ip) im
        integer (ip) iip
        integer (ip) iw
        integer (ip) ix2
        integer (ip) ix3
        integer (ip) ix4
        integer (ip) j
        integer (ip) k1
        integer (ip) l1
        integer (ip) l2
        integer (ip) m2
        integer (ip) modn
        integer (ip) na
        integer (ip) nf
        integer (ip) nl
        real (wp) wa(n)

        nf = int(fac(2), kind=ip)
        na = 0

        do k1=1,nf
            iip = int(fac(k1+2), kind=ip)
            na = 1-na
            if (iip <= 5) then
                cycle
            end if
            if (k1 == nf) then
                cycle
            end if
            na = 1-na
        end do

        half = 0.5_wp
        halfm = -0.5_wp
        modn = mod(n,2)

        if (modn /= 0) then
            nl = n-1
        else
            nl = n-2
        end if

        if (na /= 0) then
            m2 = 1-im
            do i=1,m
                m2 = m2+im
                ch(i,1) = c(m2,1)
                ch(i,n) = c(m2,n)
            end do
            do j=2,nl,2
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    ch(i,j) = half*c(m2,j)
                    ch(i,j+1) = halfm*c(m2,j+1)
                end do
            end do
        else
            do j=2,nl,2
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,j) = half*c(m2,j)
                    c(m2,j+1) = halfm*c(m2,j+1)
                end do
            end do
        end if

        l1 = 1
        iw = 1
        do k1=1,nf
            iip = int(fac(k1+2), kind=ip)
            l2 = iip*l1
            ido = n/l2
            idl1 = ido*l1

            select case (iip)
                case (2)
                    if (na == 0) then
                        call mradb2 (m,ido,l1,c,im,in,ch,1,m,wa(iw))
                    else
                        call mradb2 (m,ido,l1,ch,1,m,c,im,in,wa(iw))
                    end if
                    na = 1-na
                case (3)
                    ix2 = iw+ido
                    if (na == 0) then
                        call mradb3(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2))
                    else
                        call mradb3(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2))
                    end if
                    na = 1-na
                case(4)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    if (na == 0) then
                        call mradb4 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3))
                    else
                        call mradb4 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3))
                    end if
                    na = 1-na
                case (5)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    ix4 = ix3+ido
                    if (na == 0) then
                        call mradb5 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    else
                        call mradb5 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    end if
                    na = 1-na
                case default
                    if (na == 0) then
                        call mradbg (m,ido,iip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
                    else
                        call mradbg (m,ido,iip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
                    end if
                    if (ido == 1) then
                        na = 1-na
                    end if
            end select
            l1 = l2
            iw = iw+(iip-1)*ido
        end do

    contains

        subroutine mradb2(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,2,l1)
            real (wp) ch(in2,ido,l1,2)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,k,1) = cc(m1,1,1,k)+cc(m1,ido,2,k)
                    ch(m2,1,k,2) = cc(m1,1,1,k)-cc(m1,ido,2,k)
                end do
            end do

            if (ido < 2) then
                return
            else if (ido == 2) then
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,ido,k,1) = cc(m1,ido,1,k)+cc(m1,ido,1,k)
                        ch(m2,ido,k,2) = -(cc(m1,1,2,k)+cc(m1,1,2,k))
                    end do
                end do
            else
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+cc(m1,ic-1,2,k)
                            ch(m2,i,k,1) = cc(m1,i,1,k)-cc(m1,ic,2,k)
                            ch(m2,i-1,k,2) = wa1(i-2)*(cc(m1,i-1,1,k)-cc(m1,ic-1,2,k)) &
                                -wa1(i-1)*(cc(m1,i,1,k)+cc(m1,ic,2,k))

                            ch(m2,i,k,2) = wa1(i-2)*(cc(m1,i,1,k)+cc(m1,ic,2,k))+wa1(i-1) &
                                *(cc(m1,i-1,1,k)-cc(m1,ic-1,2,k))
                        end do
                    end do
                end do
                if (mod(ido,2) /= 1) then
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,ido,k,1) = cc(m1,ido,1,k)+cc(m1,ido,1,k)
                            ch(m2,ido,k,2) = -(cc(m1,1,2,k)+cc(m1,1,2,k))
                        end do
                    end do
                end if
            end if

        end subroutine mradb2

        subroutine mradb3 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1


            real (wp) cc(in1,ido,3,l1)
            real (wp) ch(in2,ido,l1,3)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            real (wp), parameter :: ARG= TWO_PI/3
            real (wp), parameter :: TAUI = cos(ARG)
            real (wp), parameter :: TAUR = sin(ARG)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,k,1) = cc(m1,1,1,k)+ 2.0_wp *cc(m1,ido,2,k)
                    ch(m2,1,k,2) = cc(m1,1,1,k)+( 2.0_wp *TAUR)*cc(m1,ido,2,k) &
                        -( 2.0_wp *TAUI)*cc(m1,1,3,k)
                    ch(m2,1,k,3) = cc(m1,1,1,k)+( 2.0_wp *TAUR)*cc(m1,ido,2,k) &
                        + 2.0_wp *TAUI*cc(m1,1,3,k)
                end do
            end do

            if (ido /= 1) then
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))
                            ch(m2,i,k,1) = cc(m1,i,1,k)+(cc(m1,i,3,k)-cc(m1,ic,2,k))

                            ch(m2,i-1,k,2) = wa1(i-2)* &
                                ((cc(m1,i-1,1,k)+TAUR*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))- &
                                (TAUI*(cc(m1,i,3,k)+cc(m1,ic,2,k)))) &
                                -wa1(i-1)* &
                                ((cc(m1,i,1,k)+TAUR*(cc(m1,i,3,k)-cc(m1,ic,2,k)))+ &
                                (TAUI*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))))

                            ch(m2,i,k,2) = wa1(i-2)* &
                                ((cc(m1,i,1,k)+TAUR*(cc(m1,i,3,k)-cc(m1,ic,2,k)))+ &
                                (TAUI*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))) &
                                +wa1(i-1)* &
                                ((cc(m1,i-1,1,k)+TAUR*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))- &
                                (TAUI*(cc(m1,i,3,k)+cc(m1,ic,2,k))))

                            ch(m2,i-1,k,3) = wa2(i-2)* &
                                ((cc(m1,i-1,1,k)+TAUR*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))+ &
                                (TAUI*(cc(m1,i,3,k)+cc(m1,ic,2,k)))) &
                                -wa2(i-1)* &
                                ((cc(m1,i,1,k)+TAUR*(cc(m1,i,3,k)-cc(m1,ic,2,k)))- &
                                (TAUI*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))))

                            ch(m2,i,k,3) = wa2(i-2)* &
                                ((cc(m1,i,1,k)+TAUR*(cc(m1,i,3,k)-cc(m1,ic,2,k)))- &
                                (TAUI*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))) &
                                +wa2(i-1)* &
                                ((cc(m1,i-1,1,k)+TAUR*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))+ &
                                (TAUI*(cc(m1,i,3,k)+cc(m1,ic,2,k))))
                        end do
                    end do
                end do
            end if

        end subroutine mradb3

        subroutine mradb4(m, ido, l1, cc, im1, in1, ch, im2, in2, wa1, wa2, wa3)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,4,l1)
            real (wp) ch(in2,ido,l1,4)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp), parameter :: SQRT2 = sqrt(2.0_wp)
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp) wa3(ido)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,k,3) = (cc(m1,1,1,k)+cc(m1,ido,4,k)) &
                        -(cc(m1,ido,2,k)+cc(m1,ido,2,k))
                    ch(m2,1,k,1) = (cc(m1,1,1,k)+cc(m1,ido,4,k)) &
                        +(cc(m1,ido,2,k)+cc(m1,ido,2,k))
                    ch(m2,1,k,4) = (cc(m1,1,1,k)-cc(m1,ido,4,k)) &
                        +(cc(m1,1,3,k)+cc(m1,1,3,k))
                    ch(m2,1,k,2) = (cc(m1,1,1,k)-cc(m1,ido,4,k)) &
                        -(cc(m1,1,3,k)+cc(m1,1,3,k))
                end do
            end do

            if (ido < 2) then
                return
            else if (ido == 2) then
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,ido,k,1) = (cc(m1,ido,1,k)+cc(m1,ido,3,k)) &
                            +(cc(m1,ido,1,k)+cc(m1,ido,3,k))
                        ch(m2,ido,k,2) = SQRT2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
                            -(cc(m1,1,2,k)+cc(m1,1,4,k)))
                        ch(m2,ido,k,3) = (cc(m1,1,4,k)-cc(m1,1,2,k)) &
                            +(cc(m1,1,4,k)-cc(m1,1,2,k))
                        ch(m2,ido,k,4) = -SQRT2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
                            +(cc(m1,1,2,k)+cc(m1,1,4,k)))
                    end do
                end do
            else
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,k,1) = (cc(m1,i-1,1,k)+cc(m1,ic-1,4,k)) &
                                +(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))
                            ch(m2,i,k,1) = (cc(m1,i,1,k)-cc(m1,ic,4,k)) &
                                +(cc(m1,i,3,k)-cc(m1,ic,2,k))
                            ch(m2,i-1,k,2)=wa1(i-2)*((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k)) &
                                -(cc(m1,i,3,k)+cc(m1,ic,2,k)))-wa1(i-1) &
                                *((cc(m1,i,1,k)+cc(m1,ic,4,k))+(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))
                            ch(m2,i,k,2)=wa1(i-2)*((cc(m1,i,1,k)+cc(m1,ic,4,k)) &
                                +(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))+wa1(i-1) &
                                *((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k))-(cc(m1,i,3,k)+cc(m1,ic,2,k)))
                            ch(m2,i-1,k,3)=wa2(i-2)*((cc(m1,i-1,1,k)+cc(m1,ic-1,4,k)) &
                                -(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))-wa2(i-1) &
                                *((cc(m1,i,1,k)-cc(m1,ic,4,k))-(cc(m1,i,3,k)-cc(m1,ic,2,k)))
                            ch(m2,i,k,3)=wa2(i-2)*((cc(m1,i,1,k)-cc(m1,ic,4,k)) &
                                -(cc(m1,i,3,k)-cc(m1,ic,2,k)))+wa2(i-1) &
                                *((cc(m1,i-1,1,k)+cc(m1,ic-1,4,k))-(cc(m1,i-1,3,k) &
                                +cc(m1,ic-1,2,k)))
                            ch(m2,i-1,k,4)=wa3(i-2)*((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k)) &
                                +(cc(m1,i,3,k)+cc(m1,ic,2,k)))-wa3(i-1) &
                                *((cc(m1,i,1,k)+cc(m1,ic,4,k))-(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))
                            ch(m2,i,k,4)=wa3(i-2)*((cc(m1,i,1,k)+cc(m1,ic,4,k)) &
                                -(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))+wa3(i-1) &
                                *((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k))+(cc(m1,i,3,k)+cc(m1,ic,2,k)))
                        end do
                    end do
                end do
                if (mod(ido,2) /= 1) then
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,ido,k,1) = (cc(m1,ido,1,k)+cc(m1,ido,3,k)) &
                                +(cc(m1,ido,1,k)+cc(m1,ido,3,k))
                            ch(m2,ido,k,2) = SQRT2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
                                -(cc(m1,1,2,k)+cc(m1,1,4,k)))
                            ch(m2,ido,k,3) = (cc(m1,1,4,k)-cc(m1,1,2,k)) &
                                +(cc(m1,1,4,k)-cc(m1,1,2,k))
                            ch(m2,ido,k,4) = -SQRT2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
                                +(cc(m1,1,2,k)+cc(m1,1,4,k)))
                        end do
                    end do
                end if
            end if

        end subroutine mradb4

        subroutine mradb5(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3,wa4)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,5,l1)
            real (wp) ch(in2,ido,l1,5)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp) wa3(ido)
            real (wp) wa4(ido)

            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            real (wp), parameter :: ARG = TWO_PI/5
            real (wp), parameter :: TR11=cos(ARG)
            real (wp), parameter :: TI11=sin(ARG)
            real (wp), parameter :: TR12=cos(2.0_wp*ARG)
            real (wp), parameter :: TI12=sin(2.0_wp*ARG)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,k,1) = cc(m1,1,1,k)+ 2.0_wp *cc(m1,ido,2,k)&
                        + 2.0_wp *cc(m1,ido,4,k)
                    ch(m2,1,k,2) = (cc(m1,1,1,k)+TR11* 2.0_wp *cc(m1,ido,2,k) &
                        +TR12* 2.0_wp *cc(m1,ido,4,k))-(TI11* 2.0_wp *cc(m1,1,3,k) &
                        +TI12* 2.0_wp *cc(m1,1,5,k))
                    ch(m2,1,k,3) = (cc(m1,1,1,k)+TR12* 2.0_wp *cc(m1,ido,2,k) &
                        +TR11* 2.0_wp *cc(m1,ido,4,k))-(TI12* 2.0_wp *cc(m1,1,3,k) &
                        -TI11* 2.0_wp *cc(m1,1,5,k))
                    ch(m2,1,k,4) = (cc(m1,1,1,k)+TR12* 2.0_wp *cc(m1,ido,2,k) &
                        +TR11* 2.0_wp *cc(m1,ido,4,k))+(TI12* 2.0_wp *cc(m1,1,3,k) &
                        -TI11* 2.0_wp *cc(m1,1,5,k))
                    ch(m2,1,k,5) = (cc(m1,1,1,k)+TR11* 2.0_wp *cc(m1,ido,2,k) &
                        +TR12* 2.0_wp *cc(m1,ido,4,k))+(TI11* 2.0_wp *cc(m1,1,3,k) &
                        +TI12* 2.0_wp *cc(m1,1,5,k))
                end do
            end do

            if (ido /= 1) then
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k))
                            ch(m2,i,k,1) = cc(m1,i,1,k)+(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
                                +(cc(m1,i,5,k)-cc(m1,ic,4,k))
                            ch(m2,i-1,k,2) = wa1(i-2)*((cc(m1,i-1,1,k)+TR11* &
                                (cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))+TR12 &
                                *(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(TI11*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))+TI12*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
                                -wa1(i-1)*((cc(m1,i,1,k)+TR11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
                                +TR12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))+(TI11*(cc(m1,i-1,3,k) &
                                -cc(m1,ic-1,2,k))+TI12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

                            ch(m2,i,k,2) = wa1(i-2)*((cc(m1,i,1,k)+TR11*(cc(m1,i,3,k) &
                                -cc(m1,ic,2,k))+TR12*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                +(TI11*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))+TI12 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))+wa1(i-1) &
                                *((cc(m1,i-1,1,k)+TR11*(cc(m1,i-1,3,k) &
                                +cc(m1,ic-1,2,k))+TR12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k))) &
                                -(TI11*(cc(m1,i,3,k)+cc(m1,ic,2,k))+TI12 &
                                *(cc(m1,i,5,k)+cc(m1,ic,4,k))))
                            ch(m2,i-1,k,3) = wa2(i-2) &
                                *((cc(m1,i-1,1,k)+TR12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(TI12*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))-TI11*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
                                -wa2(i-1) &
                                *((cc(m1,i,1,k)+TR12*(cc(m1,i,3,k)- &
                                cc(m1,ic,2,k))+TR11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                +(TI12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-TI11 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

                            ch(m2,i,k,3) = wa2(i-2) &
                                *((cc(m1,i,1,k)+TR12*(cc(m1,i,3,k)- &
                                cc(m1,ic,2,k))+TR11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                +(TI12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-TI11 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
                                +wa2(i-1) &
                                *((cc(m1,i-1,1,k)+TR12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(TI12*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))-TI11*(cc(m1,i,5,k)+cc(m1,ic,4,k))))

                            ch(m2,i-1,k,4) = wa3(i-2) &
                                *((cc(m1,i-1,1,k)+TR12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(TI12*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))-TI11*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
                                -wa3(i-1) &
                                *((cc(m1,i,1,k)+TR12*(cc(m1,i,3,k)- &
                                cc(m1,ic,2,k))+TR11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                -(TI12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-TI11 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

                            ch(m2,i,k,4) = wa3(i-2) &
                                *((cc(m1,i,1,k)+TR12*(cc(m1,i,3,k)- &
                                cc(m1,ic,2,k))+TR11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                -(TI12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-TI11 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
                                +wa3(i-1) &
                                *((cc(m1,i-1,1,k)+TR12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(TI12*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))-TI11*(cc(m1,i,5,k)+cc(m1,ic,4,k))))

                            ch(m2,i-1,k,5) = wa4(i-2) &
                                *((cc(m1,i-1,1,k)+TR11*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(TI11*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))+TI12*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
                                -wa4(i-1) &
                                *((cc(m1,i,1,k)+TR11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
                                +TR12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))-(TI11*(cc(m1,i-1,3,k) &
                                -cc(m1,ic-1,2,k))+TI12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

                            ch(m2,i,k,5) = wa4(i-2) &
                                *((cc(m1,i,1,k)+TR11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
                                +TR12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))-(TI11*(cc(m1,i-1,3,k) &
                                -cc(m1,ic-1,2,k))+TI12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
                                +wa4(i-1) &
                                *((cc(m1,i-1,1,k)+TR11*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(TI11*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))+TI12*(cc(m1,i,5,k)+cc(m1,ic,4,k))))
                        end do
                    end do
                end do
            end if

        end subroutine mradb5

        subroutine mradbg (m,ido,iip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

            integer (ip) idl1
            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) iip
            integer (ip) l1

            real (wp) ai1
            real (wp) ai2
            real (wp) ar1
            real (wp) ar1h
            real (wp) ar2
            real (wp) ar2h, arg
            real (wp) c1(in1,ido,l1,iip)
            real (wp) c2(in1,idl1,iip)
            real (wp) cc(in1,ido,iip,l1)
            real (wp) ch(in2,ido,l1,iip)
            real (wp) ch2(in2,idl1,iip)
            real (wp) dc2, dcp
            real (wp) ds2, dsp
            integer (ip) i
            integer (ip) ic
            integer (ip) idij
            integer (ip) idp2
            integer (ip) ik
            integer (ip) im1
            integer (ip) im2
            integer (ip) iipp2
            integer (ip) iipph
            integer (ip) is
            integer (ip) j
            integer (ip) j2
            integer (ip) jc
            integer (ip) k
            integer (ip) l
            integer (ip) lc
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) nbd
            real (wp) wa(ido)

            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)

            arg = TWO_PI/iip
            dcp = cos(arg)
            dsp = sin(arg)

            m1d = (m - 1) * im1 + 1
            m2s = 1 - im2

            idp2 = ido + 2
            nbd = (ido-1)/2
            iipp2 = iip+2
            iipph = (iip+1)/2

            if (ido >= l1) then
                do k=1,l1
                    do i=1,ido
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i,k,1) = cc(m1,i,1,k)
                        end do
                    end do
                end do
            else
                do i=1,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i,k,1) = cc(m1,i,1,k)
                        end do
                    end do
                end do
            end if

            do j=2,iipph
                jc = iipp2-j
                j2 = j+j
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,1,k,j) = cc(m1,ido,j2-2,k)+cc(m1,ido,j2-2,k)
                        ch(m2,1,k,jc) = cc(m1,1,j2-1,k)+cc(m1,1,j2-1,k)
                    end do
                end do
            end do

            if (ido /= 1) then
                if (nbd >= l1) then
                    do j=2,iipph
                        jc = iipp2-j
                        do k=1,l1
                            do i=3,ido,2
                                ic = idp2-i
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = cc(m1,i-1,2*j-1,k)+cc(m1,ic-1,2*j-2,k)
                                    ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k)-cc(m1,ic-1,2*j-2,k)
                                    ch(m2,i,k,j) = cc(m1,i,2*j-1,k)-cc(m1,ic,2*j-2,k)
                                    ch(m2,i,k,jc) = cc(m1,i,2*j-1,k)+cc(m1,ic,2*j-2,k)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2,iipph
                        jc = iipp2-j
                        do i=3,ido,2
                            ic = idp2-i
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = cc(m1,i-1,2*j-1,k)+cc(m1,ic-1,2*j-2,k)
                                    ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k)-cc(m1,ic-1,2*j-2,k)
                                    ch(m2,i,k,j) = cc(m1,i,2*j-1,k)-cc(m1,ic,2*j-2,k)
                                    ch(m2,i,k,jc) = cc(m1,i,2*j-1,k)+cc(m1,ic,2*j-2,k)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

            ar1 = 1.0_wp
            ai1 = 0.0_wp
            do l=2,iipph
                lc = iipp2-l
                ar1h = dcp*ar1-dsp*ai1
                ai1 = dcp*ai1+dsp*ar1
                ar1 = ar1h
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        c2(m1,ik,l) = ch2(m2,ik,1)+ar1*ch2(m2,ik,2)
                        c2(m1,ik,lc) = ai1*ch2(m2,ik,iip)
                    end do
                end do
                dc2 = ar1
                ds2 = ai1
                ar2 = ar1
                ai2 = ai1
                do j=3,iipph
                    jc = iipp2-j
                    ar2h = dc2*ar2-ds2*ai2
                    ai2 = dc2*ai2+ds2*ar2
                    ar2 = ar2h
                    do ik=1,idl1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            c2(m1,ik,l) = c2(m1,ik,l)+ar2*ch2(m2,ik,j)
                            c2(m1,ik,lc) = c2(m1,ik,lc)+ai2*ch2(m2,ik,jc)
                        end do
                    end do
                end do
            end do

            do j=2,iipph
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch2(m2,ik,1) = ch2(m2,ik,1)+ch2(m2,ik,j)
                    end do
                end do
            end do

            do j=2,iipph
                jc = iipp2-j
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,1,k,j) = c1(m1,1,k,j)-c1(m1,1,k,jc)
                        ch(m2,1,k,jc) = c1(m1,1,k,j)+c1(m1,1,k,jc)
                    end do
                end do
            end do

            if (ido /= 1) then
                if (nbd >= l1) then
                    do j=2,iipph
                        jc = iipp2-j
                        do k=1,l1
                            do i=3,ido,2
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = c1(m1,i-1,k,j)-c1(m1,i,k,jc)
                                    ch(m2,i-1,k,jc) = c1(m1,i-1,k,j)+c1(m1,i,k,jc)
                                    ch(m2,i,k,j) = c1(m1,i,k,j)+c1(m1,i-1,k,jc)
                                    ch(m2,i,k,jc) = c1(m1,i,k,j)-c1(m1,i-1,k,jc)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2,iipph
                        jc = iipp2-j
                        do i=3,ido,2
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = c1(m1,i-1,k,j)-c1(m1,i,k,jc)
                                    ch(m2,i-1,k,jc) = c1(m1,i-1,k,j)+c1(m1,i,k,jc)
                                    ch(m2,i,k,j) = c1(m1,i,k,j)+c1(m1,i-1,k,jc)
                                    ch(m2,i,k,jc) = c1(m1,i,k,j)-c1(m1,i-1,k,jc)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

            if (ido /= 1) then
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        c2(m1,ik,1) = ch2(m2,ik,1)
                    end do
                end do
                do j=2,iip
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            c1(m1,1,k,j) = ch(m2,1,k,j)
                        end do
                    end do
                end do
                if (l1 >= nbd ) then
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        idij = is
                        do i=3,ido,2
                            idij = idij+2
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    c1(m1,i-1,k,j) = &
                                        wa(idij-1)*ch(m2,i-1,k,j) &
                                        - wa(idij)* ch(m2,i,k,j)
                                    c1(m1,i,k,j) = &
                                        wa(idij-1)*ch(m2,i,k,j) &
                                        + wa(idij)* ch(m2,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                else
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        do k=1,l1
                            idij = is
                            do i=3,ido,2
                                idij = idij+2
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    c1(m1,i-1,k,j) = &
                                        wa(idij-1)*ch(m2,i-1,k,j)&
                                        - wa(idij)*ch(m2,i,k,j)
                                    c1(m1,i,k,j) = &
                                        wa(idij-1)*ch(m2,i,k,j)&
                                        + wa(idij)*ch(m2,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

        end subroutine mradbg

    end subroutine mrftb1


    subroutine mrftf1(m,im,n,in,c,ch,wa,fac)

        integer (ip) in
        integer (ip) m
        integer (ip) n

        real (wp) c(in,*)
        real (wp) ch(m,*)
        real (wp) fac(15)
        integer (ip) i
        integer (ip) idl1
        integer (ip) ido
        integer (ip) im
        integer (ip) iip
        integer (ip) iw
        integer (ip) ix2
        integer (ip) ix3
        integer (ip) ix4
        integer (ip) j
        integer (ip) k1
        integer (ip) kh
        integer (ip) l1
        integer (ip) l2
        integer (ip) m2
        integer (ip) modn
        integer (ip) na
        integer (ip) nf
        integer (ip) nl
        real (wp) sn
        real (wp) tsn
        real (wp) tsnm
        real (wp) wa(n)

        nf = int(fac(2), kind=ip)
        na = 1
        l2 = n
        iw = n

        do k1=1,nf
            kh = nf-k1
            iip = int(fac(kh+3), kind=ip)
            l1 = l2/iip
            ido = n/l2
            idl1 = ido*l1
            iw = iw-(iip-1)*ido
            na = 1-na
            select case (iip)
                case (2)
                    if (na == 0) then
                        call mradf2(m,ido,l1,c,im,in,ch,1,m,wa(iw))
                    else
                        call mradf2(m,ido,l1,ch,1,m,c,im,in,wa(iw))
                    end if
                case(3)
                    ix2 = iw+ido
                    if (na == 0) then
                        call mradf3(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2))
                    else
                        call mradf3(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2))
                    end if
                case (4)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    if (na == 0) then
                        call mradf4(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3))
                    else
                        call mradf4(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3))
                    end if
                case(5)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    ix4 = ix3+ido
                    if (na == 0) then
                        call mradf5(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    else
                        call mradf5(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    end if
                case default
                    if (ido == 1) then
                        na = 1-na
                    end if
                    if (na == 0) then
                        call mradfg(m,ido,iip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
                        na = 1
                    else
                        call mradfg(m,ido,iip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
                        na = 0
                    end if
            end select
            l2 = l1
        end do

        sn = 1.0_wp/n
        tsn =  2.0_wp/n
        tsnm = -tsn
        modn = mod(n,2)

        if (modn /= 0) then
            nl = n-1
        else
            nl = n-2
        end if

        if (na == 0) then
            m2 = 1-im

            do i=1,m
                m2 = m2+im
                c(m2,1) = sn*ch(i,1)
            end do

            do j=2,nl,2
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,j) = tsn*ch(i,j)
                    c(m2,j+1) = tsnm*ch(i,j+1)
                end do
            end do

            if (modn == 0) then
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,n) = sn*ch(i,n)
                end do
            end if
        else
            m2 = 1-im
            do i=1,m
                m2 = m2+im
                c(m2,1) = sn*c(m2,1)
            end do
            do j=2,nl,2
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,j) = tsn*c(m2,j)
                    c(m2,j+1) = tsnm*c(m2,j+1)
                end do
            end do

            if (modn == 0) then
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,n) = sn*c(m2,n)
                end do
            end if
        end if

    contains

        subroutine mradf2(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,l1,2)
            real (wp) ch(in2,ido,2,l1)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,1,k) = cc(m1,1,k,1)+cc(m1,1,k,2)
                    ch(m2,ido,2,k) = cc(m1,1,k,1)-cc(m1,1,k,2)
                end do
            end do

            if (ido < 2) then
                return
            else if (ido == 2) then
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,1,2,k) = -cc(m1,ido,k,2)
                        ch(m2,ido,1,k) = cc(m1,ido,k,1)
                    end do
                end do
            else
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i,1,k) = cc(m1,i,k,1)+(wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))
                            ch(m2,ic,2,k) = (wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))-cc(m1,i,k,1)
                            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+(wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))
                            ch(m2,ic-1,2,k) = cc(m1,i-1,k,1)-(wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))
                        end do
                    end do
                end do
                if (mod(ido,2) /= 1) then
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,1,2,k) = -cc(m1,ido,k,2)
                            ch(m2,ido,1,k) = cc(m1,ido,k,1)
                        end do
                    end do
                end if
            end if

        end subroutine mradf2

        subroutine mradf3(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,l1,3)
            real (wp) ch(in2,ido,3,l1)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp), parameter :: TWO_PI =2.0_wp * acos(-1.0_wp)
            real (wp), parameter :: ARG = TWO_PI/3
            real (wp), parameter :: TAUR = cos(ARG)
            real (wp), parameter :: TAUI = sin(ARG)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,1,k) = cc(m1,1,k,1)+(cc(m1,1,k,2)+cc(m1,1,k,3))
                    ch(m2,1,3,k) = TAUI*(cc(m1,1,k,3)-cc(m1,1,k,2))
                    ch(m2,ido,2,k) = cc(m1,1,k,1)+TAUR*(cc(m1,1,k,2)+cc(m1,1,k,3))
                end do
            end do

            if (ido /= 1) then
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3)))

                            ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3)))

                            ch(m2,i-1,3,k) = (cc(m1,i-1,k,1)+TAUR*((wa1(i-2)* &
                                cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
                                cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))))+(TAUI*((wa1(i-2)* &
                                cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
                                cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))

                            ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)+TAUR*((wa1(i-2)* &
                                cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
                                cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))))-(TAUI*((wa1(i-2)* &
                                cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
                                cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))

                            ch(m2,i,3,k) = (cc(m1,i,k,1)+TAUR*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))))+(TAUI*((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2))))

                            ch(m2,ic,2,k) = (TAUI*((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2))))-(cc(m1,i,k,1)+TAUR*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))))
                        end do
                    end do
                end do
            end if

        end subroutine mradf3


        subroutine mradf4(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,l1,4)
            real (wp) ch(in2,ido,4,l1)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp) wa3(ido)
            real (wp), parameter :: HALF_SQRT2 = sqrt(2.0_wp)/2

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,1,k) = (cc(m1,1,k,2)+cc(m1,1,k,4)) &
                        +(cc(m1,1,k,1)+cc(m1,1,k,3))
                    ch(m2,ido,4,k) = (cc(m1,1,k,1)+cc(m1,1,k,3)) &
                        -(cc(m1,1,k,2)+cc(m1,1,k,4))
                    ch(m2,ido,2,k) = cc(m1,1,k,1)-cc(m1,1,k,3)
                    ch(m2,1,3,k) = cc(m1,1,k,4)-cc(m1,1,k,2)
                end do
            end do

            if (ido < 2) then
                return
            else if (ido == 2) then
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,ido,1,k) = (HALF_SQRT2*(cc(m1,ido,k,2)-cc(m1,ido,k,4)))+ &
                            cc(m1,ido,k,1)
                        ch(m2,ido,3,k) = cc(m1,ido,k,1)-(HALF_SQRT2*(cc(m1,ido,k,2)- &
                            cc(m1,ido,k,4)))
                        ch(m2,1,2,k) = (-HALF_SQRT2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))- &
                            cc(m1,ido,k,3)
                        ch(m2,1,4,k) = (-HALF_SQRT2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))+ &
                            cc(m1,ido,k,3)
                    end do
                end do
            else
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,1,k) = ((wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2))+(wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4)))+(cc(m1,i-1,k,1)+(wa2(i-2)*cc(m1,i-1,k,3)+ &
                                wa2(i-1)*cc(m1,i,k,3)))

                            ch(m2,ic-1,4,k) = (cc(m1,i-1,k,1)+(wa2(i-2)*cc(m1,i-1,k,3)+ &
                                wa2(i-1)*cc(m1,i,k,3)))-((wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))+(wa3(i-2)*cc(m1,i-1,k,4)+ &
                                wa3(i-1)*cc(m1,i,k,4)))

                            ch(m2,i,1,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))+(cc(m1,i,k,1)+(wa2(i-2)*cc(m1,i,k,3)- &
                                wa2(i-1)*cc(m1,i-1,k,3)))

                            ch(m2,ic,4,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))-(cc(m1,i,k,1)+(wa2(i-2)*cc(m1,i,k,3)- &
                                wa2(i-1)*cc(m1,i-1,k,3)))

                            ch(m2,i-1,3,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))-(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))+(cc(m1,i-1,k,1)-(wa2(i-2)*cc(m1,i-1,k,3)+ &
                                wa2(i-1)*cc(m1,i,k,3)))

                            ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)-(wa2(i-2)*cc(m1,i-1,k,3)+ &
                                wa2(i-1)*cc(m1,i,k,3)))-((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))-(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))

                            ch(m2,i,3,k) = ((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))+(cc(m1,i,k,1)-(wa2(i-2)*cc(m1,i,k,3)- &
                                wa2(i-1)*cc(m1,i-1,k,3)))

                            ch(m2,ic,2,k) = ((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))-(cc(m1,i,k,1)-(wa2(i-2)*cc(m1,i,k,3)- &
                                wa2(i-1)*cc(m1,i-1,k,3)))
                        end do
                    end do
                end do
                if (mod(ido,2) /= 1) then
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,ido,1,k) = (HALF_SQRT2*(cc(m1,ido,k,2)-cc(m1,ido,k,4)))+ &
                                cc(m1,ido,k,1)
                            ch(m2,ido,3,k) = cc(m1,ido,k,1)-(HALF_SQRT2*(cc(m1,ido,k,2)- &
                                cc(m1,ido,k,4)))
                            ch(m2,1,2,k) = (-HALF_SQRT2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))- &
                                cc(m1,ido,k,3)
                            ch(m2,1,4,k) = (-HALF_SQRT2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))+ &
                                cc(m1,ido,k,3)
                        end do
                    end do
                end if
            end if

        end subroutine mradf4

        subroutine mradf5(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3,wa4)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,l1,5)
            real (wp) ch(in2,ido,5,l1)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp) wa3(ido)
            real (wp) wa4(ido)
            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            real (wp), parameter :: ARG = TWO_PI/5
            real (wp), parameter :: TR11 = cos(ARG)
            real (wp), parameter :: TI11 = sin(ARG)
            real (wp), parameter :: TR12 = cos(2.0_wp*ARG)
            real (wp), parameter :: TI12 = sin(2.0_wp*ARG)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,1,k) = cc(m1,1,k,1)+(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
                        (cc(m1,1,k,4)+cc(m1,1,k,3))
                    ch(m2,ido,2,k) = cc(m1,1,k,1)+TR11*(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
                        TR12*(cc(m1,1,k,4)+cc(m1,1,k,3))
                    ch(m2,1,3,k) = TI11*(cc(m1,1,k,5)-cc(m1,1,k,2))+TI12* &
                        (cc(m1,1,k,4)-cc(m1,1,k,3))
                    ch(m2,ido,4,k) = cc(m1,1,k,1)+TR12*(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
                        TR11*(cc(m1,1,k,4)+cc(m1,1,k,3))
                    ch(m2,1,5,k) = TI12*(cc(m1,1,k,5)-cc(m1,1,k,2))-TI11* &
                        (cc(m1,1,k,4)-cc(m1,1,k,3))
                end do
            end do

            if (ido /= 1) then
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2

                            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
                                cc(m1,i,k,5)))+((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))+(wa3(i-2)*cc(m1,i-1,k,4)+ &
                                wa3(i-1)*cc(m1,i,k,4)))

                            ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))

                            ch(m2,i-1,3,k) = cc(m1,i-1,k,1)+TR11* &
                                ( wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2) &
                                +wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5))+TR12* &
                                ( wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3) &
                                +wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))+TI11* &
                                ( wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2) &
                                -(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))+TI12* &
                                ( wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3) &
                                -(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4)))

                            ch(m2,ic-1,2,k) = cc(m1,i-1,k,1)+TR11* &
                                ( wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2) &
                                +wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5))+TR12* &
                                ( wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3) &
                                +wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))-(TI11* &
                                ( wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2) &
                                -(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))+TI12* &
                                ( wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3) &
                                -(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

                            ch(m2,i,3,k) = (cc(m1,i,k,1)+TR11*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+TR12*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4))))+(TI11*((wa4(i-2)*cc(m1,i-1,k,5)+ &
                                wa4(i-1)*cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))+TI12*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))))

                            ch(m2,ic,2,k) = (TI11*((wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
                                cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))+TI12*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))))-(cc(m1,i,k,1)+TR11*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+TR12*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4))))

                            ch(m2,i-1,5,k) = (cc(m1,i-1,k,1)+TR12*((wa1(i-2)* &
                                cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)* &
                                cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5)))+TR11*((wa2(i-2)* &
                                cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))+(wa3(i-2)* &
                                cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))))+(TI12*((wa1(i-2)* &
                                cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa4(i-2)* &
                                cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))-TI11*((wa2(i-2)* &
                                cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))-(wa3(i-2)* &
                                cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

                            ch(m2,ic-1,4,k) = (cc(m1,i-1,k,1)+TR12*((wa1(i-2)* &
                                cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)* &
                                cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5)))+TR11*((wa2(i-2)* &
                                cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))+(wa3(i-2)* &
                                cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))))-(TI12*((wa1(i-2)* &
                                cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa4(i-2)* &
                                cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))-TI11*((wa2(i-2)* &
                                cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))-(wa3(i-2)* &
                                cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

                            ch(m2,i,5,k) = (cc(m1,i,k,1)+TR12*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+TR11*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4))))+(TI12*((wa4(i-2)*cc(m1,i-1,k,5)+ &
                                wa4(i-1)*cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))-TI11*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))))

                            ch(m2,ic,4,k) = (TI12*((wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
                                cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))-TI11*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))))-(cc(m1,i,k,1)+TR12*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+TR11*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4))))
                        end do
                    end do
                end do
            end if

        end subroutine mradf5

        subroutine mradfg(m,ido,iip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

            integer (ip) idl1
            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) iip
            integer (ip) l1

            real (wp) ai1
            real (wp) ai2
            real (wp) ar1
            real (wp) ar1h
            real (wp) ar2
            real (wp) ar2h
            real (wp) arg
            real (wp) c1(in1,ido,l1,iip)
            real (wp) c2(in1,idl1,iip)
            real (wp) cc(in1,ido,iip,l1)
            real (wp) ch(in2,ido,l1,iip)
            real (wp) ch2(in2,idl1,iip)
            real (wp) dc2
            real (wp) dcp
            real (wp) ds2
            real (wp) dsp
            integer (ip) i
            integer (ip) ic
            integer (ip) idij
            integer (ip) idp2
            integer (ip) ik
            integer (ip) im1
            integer (ip) im2
            integer (ip) iipp2
            integer (ip) iipph
            integer (ip) is
            integer (ip) j
            integer (ip) j2
            integer (ip) jc
            integer (ip) k
            integer (ip) l
            integer (ip) lc
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) nbd
            real (wp), parameter :: TWO_PI= 2.0_wp * acos(-1.0_wp)
            real (wp) wa(ido)

            m1d = (m-1)*im1+1
            m2s = 1-im2
            arg = TWO_PI / iip
            dcp = cos(arg)
            dsp = sin(arg)
            iipph = (iip+1)/2
            iipp2 = iip+2
            idp2 = ido+2
            nbd = (ido-1)/2

            if (ido /= 1) then

                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch2(m2,ik,1) = c2(m1,ik,1)
                    end do
                end do

                do j=2,iip
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,1,k,j) = c1(m1,1,k,j)
                        end do
                    end do
                end do

                if ( l1 >= nbd ) then
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        idij = is
                        do i=3,ido,2
                            idij = idij+2
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = wa(idij-1)*c1(m1,i-1,k,j)+wa(idij)*c1(m1,i,k,j)
                                    ch(m2,i,k,j) = wa(idij-1)*c1(m1,i,k,j)-wa(idij)*c1(m1,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                else
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        do k=1,l1
                            idij = is
                            do i=3,ido,2
                                idij = idij+2
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = wa(idij-1)*c1(m1,i-1,k,j)+wa(idij)*c1(m1,i,k,j)
                                    ch(m2,i,k,j) = wa(idij-1)*c1(m1,i,k,j)-wa(idij)*c1(m1,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                end if

                if (nbd >= l1) then
                    do j=2,iipph
                        jc = iipp2-j
                        do k=1,l1
                            do i=3,ido,2
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    c1(m1,i-1,k,j) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
                                    c1(m1,i-1,k,jc) = ch(m2,i,k,j)-ch(m2,i,k,jc)
                                    c1(m1,i,k,j) = ch(m2,i,k,j)+ch(m2,i,k,jc)
                                    c1(m1,i,k,jc) = ch(m2,i-1,k,jc)-ch(m2,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2,iipph
                        jc = iipp2-j
                        do i=3,ido,2
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    c1(m1,i-1,k,j) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
                                    c1(m1,i-1,k,jc) = ch(m2,i,k,j)-ch(m2,i,k,jc)
                                    c1(m1,i,k,j) = ch(m2,i,k,j)+ch(m2,i,k,jc)
                                    c1(m1,i,k,jc) = ch(m2,i-1,k,jc)-ch(m2,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                end if
            else
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        c2(m1,ik,1) = ch2(m2,ik,1)
                    end do
                end do
            end if
            do j=2,iipph
                jc = iipp2-j
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        c1(m1,1,k,j) = ch(m2,1,k,j)+ch(m2,1,k,jc)
                        c1(m1,1,k,jc) = ch(m2,1,k,jc)-ch(m2,1,k,j)
                    end do
                end do
            end do

            ar1 = 1.0_wp
            ai1 = 0.0_wp
            do l=2,iipph
                lc = iipp2-l
                ar1h = dcp*ar1-dsp*ai1
                ai1 = dcp*ai1+dsp*ar1
                ar1 = ar1h
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch2(m2,ik,l) = c2(m1,ik,1)+ar1*c2(m1,ik,2)
                        ch2(m2,ik,lc) = ai1*c2(m1,ik,iip)
                    end do
                end do
                dc2 = ar1
                ds2 = ai1
                ar2 = ar1
                ai2 = ai1
                do j=3,iipph
                    jc = iipp2-j
                    ar2h = dc2*ar2-ds2*ai2
                    ai2 = dc2*ai2+ds2*ar2
                    ar2 = ar2h
                    do ik=1,idl1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch2(m2,ik,l) = ch2(m2,ik,l)+ar2*c2(m1,ik,j)
                            ch2(m2,ik,lc) = ch2(m2,ik,lc)+ai2*c2(m1,ik,jc)
                        end do
                    end do
                end do
            end do
            do j=2,iipph
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch2(m2,ik,1) = ch2(m2,ik,1)+c2(m1,ik,j)
                    end do
                end do
            end do

            if (ido >= l1) then
                do k=1,l1
                    do i=1,ido
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc(m1,i,1,k) = ch(m2,i,k,1)
                        end do
                    end do
                end do
            else
                do i=1,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc(m1,i,1,k) = ch(m2,i,k,1)
                        end do
                    end do
                end do
            end if

            do j=2,iipph
                jc = iipp2-j
                j2 = j+j
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        cc(m1,ido,j2-2,k) = ch(m2,1,k,j)
                        cc(m1,1,j2-1,k) = ch(m2,1,k,jc)
                    end do
                end do
            end do
            if (ido /= 1) then
                if (nbd >= l1) then
                    do j=2,iipph
                        jc = iipp2-j
                        j2 = j+j
                        do k=1,l1
                            do i=3,ido,2
                                ic = idp2-i
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    cc(m1,i-1,j2-1,k) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
                                    cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j)-ch(m2,i-1,k,jc)
                                    cc(m1,i,j2-1,k) = ch(m2,i,k,j)+ch(m2,i,k,jc)
                                    cc(m1,ic,j2-2,k) = ch(m2,i,k,jc)-ch(m2,i,k,j)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2,iipph
                        jc = iipp2-j
                        j2 = j+j
                        do i=3,ido,2
                            ic = idp2-i
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    cc(m1,i-1,j2-1,k) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
                                    cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j)-ch(m2,i-1,k,jc)
                                    cc(m1,i,j2-1,k) = ch(m2,i,k,j)+ch(m2,i,k,jc)
                                    cc(m1,ic,j2-2,k) = ch(m2,i,k,jc)-ch(m2,i,k,j)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

        end subroutine mradfg

    end subroutine mrftf1



    subroutine mrfti1(n,wa,fac)
        !
        !  input
        !  n, the number for which factorization and
        !  other information is needed.
        !
        !  output
        !   wa(n), trigonometric information.
        !
        !  output
        !  fac(15), factorization information. fac(1) is
        !  n, fac(2) is nf, the number of factors, and fac(3:nf+2) are the factors.
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: wa(n)
        real (wp),    intent (out) :: fac(15)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip)            :: i, ib, ido, ii, iip, ipm, is
        integer (ip)            :: j, k1, l1, l2, ld
        integer (ip)            :: nf, nfm1, nl, nq, nr, ntry
        integer (ip), parameter :: ntryh(*) = [ 4, 2, 3, 5 ]
        real (wp),    parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
        real (wp)               :: arg, argh, argld, fi
        !--------------------------------------------------------------

        ntry = 0
        nl = n
        nf = 0
        j = 0

        factorize_loop: do
            ! Increment j
            j = j+1

            ! Choose ntry
            if (j <= 4) then
                ntry = ntryh(j)
            else
                ntry = ntry+2
            end if

            inner_loop: do
                nq = nl/ntry
                nr = nl-ntry*nq
                if (nr < 0) then
                    cycle factorize_loop
                else if (nr == 0) then
                    nf = nf+1
                    fac(nf+2) = ntry
                    nl = nq

                    if (ntry == 2 .and. nf /= 1) then

                        do i=2,nf
                            ib = nf-i+2
                            fac(ib+2) = fac(ib+1)
                        end do

                        fac(3) = 2

                    end if

                    if (nl /= 1) then
                        cycle inner_loop
                    end if
                else
                    cycle factorize_loop
                end if
                exit inner_loop
            end do inner_loop
            exit factorize_loop
        end do factorize_loop

        fac(1) = n
        fac(2) = nf
        argh = TWO_PI/n
        is = 0
        nfm1 = nf-1
        l1 = 1

        do k1=1,nfm1
            iip = int(fac(k1+2), kind=ip)
            ld = 0
            l2 = l1*iip
            ido = n/l2
            ipm = iip-1

            do j=1,ipm

                ld = ld+l1
                i = is
                argld = real(ld, kind=wp) * argh
                fi = 0.0_wp
                do ii=3,ido,2
                    i = i+2
                    fi = fi + 1.0_wp
                    arg = fi*argld
                    wa(i-1) = cos(arg)
                    wa(i) = sin(arg)
                end do
                is = is+ido

            end do
            l1 = l2
        end do

    end subroutine mrfti1



    subroutine msntb1(lot,jump,n,inc,x,wsave,dsum,xh,work,ier)

        integer (ip) inc
        integer (ip) lot

        real (wp) dsum(*)
        real (wp) fnp1s4
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) k
        integer (ip) kc
        integer (ip) lj
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) lnxh
        integer (ip) m
        integer (ip) m1
        integer (ip) modn
        integer (ip) n
        integer (ip) np1
        integer (ip) ns2
        real (wp), parameter :: HALF_SQRT3 = sqrt(3.0_wp)/2
        real (wp) t1
        real (wp) t2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) xh(lot,*)
        real (wp) xhold

        ier = 0
        lj = (lot-1)*jump+1

        if (n < 2) then
            return
        else if (n == 2) then
            do m=1,lj,jump
                xhold = HALF_SQRT3*(x(m,1)+x(m,2))
                x(m,2) = HALF_SQRT3*(x(m,1)-x(m,2))
                x(m,1) = xhold
            end do
        else
            np1 = n+1
            ns2 = n/2
            do k=1,ns2
                kc = np1-k
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    t1 = x(m,k)-x(m,kc)
                    t2 = wsave(k)*(x(m,k)+x(m,kc))
                    xh(m1,k+1) = t1+t2
                    xh(m1,kc+1) = t2-t1
                end do
            end do

            modn = mod(n,2)

            if (modn /= 0) then
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    xh(m1,ns2+2) =  4.0_wp * x(m,ns2+1)
                end do
            end if

            xh(:,1) = 0.0_wp

            lnxh = lot-1 + lot*(np1-1) + 1
            lnsv = np1 + int(log(real(np1, kind=wp))/log(2.0_wp)) + 4
            lnwk = lot*np1

            call rfftmf(lot,1,np1,lot,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('msntb1',-5)
                return
            end if

            if (mod(np1,2) == 0) then
                xh(:,np1) = 2.0_wp * xh(:,np1)
            end if

            fnp1s4 = real(np1, kind=wp)/4
            m1 = 0

            do m=1,lj,jump
                m1 = m1+1
                x(m,1) = fnp1s4*xh(m1,1)
                dsum(m1) = x(m,1)
            end do

            do i=3,n,2
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    x(m,i-1) = fnp1s4*xh(m1,i)
                    dsum(m1) = dsum(m1)+fnp1s4*xh(m1,i-1)
                    x(m,i) = dsum(m1)
                end do
            end do

            if (modn == 0) then
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    x(m,n) = fnp1s4*xh(m1,n+1)
                end do
            end if
        end if

    end subroutine msntb1


    subroutine msntf1(lot,jump,n,inc,x,wsave,dsum,xh,work,ier)

        integer (ip) inc
        integer (ip) lot

        real (wp) dsum(*)
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) k
        integer (ip) kc
        integer (ip) lj
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) lnxh
        integer (ip) m
        integer (ip) m1
        integer (ip) modn
        integer (ip) n
        integer (ip) np1
        integer (ip) ns2
        real (wp) sfnp1
        real (wp) t1
        real (wp) t2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) xh(lot,*)
        real (wp) xhold

        ier = 0
        lj = (lot-1)*jump+1

        if (n < 2) then
            return
        else if (n == 2) then
            do m=1,lj,jump
                associate( sqrt3 => sqrt(3.0_wp) )
                    xhold = (x(m,1)+x(m,2))/sqrt3
                    x(m,2) = (x(m,1)-x(m,2))/sqrt3
                    x(m,1) = xhold
                end associate
            end do
        else
            np1 = n+1
            ns2 = n/2
            do k=1,ns2
                kc = np1-k
                m1 = 0
                do m=1,lj,jump
                    m1 = m1 + 1
                    t1 = x(m,k)-x(m,kc)
                    t2 = wsave(k)*(x(m,k)+x(m,kc))
                    xh(m1,k+1) = t1+t2
                    xh(m1,kc+1) = t2-t1
                end do
            end do

            modn = mod(n,2)

            if (modn /= 0) then
                m1 = 0
                do m=1,lj,jump
                    m1 = m1 + 1
                    xh(m1,ns2+2) =  4.0_wp * x(m,ns2+1)
                end do
            end if

            xh(:,1) = 0.0_wp

            lnxh = lot-1 + lot*(np1-1) + 1
            lnsv = np1 + int(log(real(np1, kind=wp))/log(2.0_wp)) + 4
            lnwk = lot*np1

            call rfftmf(lot,1,np1,lot,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,local_error_flag)

            ! Check error flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('msntf1',-5)
                return
            end if

            if (mod(np1,2) == 0) then
                xh(:,np1) = 2.0_wp * xh(:,np1)
            end if


            sfnp1 = 1.0_wp/np1
            m1 = 0

            do m=1,lj,jump
                m1 = m1+1
                x(m,1) = 0.5_wp * xh(m1,1)
                dsum(m1) = x(m,1)
            end do

            do i=3,n,2
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    x(m,i-1) = 0.5_wp * xh(m1,i)
                    dsum(m1) = dsum(m1)+ 0.5_wp * xh(m1,i-1)
                    x(m,i) = dsum(m1)
                end do
            end do

            if (modn == 0) then
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    x(m,n) = 0.5_wp * xh(m1,n+1)
                end do
            end if
        end if

    end subroutine msntf1


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



    pure subroutine mcfti1(n, wa, fnf, fac)
        !
        ! Purpose:
        !
        ! Sets up factors and tables, 64-bit float precision arithmetic.
        !
        !--------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: wa(*)
        real (wp),    intent (out) :: fnf
        real (wp),    intent (out) :: fac(*)
        !--------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------
        integer (ip) :: ido, iip, iw, k1, l1, l2, nf
        !--------------------------------------------------

        !
        !==> Get the factorization of n.
        !
        call get_factorization(n, nf, fac)
        fnf = real(nf, kind=wp)
        iw = 1
        l1 = 1
        !
        !==> Set up the trigonometric tables.
        !
        do k1 = 1, nf
            iip = int(fac(k1), kind=ip)
            l2 = l1 * iip
            ido = n / l2
            call compute_trigonometic_tables(ido, iip, wa(iw))
            iw = iw + (iip - 1) * (2*ido)
            l1 = l2
        end do

    contains

        pure subroutine get_factorization(n, nf, fac)
            !
            ! Purpose:
            !
            ! Factors of an integer for 64-bit float precision computations.
            !
            !  Parameters:
            !
            !  n, the number for which factorization and other information is needed.
            !
            !  nf, the number of factors.
            !
            !  fac(*), a list of factors of n.
            !
            !--------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------
            integer (ip), intent (in)  :: n
            integer (ip), intent (out) :: nf
            real (wp),    intent (out) :: fac(*)
            !--------------------------------------------------
            ! Dictionary: local variables
            !--------------------------------------------------
            integer (ip) :: j, nl, nq, nr, ntry
            !--------------------------------------------------

            ntry = 0
            nl = n
            nf = 0
            j = 0

            do while (1 < nl)
                j = j + 1
                select case (j)
                    case (1)
                        ntry = 4
                    case (2:3)
                        ntry = j
                    case (4)
                        ntry = 5
                    case default
                        ntry = ntry + 2
                end select

                inner_loop: do
                    nq = nl / ntry
                    nr = nl - ntry * nq

                    if ( nr /= 0 ) then
                        exit inner_loop
                    end if

                    nf = nf + 1
                    fac(nf) = real(ntry, kind=wp)
                    nl = nq

                end do inner_loop
            end do

        end subroutine get_factorization

        pure subroutine compute_trigonometic_tables(ido, iip, wa)
            !
            ! Purpose:
            !
            ! Computes trigonometric tables, 64-bit float precision arithmetic.
            !
            !--------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------
            integer (ip), intent (in)  :: ido
            integer (ip), intent (in)  :: iip
            real (wp),    intent (out) :: wa(ido,iip-1,2)
            !--------------------------------------------------
            ! Dictionary: local variables
            !--------------------------------------------------
            integer (ip)         :: i, j !! Counters
            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            real (wp)            :: argz, arg1, arg2, arg3, arg4
            !--------------------------------------------------

            argz = TWO_PI/iip
            arg1 = TWO_PI/( ido * iip)

            do j = 2, iip
                arg2 = real(j - 1, kind=wp) * arg1
                do i = 1, ido
                    arg3 = real(i - 1, kind=wp) * arg2
                    wa(i,j-1,1) = cos(arg3)
                    wa(i,j-1,2) = sin(arg3)
                end do
                if (5 < iip) then
                    arg4 = real(j - 1, kind=wp) * argz
                    wa(1,j-1,1) = cos(arg4)
                    wa(1,j-1,2) = sin(arg4)
                end if
            end do

        end subroutine compute_trigonometic_tables

    end subroutine mcfti1


    subroutine rfft1b(n, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
        !
        ! Purpose:
        !
        !  Computes the one-dimensional Fourier transform of a periodic
        !  sequence within a real array. This is referred to as the backward
        !  transform or Fourier synthesis, transforming the sequence from
        !  spectral to physical space.  This transform is normalized since a
        !  call to rfft1b followed by a call to rfft1f (or vice-versa) reproduces
        !  the original array within roundoff error.
        !
        !  Parameters:
        !
        !  integer N, the length of the sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations,
        !  in array R, of two consecutive elements within the sequence.
        !
        !  Input/real R(LENR), on input, the data to be
        !  transformed, and on output, the transformed data.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least INC*(N-1) + 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFT1I before the first call to routine
        !  RFFT1F or RFFT1B for a given transform length N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least N.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough.
        !

        integer (ip) lenr
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) inc
        integer (ip) n
        real (wp) r(lenr)
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)


        !
        !==> Check validity of calling arguments
        !
        if (lenr < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('rfft1b ', 6)
        else if (lensav < &
            n + int(log(real(n, kind=wp) )/log(2.0_wp))+4) then
            ier = 2
            call fft_error_handler('rfft1b ', 8)
        else if (lenwrk < n) then
            ier = 3
            call fft_error_handler('rfft1b ', 10)
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            call rfftb1(n,inc,r,work,wsave,wsave(n+1))
        end if

    contains

        subroutine rfftb1(n, in, c, ch, wa, fac)
            !-------------------------------------------------------------
            ! Dictionary: calling arguments
            !-------------------------------------------------------------
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: in
            real (wp),    intent (in out) :: c(in,*)
            real (wp),    intent (out)    :: ch(*)
            real (wp),    intent (in)     :: wa(n)
            real (wp),    intent (in)     :: fac(15)
            !-------------------------------------------------------------
            ! Dictionary: local variables
            !-------------------------------------------------------------
            integer (ip) :: idl1, ido, iip
            integer (ip) :: workspace_indices(4)
            integer (ip) :: j, k, l1, l2
            integer (ip) :: modn, na, nf, nl
            !-------------------------------------------------------------

            nf = int(fac(2), kind=ip)
            na = 0

            factorization_loop: do k=1,nf
                iip = int(fac(k+2), kind=ip)
                na = 1-na
                if (iip <= 5) then
                    cycle factorization_loop
                end if
                if (k == nf) then
                    cycle factorization_loop
                end if
                na = 1-na
            end do factorization_loop

            modn = mod(n,2)

            if (modn /= 0) then
                nl = n-1
            else
                nl = n-2
            end if

            if (na /= 0) then
                ch(1) = c(1,1)
                ch(n) = c(1,n)
                do j=2,nl,2
                    ch(j) = 0.5_wp*c(1,j)
                    ch(j+1) = -0.5_wp*c(1,j+1)
                end do
            else
                do j=2,nl,2
                    c(1,j) = 0.5_wp*c(1,j)
                    c(1,j+1) = -0.5_wp*c(1,j+1)
                end do
            end if

            associate( &
                iw1 => workspace_indices(1), &
                iw2 => workspace_indices(2), &
                iw3 => workspace_indices(3), &
                iw4 => workspace_indices(4) &
                )

                l1 = 1
                iw1 = 1
                do k=1,nf
                    iip = int(fac(k+2), kind=ip)
                    l2 = iip*l1
                    ido = n/l2
                    idl1 = ido*l1
                    select case (iip)
                        case (2)
                            select case (na)
                                case (0)
                                    call r1f2kb(ido,l1,c,in,ch,1,wa(iw1))
                                case default
                                    call r1f2kb(ido,l1,ch,1,c,in,wa(iw1))
                            end select
                            na = 1-na
                        case (3)
                            iw2 = iw1+ido
                            select case (na)
                                case (0)
                                    call r1f3kb(ido,l1,c,in,ch,1,wa(iw1),wa(iw2))
                                case default
                                    call r1f3kb(ido,l1,ch,1,c,in,wa(iw1),wa(iw2))
                            end select
                            na = 1-na
                        case (4)
                            iw2 = iw1+ido
                            iw3 = iw2+ido
                            select case (na)
                                case (0)
                                    call r1f4kb(ido,l1,c,in,ch,1,wa(iw1),wa(iw2),wa(iw3))
                                case default
                                    call r1f4kb(ido,l1,ch,1,c,in,wa(iw1),wa(iw2),wa(iw3))
                            end select
                            na = 1-na
                        case (5)
                            iw2 = iw1+ido
                            iw3 = iw2+ido
                            iw4 = iw3+ido
                            select case (na)
                                case (0)
                                    call r1f5kb(ido,l1,c,in,ch,1,wa(iw1),wa(iw2),wa(iw3),wa(iw4))
                                case default
                                    call r1f5kb(ido,l1,ch,1,c,in,wa(iw1),wa(iw2),wa(iw3),wa(iw4))
                            end select
                            na = 1-na
                        case default
                            select case (na)
                                case (0)
                                    call r1fgkb(ido,iip,l1,idl1,c,c,c,in,ch,ch,1,wa(iw1))
                                case default
                                    call r1fgkb(ido,iip,l1,idl1,ch,ch,ch,1,c,c,in,wa(iw1))
                            end select
                            if (ido == 1) then
                                na = 1-na
                            end if
                    end select
                    l1 = l2
                    iw1 = iw1+(iip-1)*ido
                end do
            end associate

        end subroutine rfftb1


        subroutine r1f2kb(ido,l1,cc,in1,ch,in2,wa1)
            !--------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            real (wp),    intent (in out) :: cc(in1,ido,2,l1)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,l1,2)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa1(ido)
            !--------------------------------------------------------------
            ! Dictionary: local variables
            !--------------------------------------------------------------
            integer (ip) :: i, ic, idp2
            !--------------------------------------------------------------

            ch(1,1,:,1) = cc(1,1,1,:)+cc(1,ido,2,:)
            ch(1,1,:,2) = cc(1,1,1,:)-cc(1,ido,2,:)

            if (ido < 2) then
                return
            else
                select case (ido)
                    case (2)
                        ch(1,ido,:,1) = cc(1,ido,1,:)+cc(1,ido,1,:)
                        ch(1,ido,:,2) = -(cc(1,1,2,:)+cc(1,1,2,:))
                    case default
                        idp2 = ido+2
                        do i=3,ido,2
                            ic = idp2-i

                            ch(1,i-1,:,1) = cc(1,i-1,1,:)+cc(1,ic-1,2,:)
                            ch(1,i,:,1) = cc(1,i,1,:)-cc(1,ic,2,:)

                            ch(1,i-1,:,2) = wa1(i-2)*(cc(1,i-1,1,:)-cc(1,ic-1,2,:)) &
                                -wa1(i-1)*(cc(1,i,1,:)+cc(1,ic,2,:))
                            ch(1,i,:,2) = wa1(i-2)*(cc(1,i,1,:)+cc(1,ic,2,:))+wa1(i-1) &
                                *(cc(1,i-1,1,:)-cc(1,ic-1,2,:))
                        end do
                        if (mod(ido,2) /= 1) then
                            ch(1,ido,:,1) = cc(1,ido,1,:)+cc(1,ido,1,:)
                            ch(1,ido,:,2) = -(cc(1,1,2,:)+cc(1,1,2,:))
                        end if
                end select
            end if

        end subroutine r1f2kb


        subroutine r1f3kb(ido,l1,cc,in1,ch,in2,wa1,wa2)
            !------------------------------------------------------------------
            ! Dictionary: calling arguments
            !------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            real (wp),    intent (in out) :: cc(in1,ido,3,l1)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,l1,3)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa1(ido)
            real (wp),    intent (in)     :: wa2(ido)
            !------------------------------------------------------------------
            ! Dictionary: calling arguments
            !------------------------------------------------------------------
            integer (ip)         :: i, ic, idp2
            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            real (wp), parameter :: ARG =TWO_PI/3
            real (wp), parameter :: TAUR = cos(ARG)
            real (wp), parameter :: TAUI = sin(ARG)
            !------------------------------------------------------------------

            ch(1,1,:,1) = cc(1,1,1,:) + 2.0_wp * cc(1,ido,2,:)
            ch(1,1,:,2) = cc(1,1,1,:) + ( 2.0_wp * TAUR ) * cc(1,ido,2,:) &
                - ( 2.0_wp *TAUI)*cc(1,1,3,:)
            ch(1,1,:,3) = cc(1,1,1,:) + ( 2.0_wp *TAUR)*cc(1,ido,2,:) &
                + 2.0_wp *TAUI*cc(1,1,3,:)

            select case (ido)
                case (1)
                    return
                case default
                    idp2 = ido+2
                    do i=3,ido,2
                        ic = idp2-i
                        ch(1,i-1,:,1) = cc(1,i-1,1,:)+(cc(1,i-1,3,:)+cc(1,ic-1,2,:))
                        ch(1,i,:,1) = cc(1,i,1,:)+(cc(1,i,3,:)-cc(1,ic,2,:))

                        ch(1,i-1,:,2) = wa1(i-2)* &
                            ((cc(1,i-1,1,:)+TAUR*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))- &
                            (TAUI*(cc(1,i,3,:)+cc(1,ic,2,:)))) &
                            -wa1(i-1)* &
                            ((cc(1,i,1,:)+TAUR*(cc(1,i,3,:)-cc(1,ic,2,:)))+ &
                            (TAUI*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))))

                        ch(1,i,:,2) = wa1(i-2)* &
                            ((cc(1,i,1,:)+TAUR*(cc(1,i,3,:)-cc(1,ic,2,:)))+ &
                            (TAUI*(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))) &
                            +wa1(i-1)* &
                            ((cc(1,i-1,1,:)+TAUR*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))- &
                            (TAUI*(cc(1,i,3,:)+cc(1,ic,2,:))))

                        ch(1,i-1,:,3) = wa2(i-2)* &
                            ((cc(1,i-1,1,:)+TAUR*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))+ &
                            (TAUI*(cc(1,i,3,:)+cc(1,ic,2,:)))) &
                            -wa2(i-1)* &
                            ((cc(1,i,1,:)+TAUR*(cc(1,i,3,:)-cc(1,ic,2,:)))- &
                            (TAUI*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))))

                        ch(1,i,:,3) = wa2(i-2)* &
                            ((cc(1,i,1,:)+TAUR*(cc(1,i,3,:)-cc(1,ic,2,:)))- &
                            (TAUI*(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))) &
                            +wa2(i-1)* &
                            ((cc(1,i-1,1,:)+TAUR*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))+ &
                            (TAUI*(cc(1,i,3,:)+cc(1,ic,2,:))))
                    end do
            end select

        end subroutine r1f3kb

        subroutine r1f4kb(ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)
            !------------------------------------------------------------------
            ! Dictionary: calling arguments
            !------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            real (wp),    intent (in out) :: cc(in1,ido,4,l1)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,l1,4)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa1(ido)
            real (wp),    intent (in)     :: wa2(ido)
            real (wp),    intent (in)     :: wa3(ido)
            !------------------------------------------------------------------
            ! Dictionary: local variables
            !------------------------------------------------------------------
            integer (ip)         :: i, ic, idp2
            real (wp), parameter :: SQRT2 = sqrt(2.0_wp)
            !------------------------------------------------------------------


            ch(1,1,:,3) = (cc(1,1,1,:)+cc(1,ido,4,:)) &
                -(cc(1,ido,2,:)+cc(1,ido,2,:))
            ch(1,1,:,1) = (cc(1,1,1,:)+cc(1,ido,4,:)) &
                +(cc(1,ido,2,:)+cc(1,ido,2,:))
            ch(1,1,:,4) = (cc(1,1,1,:)-cc(1,ido,4,:)) &
                +(cc(1,1,3,:)+cc(1,1,3,:))
            ch(1,1,:,2) = (cc(1,1,1,:)-cc(1,ido,4,:)) &
                -(cc(1,1,3,:)+cc(1,1,3,:))

            if (ido < 2) then
                return
            else
                select case (ido)
                    case (2)
                        ch(1,ido,:,1) = (cc(1,ido,1,:)+cc(1,ido,3,:)) &
                            +(cc(1,ido,1,:)+cc(1,ido,3,:))
                        ch(1,ido,:,2) = SQRT2*((cc(1,ido,1,:)-cc(1,ido,3,:)) &
                            -(cc(1,1,2,:)+cc(1,1,4,:)))
                        ch(1,ido,:,3) = (cc(1,1,4,:)-cc(1,1,2,:)) &
                            +(cc(1,1,4,:)-cc(1,1,2,:))
                        ch(1,ido,:,4) = -SQRT2*((cc(1,ido,1,:)-cc(1,ido,3,:)) &
                            +(cc(1,1,2,:)+cc(1,1,4,:)))
                    case default
                        idp2 = ido+2
                        do i=3,ido,2
                            ic = idp2-i
                            ch(1,i-1,:,1) = (cc(1,i-1,1,:)+cc(1,ic-1,4,:)) &
                                +(cc(1,i-1,3,:)+cc(1,ic-1,2,:))
                            ch(1,i,:,1) = (cc(1,i,1,:)-cc(1,ic,4,:)) &
                                +(cc(1,i,3,:)-cc(1,ic,2,:))
                            ch(1,i-1,:,2)=wa1(i-2)*((cc(1,i-1,1,:)-cc(1,ic-1,4,:)) &
                                -(cc(1,i,3,:)+cc(1,ic,2,:)))-wa1(i-1) &
                                *((cc(1,i,1,:)+cc(1,ic,4,:))+(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))
                            ch(1,i,:,2)=wa1(i-2)*((cc(1,i,1,:)+cc(1,ic,4,:)) &
                                +(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))+wa1(i-1) &
                                *((cc(1,i-1,1,:)-cc(1,ic-1,4,:))-(cc(1,i,3,:)+cc(1,ic,2,:)))
                            ch(1,i-1,:,3)=wa2(i-2)*((cc(1,i-1,1,:)+cc(1,ic-1,4,:)) &
                                -(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))-wa2(i-1) &
                                *((cc(1,i,1,:)-cc(1,ic,4,:))-(cc(1,i,3,:)-cc(1,ic,2,:)))
                            ch(1,i,:,3)=wa2(i-2)*((cc(1,i,1,:)-cc(1,ic,4,:)) &
                                -(cc(1,i,3,:)-cc(1,ic,2,:)))+wa2(i-1) &
                                *((cc(1,i-1,1,:)+cc(1,ic-1,4,:))-(cc(1,i-1,3,:) &
                                +cc(1,ic-1,2,:)))
                            ch(1,i-1,:,4)=wa3(i-2)*((cc(1,i-1,1,:)-cc(1,ic-1,4,:)) &
                                +(cc(1,i,3,:)+cc(1,ic,2,:)))-wa3(i-1) &
                                *((cc(1,i,1,:)+cc(1,ic,4,:))-(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))
                            ch(1,i,:,4)=wa3(i-2)*((cc(1,i,1,:)+cc(1,ic,4,:)) &
                                -(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))+wa3(i-1) &
                                *((cc(1,i-1,1,:)-cc(1,ic-1,4,:))+(cc(1,i,3,:)+cc(1,ic,2,:)))
                        end do

                        if (mod(ido,2) /= 1) then
                            ch(1,ido,:,1) = (cc(1,ido,1,:)+cc(1,ido,3,:)) &
                                +(cc(1,ido,1,:)+cc(1,ido,3,:))
                            ch(1,ido,:,2) = SQRT2*((cc(1,ido,1,:)-cc(1,ido,3,:)) &
                                -(cc(1,1,2,:)+cc(1,1,4,:)))
                            ch(1,ido,:,3) = (cc(1,1,4,:)-cc(1,1,2,:)) &
                                +(cc(1,1,4,:)-cc(1,1,2,:))
                            ch(1,ido,:,4) = -SQRT2*((cc(1,ido,1,:)-cc(1,ido,3,:)) &
                                +(cc(1,1,2,:)+cc(1,1,4,:)))
                        end if
                end select
            end if

        end subroutine r1f4kb


        subroutine r1f5kb(ido,l1,cc,in1,ch,in2,wa1,wa2,wa3,wa4)
            !--------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            real (wp),    intent (in out) :: cc(in1,ido,5,l1)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,l1,5)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa1(ido)
            real (wp),    intent (in)     :: wa2(ido)
            real (wp),    intent (in)     :: wa3(ido)
            real (wp),    intent (in)     :: wa4(ido)
            !--------------------------------------------------
            ! Dictionary: local variables
            !--------------------------------------------------
            integer (ip)         :: i, ic, idp2
            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            real (wp), parameter :: ARG= TWO_PI/5
            real (wp), parameter :: TR11=cos(ARG)
            real (wp), parameter :: TI11=sin(ARG)
            real (wp), parameter :: TR12=cos(2.0_wp*ARG)
            real (wp), parameter :: TI12=sin(2.0_wp*ARG)
            !--------------------------------------------------

            ch(1,1,:,1) = cc(1,1,1,:)+ 2.0_wp *cc(1,ido,2,:)+ 2.0_wp *cc(1,ido,4,:)
            ch(1,1,:,2) = (cc(1,1,1,:)+TR11* 2.0_wp *cc(1,ido,2,:) &
                +TR12* 2.0_wp *cc(1,ido,4,:))-(TI11* 2.0_wp *cc(1,1,3,:) &
                +TI12* 2.0_wp *cc(1,1,5,:))
            ch(1,1,:,3) = (cc(1,1,1,:)+TR12* 2.0_wp *cc(1,ido,2,:) &
                +TR11* 2.0_wp *cc(1,ido,4,:))-(TI12* 2.0_wp *cc(1,1,3,:) &
                -TI11* 2.0_wp *cc(1,1,5,:))
            ch(1,1,:,4) = (cc(1,1,1,:)+TR12* 2.0_wp *cc(1,ido,2,:) &
                +TR11* 2.0_wp *cc(1,ido,4,:))+(TI12* 2.0_wp *cc(1,1,3,:) &
                -TI11* 2.0_wp *cc(1,1,5,:))
            ch(1,1,:,5) = (cc(1,1,1,:)+TR11* 2.0_wp *cc(1,ido,2,:) &
                +TR12* 2.0_wp *cc(1,ido,4,:))+(TI11* 2.0_wp *cc(1,1,3,:) &
                +TI12* 2.0_wp *cc(1,1,5,:))

            select case (ido)
                case (1)
                    return
                case default
                    idp2 = ido+2
                    do i=3,ido,2
                        ic = idp2-i
                        ch(1,i-1,:,1) = cc(1,i-1,1,:)+(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                            +(cc(1,i-1,5,:)+cc(1,ic-1,4,:))
                        ch(1,i,:,1) = cc(1,i,1,:)+(cc(1,i,3,:)-cc(1,ic,2,:)) &
                            +(cc(1,i,5,:)-cc(1,ic,4,:))
                        ch(1,i-1,:,2) = wa1(i-2)*((cc(1,i-1,1,:)+TR11* &
                            (cc(1,i-1,3,:)+cc(1,ic-1,2,:))+TR12 &
                            *(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))-(TI11*(cc(1,i,3,:) &
                            +cc(1,ic,2,:))+TI12*(cc(1,i,5,:)+cc(1,ic,4,:)))) &
                            -wa1(i-1)*((cc(1,i,1,:)+TR11*(cc(1,i,3,:)-cc(1,ic,2,:)) &
                            +TR12*(cc(1,i,5,:)-cc(1,ic,4,:)))+(TI11*(cc(1,i-1,3,:) &
                            -cc(1,ic-1,2,:))+TI12*(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))

                        ch(1,i,:,2) = wa1(i-2)*((cc(1,i,1,:)+TR11*(cc(1,i,3,:) &
                            -cc(1,ic,2,:))+TR12*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                            +(TI11*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))+TI12 &
                            *(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))+wa1(i-1) &
                            *((cc(1,i-1,1,:)+TR11*(cc(1,i-1,3,:) &
                            +cc(1,ic-1,2,:))+TR12*(cc(1,i-1,5,:)+cc(1,ic-1,4,:))) &
                            -(TI11*(cc(1,i,3,:)+cc(1,ic,2,:))+TI12 &
                            *(cc(1,i,5,:)+cc(1,ic,4,:))))

                        ch(1,i-1,:,3) = wa2(i-2) &
                            *((cc(1,i-1,1,:)+TR12*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                            +TR11*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))-(TI12*(cc(1,i,3,:) &
                            +cc(1,ic,2,:))-TI11*(cc(1,i,5,:)+cc(1,ic,4,:)))) &
                            -wa2(i-1) &
                            *((cc(1,i,1,:)+TR12*(cc(1,i,3,:)- &
                            cc(1,ic,2,:))+TR11*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                            +(TI12*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))-TI11 &
                            *(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))

                        ch(1,i,:,3) = wa2(i-2) &
                            *((cc(1,i,1,:)+TR12*(cc(1,i,3,:)- &
                            cc(1,ic,2,:))+TR11*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                            +(TI12*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))-TI11 &
                            *(cc(1,i-1,5,:)-cc(1,ic-1,4,:)))) &
                            +wa2(i-1) &
                            *((cc(1,i-1,1,:)+TR12*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                            +TR11*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))-(TI12*(cc(1,i,3,:) &
                            +cc(1,ic,2,:))-TI11*(cc(1,i,5,:)+cc(1,ic,4,:))))

                        ch(1,i-1,:,4) = wa3(i-2) &
                            *((cc(1,i-1,1,:)+TR12*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                            +TR11*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))+(TI12*(cc(1,i,3,:) &
                            +cc(1,ic,2,:))-TI11*(cc(1,i,5,:)+cc(1,ic,4,:)))) &
                            -wa3(i-1) &
                            *((cc(1,i,1,:)+TR12*(cc(1,i,3,:)- &
                            cc(1,ic,2,:))+TR11*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                            -(TI12*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))-TI11 &
                            *(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))

                        ch(1,i,:,4) = wa3(i-2) &
                            *((cc(1,i,1,:)+TR12*(cc(1,i,3,:)- &
                            cc(1,ic,2,:))+TR11*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                            -(TI12*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))-TI11 &
                            *(cc(1,i-1,5,:)-cc(1,ic-1,4,:)))) &
                            +wa3(i-1) &
                            *((cc(1,i-1,1,:)+TR12*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                            +TR11*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))+(TI12*(cc(1,i,3,:) &
                            +cc(1,ic,2,:))-TI11*(cc(1,i,5,:)+cc(1,ic,4,:))))

                        ch(1,i-1,:,5) = wa4(i-2) &
                            *((cc(1,i-1,1,:)+TR11*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                            +TR12*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))+(TI11*(cc(1,i,3,:) &
                            +cc(1,ic,2,:))+TI12*(cc(1,i,5,:)+cc(1,ic,4,:)))) &
                            -wa4(i-1) &
                            *((cc(1,i,1,:)+TR11*(cc(1,i,3,:)-cc(1,ic,2,:)) &
                            +TR12*(cc(1,i,5,:)-cc(1,ic,4,:)))-(TI11*(cc(1,i-1,3,:) &
                            -cc(1,ic-1,2,:))+TI12*(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))

                        ch(1,i,:,5) = wa4(i-2) &
                            *((cc(1,i,1,:)+TR11*(cc(1,i,3,:)-cc(1,ic,2,:)) &
                            +TR12*(cc(1,i,5,:)-cc(1,ic,4,:)))-(TI11*(cc(1,i-1,3,:) &
                            -cc(1,ic-1,2,:))+TI12*(cc(1,i-1,5,:)-cc(1,ic-1,4,:)))) &
                            +wa4(i-1) &
                            *((cc(1,i-1,1,:)+TR11*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                            +TR12*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))+(TI11*(cc(1,i,3,:) &
                            +cc(1,ic,2,:))+TI12*(cc(1,i,5,:)+cc(1,ic,4,:))))
                    end do
            end select

        end subroutine r1f5kb


        subroutine r1fgkb(ido,iip,l1,idl1,cc,c1,c2,in1,ch,ch2,in2,wa)
            !------------------------------------------------------------------
            ! Dictionary: calling arguments
            !------------------------------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: iip
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: idl1
            real (wp),    intent (in out) :: c1(in1,ido,l1,iip)
            real (wp),    intent (in out) :: c2(in1,idl1,iip)
            real (wp),    intent (in out) :: cc(in1,ido,iip,l1)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,l1,iip)
            real (wp),    intent (in out) :: ch2(in2,idl1,iip)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido)
            !------------------------------------------------------------------
            ! Dictionary: local variables
            !------------------------------------------------------------------
            real (wp)            :: ai1, ai2, ar1, ar1h, ar2, ar2h, arg
            real (wp)            :: dc2, dcp, ds2, dsp
            integer (ip)         :: i, idij, idp2, ipp2, ipph
            integer (ip)         :: is, j, j2, jc, k, l, lc, nbd
            real (wp), parameter :: TWO_PI = acos(-1.0_wp)
            !------------------------------------------------------------------

            arg = TWO_PI /iip
            dcp = cos(arg)
            dsp = sin(arg)
            idp2 = ido+2
            nbd = (ido-1)/2
            ipp2 = iip+2
            ipph = (iip+1)/2

            ch(1,:,:,1) = cc(1,:,1,:)

            do j=2,ipph
                jc = ipp2-j
                j2 = j*2
                ch(1,1,:,j) = 2.0_wp * cc(1,ido,j2-2,:)
                ch(1,1,:,jc) = 2.0_wp * cc(1,1,j2-1,:)
            end do

            if (ido /= 1) then
                if (nbd >= l1) then
                    do j=2,ipph
                        jc = ipp2-j
                        ch(1,2:ido-1:2,:, j) = cc(1,2:ido-1:2, 2*j-1,:) + cc(1,idp2-4: &
                            idp2-1-ido:(-2), 2*j-2,:)
                        ch(1,2:ido-1:2,:, jc) = cc(1,2:ido-1:2, 2*j-1,:) - cc(1,idp2-4: &
                            idp2-1-ido:(-2), 2*j-2,:)
                        ch(1,3:ido:2,:, j) = cc(1,3:ido:2, 2*j-1,:) - cc(1,idp2-3:idp2- &
                            ido:(-2), 2*j-2,:)
                        ch(1,3:ido:2,:, jc) = cc(1,3:ido:2, 2*j-1,:) + cc(1,idp2-3:idp2- &
                            ido:(-2), 2*j-2,:)
                    end do
                else
                    do j=2,ipph
                        jc = ipp2-j
                        ch(1,2:ido-1:2,:, j) = cc(1,2:ido-1:2, 2*j-1,:) + cc(1,idp2-4: &
                            idp2-1-ido:(-2), 2*j-2,:)
                        ch(1,2:ido-1:2,:, jc) = cc(1,2:ido-1:2, 2*j-1,:) - cc(1,idp2-4: &
                            idp2-1-ido:(-2), 2*j-2,:)
                        ch(1,3:ido:2,:, j) = cc(1,3:ido:2, 2*j-1,:) - cc(1,idp2-3:idp2- &
                            ido:(-2), 2*j-2,:)
                        ch(1,3:ido:2,:, jc) = cc(1,3:ido:2, 2*j-1,:) + cc(1,idp2-3:idp2- &
                            ido:(-2), 2*j-2,:)
                    end do
                end if
            end if

            ar1 = 1.0_wp
            ai1 = 0.0_wp

            do l=2,ipph
                lc = ipp2-l
                ar1h = dcp*ar1-dsp*ai1
                ai1 = dcp*ai1+dsp*ar1
                ar1 = ar1h
                c2(1,:,l) = ch2(1,:,1)+ar1*ch2(1,:,2)
                c2(1,:,lc) = ai1*ch2(1,:,iip)
                dc2 = ar1
                ds2 = ai1
                ar2 = ar1
                ai2 = ai1

                do j=3,ipph
                    jc = ipp2-j
                    ar2h = dc2*ar2-ds2*ai2
                    ai2 = dc2*ai2+ds2*ar2
                    ar2 = ar2h
                    c2(1,:,l) = c2(1,:,l)+ar2*ch2(1,:,j)
                    c2(1,:,lc) = c2(1,:,lc)+ai2*ch2(1,:,jc)
                end do
            end do

            do j=2,ipph
                ch2(1,:,1) = ch2(1,:,1)+ch2(1,:,j)
            end do

            do j=2,ipph
                jc = ipp2-j
                ch(1,1,:,j) = c1(1,1,:,j)-c1(1,1,:,jc)
                ch(1,1,:,jc) = c1(1,1,:,j)+c1(1,1,:,jc)
            end do

            if (ido /= 1) then
                if (nbd >= l1) then
                    do j=2,ipph
                        jc = ipp2-j
                        ch(1,2:ido-1:2,:, j) = c1(1,2:ido-1:2,:, j) - c1(1,3:ido:2,:, jc)
                        ch(1,2:ido-1:2,:, jc) = c1(1,2:ido-1:2,:, j) + c1(1,3:ido:2,:, jc)
                        ch(1,3:ido:2,:, j) = c1(1,3:ido:2,:, j) + c1(1,2:ido-1:2,:, jc)
                        ch(1,3:ido:2,:, jc) = c1(1,3:ido:2,:, j) - c1(1,2:ido-1:2,:, jc)
                    end do
                else
                    do j=2,ipph
                        jc = ipp2-j
                        ch(1,2:ido-1:2,:, j) = c1(1,2:ido-1:2,:, j) - c1(1,3:ido:2,:, jc)
                        ch(1,2:ido-1:2,:, jc) = c1(1,2:ido-1:2,:, j) + c1(1,3:ido:2,:, jc)
                        ch(1,3:ido:2,:, j) = c1(1,3:ido:2,:, j) + c1(1,2:ido-1:2,:, jc)
                        ch(1,3:ido:2,:, jc) = c1(1,3:ido:2,:, j) - c1(1,2:ido-1:2,:, jc)
                    end do
                end if
            end if

            select case (ido)
                case (1)
                    return
                case default

                    c2(1,:,1) = ch2(1,:,1)
                    c1(1,1,:,2:iip) = ch(1,1,:,2:iip)

                    if ( l1 >= nbd ) then
                        is = -ido
                        do j=2,iip
                            is = is+ido
                            idij = is
                            do i=3,ido,2
                                idij = idij+2
                                c1(1,i-1,:,j) = &
                                    wa(idij-1)*ch(1,i-1,:,j) &
                                    -wa(idij)*ch(1,i,:,j)
                                c1(1,i,:,j) = &
                                    wa(idij-1)*ch(1,i,:,j) &
                                    +wa(idij)*ch(1,i-1,:,j)
                            end do
                        end do
                    else
                        is = -ido
                        do j=2,iip
                            is = is+ido
                            do k=1,l1
                                idij = is
                                c1(1,2:ido-1:2, k, j) = &
                                    wa(idij+1:ido-2+idij:2)*ch(1,2:ido-1:2, &
                                    k, j) - wa(idij+2:ido-1+idij:2)*ch(1,3:ido:2, k, j)
                                c1(1,3:ido:2, k, j) = &
                                    wa(idij+1:ido-2+idij:2)*ch(1,3:ido:2, k, j) &
                                    + wa(idij+2:ido-1+idij:2)*ch(1,2:ido-1:2, k, j)
                            end do
                        end do
                    end if
            end select

        end subroutine r1fgkb



    end subroutine rfft1b


    subroutine rfft1f(n, inc, r, lenr, wsave, lensav, work, lenwrk, ier)
        !
        ! rfft1f: 1d forward fast fourier transform
        !
        !  Purpose:
        !
        !  rfft1f computes the one-dimensional fourier transform of a periodic
        !  sequence within a real array. This is referred to as the forward
        !  transform or fourier analysis, transforming the sequence from physical
        !  to spectral space. This transform is normalized since a call to
        !  rfft1f followed by a call to rfft1b (or vice-versa) reproduces the
        !  original array within roundoff error.
        !
        !  Parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr), on input, contains the sequence
        !  to be transformed, and on output, the transformed data.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1) + 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to rfft1i before the first call to routine rfft1f
        !  or rfft1b for a given transform length n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least n.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough:
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough.
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: inc
        real (wp),    intent (in out) :: r(lenr)
        integer (ip), intent (in)     :: lenr
        real (wp),    intent (out)    :: work(lenwrk)
        integer (ip), intent (in)     :: lenwrk
        real (wp),    intent (in)     :: wsave(lensav)
        integer (ip), intent (in)     :: lensav
        integer (ip), intent (out)    :: ier
        !--------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (lenr < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('rfft1f ', 6)
        else if (lensav < n + int(log(real(n, kind=wp) )/log(2.0_wp)) +4) then
            ier = 2
            call fft_error_handler('rfft1f ', 8)
        else if (lenwrk < n) then
            ier = 3
            call fft_error_handler('rfft1f ', 10)
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            call rfftf1(n,inc,r,work,wsave,wsave(n+1))
        end if


    contains


        subroutine rfftf1(n, in, c, ch, wa, fac)
            !--------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: in
            real (wp),    intent (in out) :: c(in,*)
            real (wp),    intent (out)    :: ch(:)
            real (wp),    intent (in)     :: wa(n)
            real (wp),    intent (in)     :: fac(15)
            !--------------------------------------------------------------
            ! Dictionary: calling arguments
            !--------------------------------------------------------------
            integer (ip) :: idl1, ido, iip
            integer (ip) :: workspace_indices(4)
            integer (ip) :: j,  k1, kh, l1, l2
            integer (ip) :: modn, na, nf, nl
            real (wp)    :: sn, tsn, tsnm
            !--------------------------------------------------------------

            nf = int(fac(2), kind=ip)
            na = 1
            l2 = n

            associate( &
                iw1 => workspace_indices(1), &
                iw2 => workspace_indices(2), &
                iw3 => workspace_indices(3), &
                iw4 => workspace_indices(4) &
                )

                iw1 = n

                do k1=1,nf
                    kh = nf-k1
                    iip = int(fac(kh+3), kind=ip)
                    l1 = l2/iip
                    ido = n/l2
                    idl1 = ido*l1
                    iw1 = iw1-(iip-1)*ido
                    na = 1-na
                    select case (iip)
                        case (2)
                            select case (na)
                                case (0)
                                    call r1f2kf(ido,l1,c,in,ch,1,wa(iw1))
                                case default
                                    call r1f2kf(ido,l1,ch,1,c,in,wa(iw1))
                            end select
                        case (3)
                            iw2 = iw1+ido
                            select case (na)
                                case (0)
                                    call r1f3kf(ido,l1,c,in,ch,1,wa(iw1),wa(iw2))
                                case default
                                    call r1f3kf(ido,l1,ch,1,c,in,wa(iw1),wa(iw2))
                            end select
                        case (4)
                            iw2 = iw1+ido
                            iw3 = iw2+ido
                            select case (na)
                                case (0)
                                    call r1f4kf(ido,l1,c,in,ch,1,wa(iw1),wa(iw2),wa(iw3))
                                case default
                                    call r1f4kf(ido,l1,ch,1,c,in,wa(iw1),wa(iw2),wa(iw3))
                            end select
                        case (5)
                            iw2 = iw1+ido
                            iw3 = iw2+ido
                            iw4 = iw3+ido
                            select case (na)
                                case (0)
                                    call r1f5kf(ido,l1,c,in,ch,1,wa(iw1),wa(iw2),wa(iw3),wa(iw4))
                                case default
                                    call r1f5kf(ido,l1,ch,1,c,in,wa(iw1),wa(iw2),wa(iw3),wa(iw4))
                            end select
                        case default
                            if (ido == 1) then
                                na = 1-na
                            end if
                            select case (na)
                                case (0)
                                    call r1fgkf(ido,iip,l1,idl1,c,c,c,in,ch,ch,1,wa(iw1))
                                    na = 1
                                case default
                                    call r1fgkf(ido,iip,l1,idl1,ch,ch,ch,1,c,c,in,wa(iw1))
                                    na = 0
                            end select
                    end select
                    l2 = l1
                end do

            end associate

            sn = 1.0_wp/n
            tsn =  2.0_wp/n
            tsnm = -tsn
            modn = mod(n,2)

            if (modn /= 0) then
                nl = n-1
            else
                nl = n-2
            end if

            if (na == 0) then
                c(1,1) = sn*ch(1)
                do j=2,nl,2
                    c(1,j) = tsn*ch(j)
                    c(1,j+1) = tsnm*ch(j+1)
                end do
                if (modn == 0) then
                    c(1,n) = sn*ch(n)
                end if
            else
                c(1,1) = sn*c(1,1)
                do j=2,nl,2
                    c(1,j) = tsn*c(1,j)
                    c(1,j+1) = tsnm*c(1,j+1)
                end do
                if (modn == 0) then
                    c(1,n) = sn*c(1,n)
                end if
            end if

        end subroutine rfftf1


        subroutine r1f2kf(ido,l1,cc,in1,ch,in2,wa1)
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            real (wp),    intent (in out) :: cc(in1,ido,l1,2)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,2,l1)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa1(ido)
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip) :: i, ic, idp2
            !-----------------------------------------------

            ch(1,1,1,:) = cc(1,1,:,1)+cc(1,1,:,2)
            ch(1,ido,2,:) = cc(1,1,:,1)-cc(1,1,:,2)

            if (ido < 2) then
                return
            end if

            select case (ido)
                case (2)
                    ch(1,1,2,:) = -cc(1,ido,:,2)
                    ch(1,ido,1,:) = cc(1,ido,:,1)
                case default
                    idp2 = ido+2
                    do i=3,ido,2
                        ic = idp2-i

                        ch(1,i,1,:) = &
                            cc(1,i,:,1)+(wa1(i-2)*cc(1,i,:,2) &
                            -wa1(i-1)*cc(1,i-1,:,2))

                        ch(1,ic,2,:) = &
                            (wa1(i-2)*cc(1,i,:,2) &
                            -wa1(i-1)*cc(1,i-1,:,2))-cc(1,i,:,1)

                        ch(1,i-1,1,:) = &
                            cc(1,i-1,:,1)+(wa1(i-2)*cc(1,i-1,:,2) &
                            +wa1(i-1)*cc(1,i,:,2))

                        ch(1,ic-1,2,:) = &
                            cc(1,i-1,:,1)-(wa1(i-2)*cc(1,i-1,:,2) &
                            +wa1(i-1)*cc(1,i,:,2))
                    end do

                    if (mod(ido,2) /= 1) then
                        ch(1,1,2,:) = -cc(1,ido,:,2)
                        ch(1,ido,1,:) = cc(1,ido,:,1)
                    end if
            end select

        end subroutine r1f2kf


        subroutine r1f3kf(ido,l1,cc,in1,ch,in2,wa1,wa2)
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            real (wp),    intent (in out) :: cc(in1,ido,l1,3)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,3,l1)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa1(ido)
            real (wp),    intent (in)     :: wa2(ido)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip)         :: i, ic, idp2
            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            real (wp), parameter :: ARG = TWO_PI/3
            real (wp), parameter :: TAUR = cos(ARG)
            real (wp), parameter :: TAUI = sin(ARG)
            !-----------------------------------------------

            ch(1,1,1,:) = cc(1,1,:,1)+(cc(1,1,:,2)+cc(1,1,:,3))
            ch(1,1,3,:) = TAUI*(cc(1,1,:,3)-cc(1,1,:,2))
            ch(1,ido,2,:) = cc(1,1,:,1)+TAUR*(cc(1,1,:,2)+cc(1,1,:,3))

            select case (ido)
                case (1)
                    return
                case default
                    idp2 = ido+2
                    do i=3,ido,2
                        ic = idp2-i

                        ch(1,i-1,1,:) = cc(1,i-1,:,1)+((wa1(i-2)*cc(1,i-1,:,2)+ &
                            wa1(i-1)*cc(1,i,:,2))+(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                            cc(1,i,:,3)))

                        ch(1,i,1,:) = cc(1,i,:,1)+((wa1(i-2)*cc(1,i,:,2)- &
                            wa1(i-1)*cc(1,i-1,:,2))+(wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                            cc(1,i-1,:,3)))

                        ch(1,i-1,3,:) = (cc(1,i-1,:,1)+TAUR*((wa1(i-2)* &
                            cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2))+(wa2(i-2)* &
                            cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3))))+(TAUI*((wa1(i-2)* &
                            cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2))-(wa2(i-2)* &
                            cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3))))

                        ch(1,ic-1,2,:) = (cc(1,i-1,:,1)+TAUR*((wa1(i-2)* &
                            cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2))+(wa2(i-2)* &
                            cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3))))-(TAUI*((wa1(i-2)* &
                            cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2))-(wa2(i-2)* &
                            cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3))))

                        ch(1,i,3,:) = (cc(1,i,:,1)+TAUR*((wa1(i-2)*cc(1,i,:,2)- &
                            wa1(i-1)*cc(1,i-1,:,2))+(wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                            cc(1,i-1,:,3))))+(TAUI*((wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                            cc(1,i,:,3))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2))))

                        ch(1,ic,2,:) = (TAUI*((wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                            cc(1,i,:,3))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2))))-(cc(1,i,:,1)+TAUR*((wa1(i-2)*cc(1,i,:,2)- &
                            wa1(i-1)*cc(1,i-1,:,2))+(wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                            cc(1,i-1,:,3))))
                    end do
            end select

        end subroutine r1f3kf

        subroutine r1f4kf(ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            real (wp),    intent (in out) :: cc(in1,ido,l1,4)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,4,l1)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa1(ido)
            real (wp),    intent (in)     :: wa2(ido)
            real (wp),    intent (in)     :: wa3(ido)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip)         :: i, ic, idp2
            real (wp), parameter :: HALF_SQRT2=sqrt(2.0_wp)/2
            !-----------------------------------------------

            ch(1,1,1,:) = (cc(1,1,:,2)+cc(1,1,:,4))+(cc(1,1,:,1)+cc(1,1,:,3))
            ch(1,ido,4,:) = (cc(1,1,:,1)+cc(1,1,:,3))-(cc(1,1,:,2)+cc(1,1,:,4))
            ch(1,ido,2,:) = cc(1,1,:,1)-cc(1,1,:,3)
            ch(1,1,3,:) = cc(1,1,:,4)-cc(1,1,:,2)

            if (ido < 2) then
                return
            end if

            select case (ido)
                case (2)
                    ch(1,ido,1,:) = (HALF_SQRT2*(cc(1,ido,:,2)-cc(1,ido,:,4)))+cc(1,ido,:,1)
                    ch(1,ido,3,:) = cc(1,ido,:,1)-(HALF_SQRT2*(cc(1,ido,:,2)-cc(1,ido,:,4)))
                    ch(1,1,2,:) = (-HALF_SQRT2*(cc(1,ido,:,2)+cc(1,ido,:,4)))-cc(1,ido,:,3)
                    ch(1,1,4,:) = (-HALF_SQRT2*(cc(1,ido,:,2)+cc(1,ido,:,4)))+cc(1,ido,:,3)
                case default
                    idp2 = ido+2
                    do i=3,ido,2
                        ic = idp2-i
                        ch(1,i-1,1,:) = ((wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2))+(wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                            cc(1,i,:,4)))+(cc(1,i-1,:,1)+(wa2(i-2)*cc(1,i-1,:,3)+ &
                            wa2(i-1)*cc(1,i,:,3)))
                        ch(1,ic-1,4,:) = (cc(1,i-1,:,1)+(wa2(i-2)*cc(1,i-1,:,3)+ &
                            wa2(i-1)*cc(1,i,:,3)))-((wa1(i-2)*cc(1,i-1,:,2)+ &
                            wa1(i-1)*cc(1,i,:,2))+(wa3(i-2)*cc(1,i-1,:,4)+ &
                            wa3(i-1)*cc(1,i,:,4)))
                        ch(1,i,1,:) = ((wa1(i-2)*cc(1,i,:,2)-wa1(i-1)* &
                            cc(1,i-1,:,2))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4)))+(cc(1,i,:,1)+(wa2(i-2)*cc(1,i,:,3)- &
                            wa2(i-1)*cc(1,i-1,:,3)))
                        ch(1,ic,4,:) = ((wa1(i-2)*cc(1,i,:,2)-wa1(i-1)* &
                            cc(1,i-1,:,2))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4)))-(cc(1,i,:,1)+(wa2(i-2)*cc(1,i,:,3)- &
                            wa2(i-1)*cc(1,i-1,:,3)))
                        ch(1,i-1,3,:) = ((wa1(i-2)*cc(1,i,:,2)-wa1(i-1)* &
                            cc(1,i-1,:,2))-(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4)))+(cc(1,i-1,:,1)-(wa2(i-2)*cc(1,i-1,:,3)+ &
                            wa2(i-1)*cc(1,i,:,3)))
                        ch(1,ic-1,2,:) = (cc(1,i-1,:,1)-(wa2(i-2)*cc(1,i-1,:,3)+ &
                            wa2(i-1)*cc(1,i,:,3)))-((wa1(i-2)*cc(1,i,:,2)-wa1(i-1)* &
                            cc(1,i-1,:,2))-(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4)))
                        ch(1,i,3,:) = ((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                            cc(1,i,:,4))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2)))+(cc(1,i,:,1)-(wa2(i-2)*cc(1,i,:,3)- &
                            wa2(i-1)*cc(1,i-1,:,3)))
                        ch(1,ic,2,:) = ((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                            cc(1,i,:,4))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2)))-(cc(1,i,:,1)-(wa2(i-2)*cc(1,i,:,3)- &
                            wa2(i-1)*cc(1,i-1,:,3)))
                    end do

                    if (mod(ido,2) /= 1) then
                        ch(1,ido,1,:) = (HALF_SQRT2*(cc(1,ido,:,2)-cc(1,ido,:,4)))+cc(1,ido,:,1)
                        ch(1,ido,3,:) = cc(1,ido,:,1)-(HALF_SQRT2*(cc(1,ido,:,2)-cc(1,ido,:,4)))
                        ch(1,1,2,:) = (-HALF_SQRT2*(cc(1,ido,:,2)+cc(1,ido,:,4)))-cc(1,ido,:,3)
                        ch(1,1,4,:) = (-HALF_SQRT2*(cc(1,ido,:,2)+cc(1,ido,:,4)))+cc(1,ido,:,3)
                    end if
            end select

        end subroutine r1f4kf



        subroutine r1f5kf(ido,l1,cc,in1,ch,in2,wa1,wa2,wa3,wa4)
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: l1
            real (wp),    intent (in out) :: cc(in1,ido,l1,5)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,5,l1)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa1(ido)
            real (wp),    intent (in)     :: wa2(ido)
            real (wp),    intent (in)     :: wa3(ido)
            real (wp),    intent (in)     :: wa4(ido)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            integer (ip)         :: i, ic, idp2
            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            real (wp), parameter :: ARG= TWO_PI/5
            real (wp), parameter :: TR11=cos(ARG)
            real (wp), parameter :: TI11=sin(ARG)
            real (wp), parameter :: TR12=cos(2.0_wp *ARG)
            real (wp), parameter :: TI12=sin(2.0_wp *ARG)
            !-----------------------------------------------

            ch(1,1,1,:) = cc(1,1,:,1)+(cc(1,1,:,5)+cc(1,1,:,2))+ &
                (cc(1,1,:,4)+cc(1,1,:,3))
            ch(1,ido,2,:) = cc(1,1,:,1)+TR11*(cc(1,1,:,5)+cc(1,1,:,2))+ &
                TR12*(cc(1,1,:,4)+cc(1,1,:,3))
            ch(1,1,3,:) = TI11*(cc(1,1,:,5)-cc(1,1,:,2))+TI12* &
                (cc(1,1,:,4)-cc(1,1,:,3))
            ch(1,ido,4,:) = cc(1,1,:,1)+TR12*(cc(1,1,:,5)+cc(1,1,:,2))+ &
                TR11*(cc(1,1,:,4)+cc(1,1,:,3))
            ch(1,1,5,:) = TI12*(cc(1,1,:,5)-cc(1,1,:,2))-TI11* &
                (cc(1,1,:,4)-cc(1,1,:,3))

            select case (ido)
                case (1)
                    return
                case default
                    idp2 = ido+2
                    do i=3,ido,2
                        ic = idp2-i
                        ch(1,i-1,1,:) = cc(1,i-1,:,1)+((wa1(i-2)*cc(1,i-1,:,2)+ &
                            wa1(i-1)*cc(1,i,:,2))+(wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)* &
                            cc(1,i,:,5)))+((wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                            cc(1,i,:,3))+(wa3(i-2)*cc(1,i-1,:,4)+ &
                            wa3(i-1)*cc(1,i,:,4)))

                        ch(1,i,1,:) = cc(1,i,:,1)+((wa1(i-2)*cc(1,i,:,2)- &
                            wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                            cc(1,i-1,:,5)))+((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                            cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4)))

                        ch(1,i-1,3,:) = cc(1,i-1,:,1)+TR11* &
                            ( wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2) &
                            +wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)*cc(1,i,:,5))+TR12* &
                            ( wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3) &
                            +wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)*cc(1,i,:,4))+TI11* &
                            ( wa1(i-2)*cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2) &
                            -(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)*cc(1,i-1,:,5)))+TI12* &
                            ( wa2(i-2)*cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3) &
                            -(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)*cc(1,i-1,:,4)))

                        ch(1,ic-1,2,:) = cc(1,i-1,:,1)+TR11* &
                            ( wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2) &
                            +wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)*cc(1,i,:,5))+TR12* &
                            ( wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3) &
                            +wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)*cc(1,i,:,4))-(TI11* &
                            ( wa1(i-2)*cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2) &
                            -(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)*cc(1,i-1,:,5)))+TI12* &
                            ( wa2(i-2)*cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3) &
                            -(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)*cc(1,i-1,:,4))))

                        ch(1,i,3,:) = (cc(1,i,:,1)+TR11*((wa1(i-2)*cc(1,i,:,2)- &
                            wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                            cc(1,i-1,:,5)))+TR12*((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                            cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4))))+(TI11*((wa4(i-2)*cc(1,i-1,:,5)+ &
                            wa4(i-1)*cc(1,i,:,5))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2)))+TI12*((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                            cc(1,i,:,4))-(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                            cc(1,i,:,3))))

                        ch(1,ic,2,:) = (TI11*((wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)* &
                            cc(1,i,:,5))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2)))+TI12*((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                            cc(1,i,:,4))-(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                            cc(1,i,:,3))))-(cc(1,i,:,1)+TR11*((wa1(i-2)*cc(1,i,:,2)- &
                            wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                            cc(1,i-1,:,5)))+TR12*((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                            cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4))))

                        ch(1,i-1,5,:) = (cc(1,i-1,:,1)+TR12*((wa1(i-2)* &
                            cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2))+(wa4(i-2)* &
                            cc(1,i-1,:,5)+wa4(i-1)*cc(1,i,:,5)))+TR11*((wa2(i-2)* &
                            cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3))+(wa3(i-2)* &
                            cc(1,i-1,:,4)+wa3(i-1)*cc(1,i,:,4))))+(TI12*((wa1(i-2)* &
                            cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2))-(wa4(i-2)* &
                            cc(1,i,:,5)-wa4(i-1)*cc(1,i-1,:,5)))-TI11*((wa2(i-2)* &
                            cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3))-(wa3(i-2)* &
                            cc(1,i,:,4)-wa3(i-1)*cc(1,i-1,:,4))))

                        ch(1,ic-1,4,:) = (cc(1,i-1,:,1)+TR12*((wa1(i-2)* &
                            cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2))+(wa4(i-2)* &
                            cc(1,i-1,:,5)+wa4(i-1)*cc(1,i,:,5)))+TR11*((wa2(i-2)* &
                            cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3))+(wa3(i-2)* &
                            cc(1,i-1,:,4)+wa3(i-1)*cc(1,i,:,4))))-(TI12*((wa1(i-2)* &
                            cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2))-(wa4(i-2)* &
                            cc(1,i,:,5)-wa4(i-1)*cc(1,i-1,:,5)))-TI11*((wa2(i-2)* &
                            cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3))-(wa3(i-2)* &
                            cc(1,i,:,4)-wa3(i-1)*cc(1,i-1,:,4))))

                        ch(1,i,5,:) = (cc(1,i,:,1)+TR12*((wa1(i-2)*cc(1,i,:,2)- &
                            wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                            cc(1,i-1,:,5)))+TR11*((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                            cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4))))+(TI12*((wa4(i-2)*cc(1,i-1,:,5)+ &
                            wa4(i-1)*cc(1,i,:,5))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2)))-TI11*((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                            cc(1,i,:,4))-(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                            cc(1,i,:,3))))

                        ch(1,ic,4,:) = (TI12*((wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)* &
                            cc(1,i,:,5))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                            cc(1,i,:,2)))-TI11*((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                            cc(1,i,:,4))-(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                            cc(1,i,:,3))))-(cc(1,i,:,1)+TR12*((wa1(i-2)*cc(1,i,:,2)- &
                            wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                            cc(1,i-1,:,5)))+TR11*((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                            cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                            cc(1,i-1,:,4))))
                    end do
            end select

        end subroutine r1f5kf


        subroutine r1fgkf(ido,iip,l1,idl1,cc,c1,c2,in1,ch,ch2,in2,wa)
            !-----------------------------------------------
            ! Dictionary: calling arguments
            !-----------------------------------------------
            integer (ip), intent (in)     :: ido
            integer (ip), intent (in)     :: iip
            integer (ip), intent (in)     :: l1
            integer (ip), intent (in)     :: idl1
            real (wp),    intent (out)    :: cc(in1,ido,iip,l1)
            real (wp),    intent (in out) :: c1(in1,ido,l1,iip)
            real (wp),    intent (in out) :: c2(in1,idl1,iip)
            integer (ip), intent (in)     :: in1
            real (wp),    intent (in out) :: ch(in2,ido,l1,iip)
            real (wp),    intent (in out) :: ch2(in2,idl1,iip)
            integer (ip), intent (in)     :: in2
            real (wp),    intent (in)     :: wa(ido)
            !-----------------------------------------------
            ! Dictionary: local variables
            !-----------------------------------------------
            real (wp)            :: ai1, ai2, ar1, ar1h
            real (wp)            :: ar2, ar2h, arg
            real (wp)            :: dc2, dcp, ds2, dsp
            integer (ip)         :: i, idij, idp2
            integer (ip)         :: ipp2, ipph, is
            integer (ip)         :: j, jc
            integer (ip)         :: k, l, lc, nbd
            real (wp), parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
            !-----------------------------------------------

            arg = TWO_PI/iip
            dcp = cos(arg)
            dsp = sin(arg)
            ipph = (iip+1)/2
            ipp2 = iip+2
            idp2 = ido+2
            nbd = (ido-1)/2


            select case (ido)
                case (1)
                    c2(1,:,1) = ch2(1,:,1)
                case default
                    ch2(1,:,1) = c2(1,:,1)
                    ch(1,1,:,2:iip) = c1(1,1,:,2:iip)
                    if ( nbd <= l1 ) then
                        is = -ido
                        do j=2,iip
                            is = is+ido
                            idij = is
                            do i=3,ido,2
                                idij = idij+2
                                ch(1,i-1,:,j) = wa(idij-1)*c1(1,i-1,:,j)+wa(idij)*c1(1,i,:,j)
                                ch(1,i,:,j) = wa(idij-1)*c1(1,i,:,j)-wa(idij)*c1(1,i-1,:,j)
                            end do
                        end do
                    else
                        is = -ido
                        do j=2,iip
                            is = is+ido
                            do k=1,l1
                                idij = is
                                ch(1,2:ido-1:2,k,j) = &
                                    wa(idij+1:ido-2+idij:2)*c1(1,2:ido-1:2,k,j)&
                                    +wa(idij+2:ido-1+idij:2)*c1(1,3:ido:2,k,j)
                                ch(1,3:ido:2,k,j) = &
                                    wa(idij+1:ido-2+idij:2)*c1(1,3:ido:2,k,j)&
                                    -wa(idij+2:ido-1+idij:2)*c1(1,2:ido-1:2,k,j)
                            end do
                        end do
                    end if

                    if (nbd >= l1) then
                        do j=2,ipph
                            jc = ipp2-j
                            c1(1,2:ido-1:2,:, j)=ch(1,2:ido-1:2,:, j)+ch(1,2:ido-1:2,:, jc)
                            c1(1,2:ido-1:2,:, jc) = ch(1,3:ido:2,:, j) - ch(1,3:ido:2,:, jc)
                            c1(1,3:ido:2,:, j) = ch(1,3:ido:2,:, j) + ch(1,3:ido:2,:, jc)
                            c1(1,3:ido:2,:, jc) = ch(1,2:ido-1:2,:, jc) - ch(1,2:ido-1:2,:, j)
                        end do
                    else
                        do j=2,ipph
                            jc = ipp2-j
                            c1(1,2:ido-1:2,:, j) = ch(1,2:ido-1:2,:, j) + ch(1,2:ido-1:2,:, jc)
                            c1(1,2:ido-1:2,:, jc) = ch(1,3:ido:2,:, j) - ch(1,3:ido:2,:, jc)
                            c1(1,3:ido:2,:, j) = ch(1,3:ido:2,:, j) + ch(1,3:ido:2,:, jc)
                            c1(1,3:ido:2,:, jc) = ch(1,2:ido-1:2,:, jc) - ch(1,2:ido-1:2,:, j)
                        end do
                    end if
            end select

            do j=2,ipph
                jc = ipp2-j
                c1(1,1,:,j) = ch(1,1,:,j)+ch(1,1,:,jc)
                c1(1,1,:,jc) = ch(1,1,:,jc)-ch(1,1,:,j)
            end do

            ar1 = 1.0_wp
            ai1 = 0.0_wp
            do l=2,ipph
                lc = ipp2-l
                ar1h = dcp*ar1-dsp*ai1
                ai1 = dcp*ai1+dsp*ar1
                ar1 = ar1h
                ch2(1,:,l) = c2(1,:,1)+ar1*c2(1,:,2)
                ch2(1,:,lc) = ai1*c2(1,:,iip)
                dc2 = ar1
                ds2 = ai1
                ar2 = ar1
                ai2 = ai1
                do j=3,ipph
                    jc = ipp2-j
                    ar2h = dc2*ar2-ds2*ai2
                    ai2 = dc2*ai2+ds2*ar2
                    ar2 = ar2h
                    ch2(1,:,l) = ch2(1,:,l)+ar2*c2(1,:,j)
                    ch2(1,:,lc) = ch2(1,:,lc)+ai2*c2(1,:,jc)
                end do
            end do

            do j=2,ipph
                ch2(1,:,1) = ch2(1,:,1)+c2(1,:,j)
            end do

            cc(1,:,1,:) = ch(1,:,:,1)

            cc(1,ido,2:(ipph-1)*2:2,:) = transpose(ch(1,1,:,2:ipph))

            cc(1,1,3:ipph*2-1:2,:) = transpose(ch(1,1,:,ipp2-2:ipp2-ipph:(-1)))

            select case (ido)
                case (1)
                    return
                case default
                    if (nbd >= l1) then

                        cc(1,2:ido-1:2, 3:ipph*2-1:2,:) = &
                            reshape(source = &
                            ch(1,2:ido-1:2,:,2:ipph)+ch(1,2:ido-1:2,:,ipp2-2:ipp2-ipph:(-1)), &
                            shape = [(ido -1)/2, ipph-1, l1], &
                            order = [1, 3, 2])

                        cc(1,idp2-4:idp2-1-ido:(-2), 2:(ipph-1)*2:2,:) = &
                            reshape(source = &
                            ch(1,2:ido-1:2,:, 2:ipph)-ch(1,2:ido-1:2,:,ipp2-2:ipp2-ipph:(-1)), &
                            shape = [(ido-1)/2, ipph-1, l1], &
                            order = [1, 3, 2])

                        cc(1,3:ido:2, 3:ipph*2-1:2,:) = &
                            reshape(source = &
                            ch(1,3:ido:2,:,2:ipph)+ch(1,3:ido:2,:,ipp2-2:ipp2-ipph:(-1)), &
                            shape = [(ido-1)/2, ipph-1, l1], &
                            order = [1, 3, 2])

                        cc(1,idp2-3:idp2-ido:(-2), 2:(ipph-1)*2:2,:) = &
                            reshape(source = &
                            ch(1,3:ido:2,:,ipp2-2:ipp2-ipph:(-1))-ch(1,3:ido:2,:, 2:ipph), shape &
                            = [(ido-1)/2, ipph-1, l1], order = [1, 3, 2])
                    else
                        cc(1,2:ido-1:2, 3:ipph*2-1:2,:) = &
                            reshape(source = &
                            ch(1,2:ido-1:2,:, 2:ipph)+ch(1,2:ido-1:2,:,ipp2-2:ipp2-ipph:(-1)), &
                            shape = [(ido-1)/2, ipph-1, l1], &
                            order = [1, 3, 2])

                        cc(1,idp2-4:idp2-1-ido:(-2), 2:(ipph-1)*2:2,:) = &
                            reshape(source = &
                            ch(1,2:ido-1:2,:, 2:ipph)-ch(1,2:ido-1:2,:,ipp2-2:ipp2-ipph:(-1)), &
                            shape = [(ido-1)/2, ipph-1, l1], &
                            order = [1, 3, 2])

                        cc(1,3:ido:2, 3:ipph*2-1:2,:) = &
                            reshape(source = &
                            ch(1,3:ido:2,:, 2:ipph) +ch(1,3:ido:2,:, ipp2-2:ipp2-ipph:(-1)), &
                            shape = [(ido-1)/2, ipph-1, l1], &
                            order = [1, 3, 2])

                        cc(1,idp2-3:idp2-ido:(-2), 2:(ipph-1)*2:2,:) = &
                            reshape(source = &
                            ch(1,3:ido:2,:,ipp2-2:ipp2-ipph:(-1))-ch(1,3:ido:2,:,2:ipph), &
                            shape = [(ido-1)/2, ipph-1, l1], order = [1, 3, 2])
                    end if
            end select

        end subroutine r1fgkf

    end subroutine rfft1f



    subroutine rfft1i(n, wsave, lensav, ier)
        !
        !! RFFT1I: initialization for RFFT1B and RFFT1F.
        !
        !  Purpose:
        !
        !  RFFT1I initializes array wsave for use in its companion routines
        !  RFFT1B and RFFT1F.  The prime factorization of N together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  Separate wsave arrays are required for different
        !  values of N.
        !
        !  Parameters:
        !
        !  integer N, the length of the sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  real wsave(LENSAV), containing the prime factors of
        !  N and also containing certain trigonometric values which will be used in
        !  routines RFFT1B or RFFT1F.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough.
        !


        integer (ip) lensav

        integer (ip) ier
        integer (ip) n
        real (wp) wsave(lensav)

        ier = 0

        if (lensav < n + int(log(real(n, kind=wp) )/log(2.0_wp)) +4) then
            ier = 2
            call fft_error_handler('rfft1i ', 3)
        end if

        if (n /= 1) then
            call rffti1(n,wsave(1),wsave(n+1))
        end if

    end subroutine rfft1i


    subroutine rfft2b ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier)
        !
        ! rfft2b: 64-bit float precision backward fast Fourier transform, 2D.
        !
        ! Purpose:
        !
        !  Computes the two-dimensional discrete Fourier transform of the
        !  complex Fourier coefficients a real periodic array.  This transform is
        !  known as the backward transform or Fourier synthesis, transforming from
        !  spectral to physical space.  Routine RFFT2B is normalized: a call to
        !  RFFT2B followed by a call to RFFT2F (or vice-versa) reproduces the
        !  original array within roundoff error.
        !
        !  Parameters:
        !
        !  integer LDIM, the first dimension of the 2D real
        !  array R, which must be at least 2*(L/2+1).
        !
        !  integer L, the number of elements to be transformed
        !  in the first dimension of the two-dimensional real array R.  The value of
        !  L must be less than or equal to that of LDIM.  The transform is most
        !  efficient when L is a product of small primes.
        !
        !  integer M, the number of elements to be transformed
        !  in the second dimension of the two-dimensional real array R.  The transform
        !  is most efficient when M is a product of small primes.
        !
        !  Input/real R(LDIM,M), the real array of two
        !  dimensions.  On input, R contains the L/2+1-by-M complex subarray of
        !  spectral coefficients, on output, the physical coefficients.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFT2I before the first call to routine RFFT2F
        !  or RFFT2B with lengths L and M.  wsave's contents may be re-used for
        !  subsequent calls to RFFT2F and RFFT2B with the same transform lengths
        !  L and M.
        !
        !  integer LENSAV, the number of elements in the wsave
        !  array.  LENSAV must be at least L + M + INT(LOG(REAL(L)))
        !  + INT(LOG(REAL(M))) + 8.
        !
        !  Workspace, real (wp) WORK(LENWRK).  WORK provides workspace, and
        !  its contents need not be saved between calls to routines RFFT2B and RFFT2F.
        !
        !  integer  LENWRK, the number of elements in the WORK
        !  array.  LENWRK must be at least LDIM*M.
        !
        !  integer IER, the error flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  6, input parameter LDIM < 2*(L/2+1);
        !  20, input error returned by lower level routine.
        !
        integer (ip) ldim
        integer (ip) lensav
        integer (ip) lenwrk
        integer (ip) m

        
        integer (ip) ier
        integer (ip) local_error_flag
        
        integer (ip) l
        integer (ip) ldh
        integer (ip) ldw
        integer (ip) ldx
        integer (ip) lwsav
        integer (ip) mmsav
        integer (ip) modl
        integer (ip) modm
        integer (ip) mwsav
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) r(ldim,m)

        ier = 0
        !
        !==> verify lensav
        !
        lwsav = l+int(log(real(l, kind=wp))/log(2.0_wp))+4
        mwsav = get_1d_saved_workspace_size(m)
        mmsav = m+int(log(real(m, kind=wp))/log(2.0_wp))+4
        modl = mod(l,2)
        modm = mod(m,2)

        if (lensav < lwsav+mwsav+mmsav) then
            ier = 2
            call fft_error_handler('rfft2f', 6)
            return
        end if
        !
        ! verify lenwrk
        !
        if (lenwrk < (l+1)*m) then
            ier = 3
            call fft_error_handler('rfft2f', 8)
            return
        end if
        !
        ! verify ldim is as big as l
        !
        if (ldim < l) then
            ier = 5
            call fft_error_handler('rfft2f', -6)
            return
        end if
        !
        !==> Transform second dimension of array
        !
        associate( mj => 2*((m+1)/2)-1 )

            r(1,2:mj) = 2.0_wp * r(1,2:mj)

        end associate

        r(1,3:m:2) = -r(1,3:m:2)

        call rfftmb(1,1,m,ldim,r,m*ldim, wsave(lwsav+mwsav+1),mmsav,work,lenwrk,local_error_flag)

        ldh = int((l+1)/2)

        if ( 1 < ldh ) then
            ldw = ldh+ldh
            !
            !  r and work are switched because the the first dimension
            !  of the input to complex cfftmf must be even.
            !
            call r2w(ldim,ldw,l,m,r,work)

            call cfftmb(ldh-1,1,m,ldh,work(2),ldh*m, &
                wsave(lwsav+1),mwsav,r,l*m, local_error_flag)

            if (local_error_flag/=0) then
                ier=20
                call fft_error_handler('rfft2b',-5)
                return
            end if

            call w2r(ldim,ldw,l,m,r,work)
        end if

        if (modl == 0) then

            associate( mj => 2*((m+1)/2)-1 )

                r(l,2:mj) = 2.0_wp * r(l,2:mj)

            end associate

            r(l,3:m:2) = -r(l,3:m:2)

            call rfftmb(1,1,m,ldim,r(l,1),m*ldim, &
                wsave(lwsav+mwsav+1),mmsav,work,lenwrk,local_error_flag)
        end if
        !
        !==> Transform first dimension of array
        !
        ldx = 2*int((l+1)/2)-1

        r(2:ldx,1:m) = 2.0_wp * r(2:ldx,1:m)
        r(3:ldx:2,1:m) = -r(3:ldx:2,1:m)

        associate( &
            arg_1 => m*ldim, &
            arg_2 => l+int(log(real(l, kind=wp) )/log(2.0_wp))+4 &
            )

            call rfftmb(m,ldim,l,1,r,arg_1,wsave(1), arg_2,work,lenwrk,local_error_flag)

        end associate

        if (local_error_flag /= 0) then
            ier=20
            call fft_error_handler('rfft2f',-5)
        end if

    end subroutine rfft2b


    subroutine rfft2f(ldim, l, m, r, wsave, lensav, work, lenwrk, ier)
        !
        ! rfft2f: 64-bit float precision forward fast Fourier transform, 2D.
        !
        ! Purpose:
        !
        !  RFFT2F computes the two-dimensional discrete Fourier transform of a
        !  real periodic array.  This transform is known as the forward transform
        !  or Fourier analysis, transforming from physical to spectral space.
        !  Routine rfft2f is normalized: a call to rfft2f followed by a call to
        !  rfft2b (or vice-versa) reproduces the original array within roundoff
        !  error.
        !
        !  Parameters:
        !
        !  integer LDIM, the first dimension of the 2D real
        !  array R, which must be at least 2*(L/2+1).
        !
        !  integer L, the number of elements to be transformed
        !  in the first dimension of the two-dimensional real array R.  The value
        !  of L must be less than or equal to that of LDIM.  The transform is most
        !  efficient when L is a product of small primes.
        !
        !  integer M, the number of elements to be transformed
        !  in the second dimension of the two-dimensional real array R.  The
        !  transform is most efficient when M is a product of small primes.
        !
        !  Input/real R(LDIM,M), the real array of two
        !  dimensions.  On input, containing the L-by-M physical data to be
        !  transformed.  On output, the spectral coefficients.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFT2I before the first call to routine RFFT2F
        !  or RFFT2B with lengths L and M.  wsave's contents may be re-used for
        !  subsequent calls to RFFT2F and RFFT2B with the same transform lengths.
        !
        !  integer LENSAV, the number of elements in the wsave
        !  array.  LENSAV must be at least L + M + INT(LOG(REAL(L)))
        !  + INT(LOG(REAL(M))) + 8.
        !
        !  Workspace, real (wp) WORK(LENWRK), provides workspace, and its
        !  contents need not be saved between calls to routines RFFT2F and RFFT2B.
        !
        !  integer LENWRK, the number of elements in the WORK
        !  array.  LENWRK must be at least LDIM*M.
        !
        !  integer IER, the error flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  6, input parameter LDIM < 2*(L+1);
        !  20, input error returned by lower level routine.
        !


        integer (ip) ldim
        integer (ip) lensav
        integer (ip) lenwrk
        integer (ip) m

        
        integer (ip) ier
        integer (ip) local_error_flag
        
        integer (ip) l
        integer (ip) ldh
        integer (ip) ldw
        integer (ip) ldx
        integer (ip) lwsav
        integer (ip) mmsav
        integer (ip) modl
        integer (ip) modm
        integer (ip) mwsav
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) r(ldim,m)

        ier = 0
        !
        !==> verify lensav
        !
        lwsav = l+int(log(real (l, kind=wp))/log(2.0_wp))+4
        mwsav = get_1d_saved_workspace_size(m)
        mmsav = m+int(log(real(m, kind=wp))/log(2.0_wp))+4

        if (lensav < lwsav+mwsav+mmsav) then
            ier = 2
            call fft_error_handler('rfft2f', 6)
            return
        end if
        !
        !==>  verify lenwrk
        !
        if (lenwrk < (l+1)*m) then
            ier = 3
            call fft_error_handler('rfft2f', 8)
            return
        end if
        !
        !==>  verify ldim is as big as l
        !
        if (ldim < l) then
            ier = 5
            call fft_error_handler('rfft2f', -6)
            return
        end if
        !
        !==>  Transform first dimension of array
        !
        associate( &
            arg_1 => m*ldim, &
            arg_2 => l+int(log(real(l, kind=wp))/log(2.0_wp))+4 &
            )

            call rfftmf(m,ldim,l,1,r,arg_1,wsave(1), arg_2,work,size(work),local_error_flag)

        end associate

        if (local_error_flag /= 0) then
            ier=20
            call fft_error_handler('rfft2f',-5)
            return
        end if

        ldx = 2*int((l+1)/2)-1

        r(2:ldx,1:m) = 0.5_wp * r(2:ldx,1:m)

        r(3:ldx:2,1:m) = -r(3:ldx:2,1:m)

        !==>  Reshuffle to add in nyquist imaginary components
        !
        modl = mod(l,2)
        modm = mod(m,2)
        !
        !==>  Transform second dimension of array
        !
        call rfftmf(1,1,m,ldim,r,m*ldim, &
            wsave(lwsav+mwsav+1),mmsav,work,size(work),local_error_flag)

        associate( mj => 2*((m+1)/2)-1 )

            r(1,2:mj) = 0.5_wp * r(1,2:mj)

        end associate

        r(1,3:m:2) = -r(1,3:m:2)

        ldh = int((l+1)/2)

        if ( 1 < ldh ) then
            ldw = 2*ldh
            !
            !==> r and work are switched because the the first dimension
            !    of the input to complex cfftmf must be even.
            !
            call r2w(ldim,ldw,l,m,r,work)
            call cfftmf(ldh-1,1,m,ldh,work(2),ldh*m, &
                wsave(lwsav+1),mwsav,r,l*m, local_error_flag)

            if (local_error_flag /= 0) then
                ier=20
                call fft_error_handler('rfft2f',-5)
                return
            end if

            call w2r(ldim,ldw,l,m,r,work)
        end if

        if (modl == 0) then

            call rfftmf(1,1,m,ldim,r(l,1),m*ldim, &
                wsave(lwsav+mwsav+1),mmsav,work,size(work),local_error_flag)

            associate( mj => 2*((m+1)/2)-1 )

                r(l,2:mj) = 0.5_wp * r(l,2:mj)

            end associate

            r(l,3:m:2) = -r(l,3:m:2)

        end if

        if (local_error_flag /= 0) then
            ier=20
            call fft_error_handler('rfft2f',-5)
        end if

    end subroutine rfft2f



    subroutine rfft2i(l, m, wsave, lensav, ier)
        ! RFFT2I: initialization for RFFT2B and RFFT2F.
        !
        !  Purpose:
        !  RFFT2I initializes real array wsave for use in its companion routines
        !  RFFT2F and RFFT2B for computing the two-dimensional fast Fourier
        !  transform of real data.  Prime factorizations of L and M, together with
        !  tabulations of the trigonometric functions, are computed and stored in
        !  array wsave.  RFFT2I must be called prior to the first call to RFFT2F
        !  or RFFT2B.  Separate wsave arrays are required for different values of
        !  L or M.
        !
        !
        !  integer L, the number of elements to be transformed
        !  in the first dimension.  The transform is most efficient when L is a
        !  product of small primes.
        !
        !  integer M, the number of elements to be transformed
        !  in the second dimension.  The transform is most efficient when M is a
        !  product of small primes.
        !
        !  integer LENSAV, the number of elements in the wsave
        !  array.  LENSAV must be at least L + M + INT(LOG(REAL(L)))
        !  + INT(LOG(REAL(M))) + 8.
        !
        !  real wsave(LENSAV), containing the prime factors
        !  of L and M, and also containing certain trigonometric values which
        !  will be used in routines RFFT2B or RFFT2F.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  20, input error returned by lower level routine.
        !


        integer (ip) lensav

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) l
        integer (ip) lwsav
        integer (ip) m
        integer (ip) mmsav
        integer (ip) mwsav
        real (wp) wsave(lensav)
        !
        ! initialize ier
        !
        ier = 0
        !
        ! verify lensav
        !
        lwsav = l+int(log(real(l, kind=wp) )/log(2.0_wp))+4
        mwsav = 2*m+int(log(real(m, kind=wp) )/log(2.0_wp))+4
        mmsav = m+int(log(real(m, kind=wp) )/log(2.0_wp))+4

        if (lensav < lwsav+mwsav+mmsav) then
            ier = 2
            call fft_error_handler('rfft2i', 4)
            return
        end if

        call rfftmi(l, wsave(1), lwsav, local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('rfft2i',-5)
            return
        end if

        call cfftmi(m, wsave(lwsav+1),mwsav,local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('rfft2i',-5)
            return
        end if

        call rfftmi(m,wsave(lwsav+mwsav+1),mmsav, local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('rfft2i',-5)
            return
        end if

    end subroutine rfft2i


    subroutine rffti1(n, wa, fac)
        !
        !  Parameters:
        !
        !  input
        !
        ! n, the number for which factorization
        !  and other information is needed.
        !
        !  output
        ! wa(n), trigonometric information.
        !
        !  output
        !
        !  fac(15), factorization information.
        !  fac(1) is n, fac(2) is nf, the number of factors, and fac(3:nf+2) are the
        !  factors.
        !
        !--------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: fac(15)
        real (wp),    intent (out) :: wa(n)
        !--------------------------------------------------------------
        ! Dictionary: local variables
        !--------------------------------------------------------------
        integer (ip)            :: i, ib, ido, ii, iip, ipm, is
        integer (ip)            :: j, k1, l1, l2, ld
        integer (ip)            :: nf, nfm1, nl, nq, nr, ntry
        integer (ip), parameter :: ntryh(*)=[ 4, 2, 3, 5]
        real (wp),    parameter :: TWO_PI = 2.0_wp * acos(-1.0_wp)
        real (wp)               :: arg,  argh, argld, fi
        !--------------------------------------------------------------

        ntry = 0
        nl = n
        nf = 0
        j = 0

        factorize_loop: do
            ! Increment j
            j = j+1

            ! Choose ntry
            if (j <= 4) then
                ntry = ntryh(j)
            else
                ntry = ntry+2
            end if

            inner_loop: do
                nq = nl/ntry
                nr = nl-ntry*nq
                if (nr < 0) then
                    cycle factorize_loop
                else if (nr == 0) then
                    nf = nf+1
                    fac(nf+2) = ntry
                    nl = nq

                    if (ntry == 2 .and. nf /= 1) then
                        do i=2,nf
                            ib = nf-i+2
                            fac(ib+2) = fac(ib+1)
                        end do
                        fac(3) = 2
                    end if

                    if (nl /= 1) then
                        cycle inner_loop
                    end if
                else
                    cycle factorize_loop
                end if
                exit inner_loop
            end do inner_loop
            exit factorize_loop
        end do factorize_loop

        fac(1) = n
        fac(2) = nf
        argh = TWO_PI/n
        is = 0
        nfm1 = nf-1
        l1 = 1

        if (nfm1 /= 0) then
            do k1=1,nfm1
                iip = int(fac(k1+2), kind=ip)
                ld = 0
                l2 = l1*iip
                ido = n/l2
                ipm = iip-1
                do j=1,ipm
                    ld = ld+l1
                    i = is
                    argld = real(ld, kind=wp) * argh
                    fi = 0.0_wp
                    do ii=3,ido,2
                        i = i+2
                        fi = fi + 1.0_wp
                        arg = fi*argld
                        wa(i-1) = cos(arg)
                        wa(i) = sin(arg)
                    end do
                    is = is+ido
                end do
                l1 = l2
            end do
        end if

    end subroutine rffti1



    subroutine rfftmb(lot, jump, n, inc, r, lenr, wsave, &
        lensav, work, lenwrk, ier)
        !
        ! rfftmb: 64-bit float precision backward fft, 1d, multiple vectors.
        !
        ! Purpose:
        !
        !  Computes the one-dimensional Fourier transform of multiple
        !  periodic sequences within a real array.  This transform is referred
        !  to as the backward transform or Fourier synthesis, transforming the
        !  sequences from spectral to physical space.
        !
        !  This transform is normalized since a call to RFFTMB followed
        !  by a call to RFFTMF (or vice-versa) reproduces the original
        !  array  within roundoff error.
        !
        !  Parameters:
        !
        !  integer LOT, the number of sequences to be transformed
        !  within array R.
        !
        !  integer JUMP, the increment between the locations, in
        !  array R, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations, in
        !  array R, of two consecutive elements within the same sequence.
        !
        !  Input/real R(LENR), real array containing LOT
        !  sequences, each having length N.  R can have any number of dimensions,
        !  but the total number of locations must be at least LENR.  On input, the
        !  spectral data to be transformed, on output the physical data.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least (LOT-1)*JUMP + INC*(N-1) + 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFTMI before the first call to routine RFFTMF
        !  or RFFTMB for a given transform length N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must  be at least N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least LOT*N.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  4, input parameters INC, JUMP, N, LOT are not consistent.
        !
        integer (ip) lenr
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) inc
        integer (ip) jump
        integer (ip) lot
        integer (ip) n
        real (wp) r(lenr)
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)

        !
        !==> Validity of calling arguments
        !
        if (lenr < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('rfftmb ', 6)
            return
        else if (lensav < n+int(log(real(n, kind=wp))/log(2.0_wp))+4) then
            ier = 2
            call fft_error_handler('rfftmb ', 8)
            return
        else if (lenwrk < lot*n) then
            ier = 3
            call fft_error_handler('rfftmb ', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('rfftmb ', -1)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            call mrftb1(lot,jump,n,inc,r,work,wsave,wsave(n+1))
        end if

    end subroutine rfftmb


    subroutine rfftmf(lot, jump, n, inc, r, lenr, wsave, lensav, &
        work, lenwrk, ier)
        !
        ! RFFTMF: 64-bit float precision forward FFT, 1D, multiple vectors.
        !
        !  Purpose:
        !
        !  RFFTMF computes the one-dimensional Fourier transform of multiple
        !  periodic sequences within a real array.  This transform is referred
        !  to as the forward transform or Fourier analysis, transforming the
        !  sequences from physical to spectral space.
        !
        !  This transform is normalized since a call to RFFTMF followed
        !  by a call to RFFTMB (or vice-versa) reproduces the original array
        !  within roundoff error.
        !
        !  Parameters:
        !
        !  integer LOT, the number of sequences to be transformed
        !  within array R.
        !
        !  integer JUMP, the increment between the locations, in
        !  array R, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations,
        !  in array R, of two consecutive elements within the same sequence.
        !
        !  Input/real R(LENR), real array containing LOT
        !  sequences, each having length N.  R can have any number of dimensions, but
        !  the total number of locations must be at least LENR.  On input, the
        !  physical data to be transformed, on output the spectral data.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least (LOT-1)*JUMP + INC*(N-1) + 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFTMI before the first call to routine RFFTMF
        !  or RFFTMB for a given transform length N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least LOT*N.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  4, input parameters INC, JUMP, N, LOT are not consistent.
        !


        integer (ip) lenr
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) inc
        integer (ip) jump
        integer (ip) lot
        integer (ip) n
        real (wp) r(lenr)
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lenr < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('rfftmf ', 6)
            return
        else if (lensav < n + int(log(real(n, kind=wp) )/log(2.0_wp)) +4) then
            ier = 2
            call fft_error_handler('rfftmf ', 8)
            return
        else if (lenwrk < lot*n) then
            ier = 3
            call fft_error_handler('rfftmf ', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('rfftmf ', -1)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            call mrftf1(lot,jump,n,inc,r,work,wsave,wsave(n+1))
        end if

    end subroutine rfftmf

    subroutine rfftmi(n, wsave, lensav, ier)
        !
        ! rfftmi: initialization for rfftmb and rfftmf.
        !
        !  Purpose:
        !
        !  rfftmi initializes array wsave for use in its companion routines
        !  rfftmb and rfftmf.  the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  separate wsave arrays are required for different
        !  values of n.
        !
        !  INPUT
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n + int(log(real(n))) + 4.
        !
        !  OUTPUT
        !  real wsave(lensav), work array containing the prime
        !  factors of n and also containing certain trigonometric
        !  values which will be used in routines rfftmb or rfftmf.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough.
        !


        integer (ip) lensav

        integer (ip) ier
        integer (ip) n
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lensav < n + int(log(real(n, kind=wp) )/log(2.0_wp)) +4) then
            ier = 2
            call fft_error_handler('rfftmi ', 3)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            call mrfti1(n,wsave(1),wsave(n+1))
        end if

    end subroutine rfftmi


    subroutine sinq1b(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! sinq1b: 64-bit float precision backward sine quarter wave transform, 1D.
        !
        !  Purpose:
        !
        !  Computes the one-dimensional Fourier transform of a sequence
        !  which is a sine series with odd wave numbers.  This transform is
        !  referred to as the backward transform or Fourier synthesis,
        !  transforming the sequence from spectral to physical space.
        !
        !  This transform is normalized since a call to sinq1b followed
        !  by a call to sinq1f (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  Parameters:
        !
        !  integer N, the length of the sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations, in
        !  array R, of two consecutive elements within the sequence.
        !
        !  Input/real R(LENR), on input, the sequence to be
        !  transformed.  On output, the transformed sequence.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least INC*(N-1)+ 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to SINQ1I before the first call to routine SINQ1F
        !  or SINQ1B for a given transform length N.  wsave's contents may be
        !  re-used for subsequent calls to SINQ1F and SINQ1B with the same N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least N.
        !
        !  integer IER, the error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  20, input error returned by lower level routine.
        !

        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) n
        integer (ip) ns2
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)
        real (wp) xhold

        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('sinq1b', 6)
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('sinq1b', 8)
        else if (lenwrk < n) then
            ier = 3
            call fft_error_handler('sinq1b', 10)
        else
            ier = 0
        end if

        if (1 >= n) then
            !
            !   x(1,1) = 4.*x(1,1) line disabled by dick valent 08/26/2010
            !
            return
        else
            ns2 = n/2
            x(1,2:n:2) = -x(1,2:n:2)

            call cosq1b(n,inc,x,lenx,wsave,lensav,work,lenwrk,local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sinq1b',-5)
                return
            end if

            do k=1,ns2
                kc = n-k
                xhold = x(1,k)
                x(1,k) = x(1,kc+1)
                x(1,kc+1) = xhold
            end do
        end if

    end subroutine sinq1b


    subroutine sinq1f(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! sinq1f: 64-bit float precision forward sine quarter wave transform, 1D.
        !
        !  Purpose:
        !
        !  Computes the one-dimensional Fourier transform of a sequence
        !  which is a sine series of odd wave numbers.  This transform is
        !  referred to as the forward transform or Fourier analysis, transforming
        !  the sequence from physical to spectral space.
        !
        !  This transform is normalized since a call to sinq1f followed
        !  by a call to sinq1b (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  Parameters:
        !
        !  integer N, the length of the sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations,
        !  in array R, of two consecutive elements within the sequence.
        !
        !  Input/real R(LENR), on input, the sequence to be
        !  transformed.  On output, the transformed sequence.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least INC*(N-1)+ 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to SINQ1I before the first call to routine SINQ1F
        !  or SINQ1B for a given transform length N.  wsave's contents may be re-used
        !  for subsequent calls to SINQ1F and SINQ1B with the same N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least N.
        !
        !  integer IER, the error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  20, input error returned by lower level routine.
        !
        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) n
        integer (ip) ns2
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)
        real (wp) xhold

        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('sinq1f', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('sinq1f', 8)
            return
        else if (lenwrk < n) then
            ier = 3
            call fft_error_handler('sinq1f', 10)
            return
        else
            ier = 0
        end if

        if (n /= 1) then
            ns2 = n/2
            do k=1,ns2
                kc = n-k
                xhold = x(1,k)
                x(1,k) = x(1,kc+1)
                x(1,kc+1) = xhold
            end do

            call cosq1f(n,inc,x,lenx,wsave,lensav,work,lenwrk,local_error_flag)

            ! check error flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sinq1f',-5)
                return
            end if

            x(1,2:n:2) = -x(1,2:n:2)
        end if

    end subroutine sinq1f


    subroutine sinq1i(n, wsave, lensav, ier)
        !
        !  sinq1i: initialization for sinq1b and sinq1f.
        !
        !  Purpose:
        !
        !  Initializes array wsave for use in its companion routines
        !  sinq1f and sinq1b. The prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave. Separate wsave arrays are required for different
        !  values of n.
        !
        !  Parameters:
        !
        !  integer N, the length of the sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
        !
        !  real wsave(LENSAV), containing the prime factors
        !  of N and also containing certain trigonometric values which will be used
        ! in routines SINQ1B or SINQ1F.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  20, input error returned by lower level routine.
        !

        integer (ip) lensav
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) n
        real (wp) wsave(lensav)

        if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('sinq1i', 3)
            return
        else
            ier = 0
        end if

        call cosq1i(n, wsave, lensav, local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('sinq1i',-5)
        end if

    end subroutine sinq1i

    subroutine sinqmb(lot, jump, n, inc, x, lenx, wsave, lensav, &
        work, lenwrk, ier)
        !
        ! SINQMB: 64-bit float precision backward sine quarter wave, multiple vectors.
        !
        !  Purpose:
        !
        !  SINQMB computes the one-dimensional Fourier transform of multiple
        !  sequences within a real array, where each of the sequences is a
        !  sine series with odd wave numbers.  This transform is referred to as
        !  the backward transform or Fourier synthesis, transforming the
        !  sequences from spectral to physical space.
        !
        !  This transform is normalized since a call to SINQMB followed
        !  by a call to SINQMF (or vice-versa) reproduces the original
        !  array  within roundoff error.
        !
        !  Parameters:
        !
        !  integer LOT, the number of sequences to be transformed
        !  within array R.
        !
        !  integer JUMP, the increment between the locations, in
        !  array R, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations, in
        !  array R, of two consecutive elements within the same sequence.
        !
        !  Input/real R(LENR), containing LOT sequences, each
        !  having length N.  R can have any number of dimensions, but the total
        !  number of locations must be at least LENR.  On input, R contains the data
        !  to be transformed, and on output the transformed data.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to SINQMI before the first call to routine SINQMF
        !  or SINQMB for a given transform length N.  wsave's contents may be re-used
        !  for subsequent calls to SINQMF and SINQMB with the same N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least LOT*N.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  4, input parameters INC,JUMP,N,LOT are not consistent;
        !  20, input error returned by lower level routine.
        !

        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) lj
        integer (ip) lot
        integer (ip) m
        integer (ip) n
        integer (ip) ns2
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)
        real (wp) xhold

        if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('sinqmb', 6)
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('sinqmb', 8)
        else if (lenwrk < lot*n) then
            ier = 3
            call fft_error_handler('sinqmb', 10)
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('sinqmb', -1)
        else
            ier = 0
        end if

        lj = (lot-1)*jump+1

        if (1 >= n ) then
            x(1:lj:jump,1) =  4.0_wp * x(1:lj:jump,1)
        else
            ns2 = n/2
            x(2:n:2,1:lj:jump) = -x(1:lj:jump,2:n:2)

            call cosqmb(lot,jump,n,inc,x,lenx,wsave,lensav,work,lenwrk,local_error_flag)

            ! Check error flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sinqmb',-5)
                return
            end if

            do k=1,ns2
                kc = n-k
                do m=1,lj,jump
                    xhold = x(m,k)
                    x(m,k) = x(m,kc+1)
                    x(m,kc+1) = xhold
                end do
            end do
        end if

    end subroutine sinqmb

    subroutine sinqmf(lot, jump, n, inc, x, lenx, wsave, lensav, &
        work, lenwrk, ier)
        !
        ! SINQMF: 64-bit float precision forward sine quarter wave, multiple vectors.
        !
        !  Purpose:
        !
        !  SINQMF computes the one-dimensional Fourier transform of multiple
        !  sequences within a real array, where each sequence is a sine series
        !  with odd wave numbers.  This transform is referred to as the forward
        !  transform or Fourier synthesis, transforming the sequences from
        !  spectral to physical space.
        !
        !  This transform is normalized since a call to SINQMF followed
        !  by a call to SINQMB (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  Parameters:
        !
        !  integer LOT, the number of sequences to be transformed
        !  within array R.
        !
        !  integer JUMP, the increment between the locations,
        !  in array R, of the first elements of two consecutive sequences to
        !  be transformed.
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations,
        !  in array R, of two consecutive elements within the same sequence.
        !
        !  Input/real R(LENR), containing LOT sequences, each
        !  having length N.  R can have any number of dimensions, but the total
        !  number of locations must be at least LENR.  On input, R contains the data
        !  to be transformed, and on output the transformed data.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to SINQMI before the first call to routine SINQMF
        !  or SINQMB for a given transform length N.  wsave's contents may be re-used
        !  for subsequent calls to SINQMF and SINQMB with the same N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least LOT*N.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  4, input parameters INC,JUMP,N,LOT are not consistent;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) jump
        integer (ip) k
        integer (ip) kc
        integer (ip) lenx
        integer (ip) lj
        integer (ip) lot
        integer (ip) m
        integer (ip) n
        integer (ip) ns2
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)
        real (wp) xhold


        if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('sinqmf', 6)
            return
        else if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('sinqmf', 8)
            return
        else if (lenwrk < lot*n) then
            ier = 3
            call fft_error_handler('sinqmf', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('sinqmf', -1)
            return
        else
            ier = 0
        end if

        if (n /= 1) then

            ns2 = n/2
            lj = (lot-1)*jump+1

            do k=1,ns2
                kc = n-k
                do m=1,lj,jump
                    xhold = x(m,k)
                    x(m,k) = x(m,kc+1)
                    x(m,kc+1) = xhold
                end do
            end do

            call cosqmf(lot,jump,n,inc,x,lenx,wsave,lensav,work,lenwrk,local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sinqmf',-5)
                return
            end if

            x(1:lj:jump,2:n:2) = -x(1:lj:jump,2:n:2)
        end if

    end subroutine sinqmf



    subroutine sinqmi(n, wsave, lensav, ier)


        !
        !! SINQMI: initialization for SINQMB and SINQMF.
        !
        !  Purpose:
        !
        !  SINQMI initializes array wsave for use in its companion routines
        !  SINQMF and SINQMB.  The prime factorization of N together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  Separate wsave arrays are required for different
        !  values of N.
        !
        !  Parameters:
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
        !
        !  real wsave(LENSAV), containing the prime factors
        !  of N and also containing certain trigonometric values which will be used
        !  in routines SINQMB or SINQMF.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  20, input error returned by lower level routine.
        !


        integer (ip) lensav
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) n
        real (wp) wsave(lensav)

        ier = 0

        if (lensav < get_1d_saved_workspace_size(n)) then
            ier = 2
            call fft_error_handler('sinqmi', 3)
            return
        end if

        call cosqmi(n, wsave, lensav, local_error_flag)

        ! Check error flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('sinqmi',-5)
        end if

    end subroutine sinqmi

    subroutine sint1b(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)


        !
        !! SINT1B: 64-bit float precision backward sine transform, 1D.
        !
        !  Purpose:
        !
        !  SINT1B computes the one-dimensional Fourier transform of an odd
        !  sequence within a real array.  This transform is referred to as
        !  the backward transform or Fourier synthesis, transforming the
        !  sequence from spectral to physical space.
        !
        !  This transform is normalized since a call to SINT1B followed
        !  by a call to SINT1F (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  Parameters:
        !
        !  integer N, the length of the sequence to be
        !  transformed.  The transform is most efficient when N+1 is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations, in
        !  array R, of two consecutive elements within the sequence.
        !
        !  Input/real R(LENR), on input, contains the sequence
        !  to be transformed, and on output, the transformed sequence.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least INC*(N-1)+ 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to SINT1I before the first call to routine SINT1F
        !  or SINT1B for a given transform length N.  wsave's contents may be re-used
        !  for subsequent calls to SINT1F and SINT1B with the same N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least 2*N+2.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  20, input error returned by lower level routine.
        !

        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) lenx
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Check validity of input arguments
        !
        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('sint1b', 6)
            return
        else if ( lensav < n / 2 + n + int(log(real(n, kind=wp) ) &
            / log(2.0_wp ) ) + 4 ) then
            ier = 2
            call fft_error_handler('sint1b', 8)
            return
        else if (lenwrk < (2*n+2)) then
            ier = 3
            call fft_error_handler('sint1b', 10)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        call sintb1(n,inc,x,wsave,work,work(n+2),local_error_flag)

        ! Check error flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('sint1b',-5)
        end if

    end subroutine sint1b


    subroutine sint1f(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        !! SINT1F: 64-bit float precision forward sine transform, 1D.
        !
        !  Purpose:
        !
        !  SINT1F computes the one-dimensional Fourier transform of an odd
        !  sequence within a real array.  This transform is referred to as the
        !  forward transform or Fourier analysis, transforming the sequence
        !  from physical to spectral space.
        !
        !  This transform is normalized since a call to SINT1F followed
        !  by a call to SINT1B (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  Parameters:
        !
        !  integer N, the length of the sequence to be
        !  transformed.  The transform is most efficient when N+1 is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations,
        !  in array R, of two consecutive elements within the sequence.
        !
        !  Input/real R(LENR), on input, contains the sequence
        !  to be transformed, and on output, the transformed sequence.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least INC*(N-1)+ 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to SINT1I before the first call to routine SINT1F
        !  or SINT1B for a given transform length N.  wsave's contents may be re-used
        !  for subsequent calls to SINT1F and SINT1B with the same N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least 2*N+2.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) lenx
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Check validity of input arguments
        !
        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('sint1f', 6)
            return
        else if (lensav < n/2 + n + int(log(real(n, kind=wp) ) &
            /log(2.0_wp)) +4) then
            ier = 2
            call fft_error_handler('sint1f', 8)
            return
        else if (lenwrk < (2*n+2)) then
            ier = 3
            call fft_error_handler('sint1f', 10)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        call sintf1(n,inc,x,wsave,work,work(n+2),local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('sint1f',-5)
        end if

    end subroutine sint1f



    subroutine sint1i(n, wsave, lensav, ier)
        !
        !! SINT1I: initialization for SINT1B and SINT1F.
        !
        !  Purpose:
        !
        !  SINT1I initializes array wsave for use in its companion routines
        !  SINT1F and SINT1B.  The prime factorization of N together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  Separate wsave arrays are required for different
        !  values of N.
        !
        !  Parameters:
        !
        !  integer N, the length of the sequence to be
        !  transformed.  The transform is most efficient when N+1 is a product
        !  of small primes.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
        !
        !  real wsave(LENSAV), containing the prime factors
        !  of N and also containing certain trigonometric values which will be used
        !  in routines SINT1B or SINT1F.
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  20, input error returned by lower level routine.
        !


        integer (ip) lensav

        real (wp) dt
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) lnsv
        integer (ip) n
        integer (ip) np1
        integer (ip) ns2
        real (wp), parameter :: pi = acos(-1.0_wp)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lensav < n/2 + n + int(log(real(n, kind=wp) ) &
            /log(2.0_wp)) +4) then
            ier = 2
            call fft_error_handler('sint1i', 3)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n > 1) then

            ns2 = n/2
            np1 = n+1
            dt = pi/np1

            do k=1,ns2
                wsave(k) = 2.0_wp *sin(k*dt)
            end do

            lnsv = np1 + int(log(real(np1, kind=wp))/log(2.0_wp)) +4

            call rfft1i(np1, wsave(ns2+1), lnsv, local_error_flag)

            ! Check error flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sint1i',-5)
            end if
        end if

    end subroutine sint1i



    subroutine sintb1(n, inc, x, wsave, xh, work, ier)

        integer (ip) inc
        real (wp) dsum
        real (wp) fnp1s4
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) lnxh
        integer (ip) modn
        integer (ip) n
        integer (ip) np1
        integer (ip) ns2
        real (wp), parameter :: HALF_SQRT3 = sqrt(3.0_wp)/2
        real (wp) t1
        real (wp) t2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) xh(*)
        real (wp) xhold

        ier = 0

        if (n < 2) then
            return
        else if (n == 2) then
            xhold = HALF_SQRT3*(x(1,1)+x(1,2))
            x(1,2) = HALF_SQRT3*(x(1,1)-x(1,2))
            x(1,1) = xhold
        else
            np1 = n+1
            ns2 = n/2
            do k=1,ns2
                kc = np1-k
                t1 = x(1,k)-x(1,kc)
                t2 = wsave(k)*(x(1,k)+x(1,kc))
                xh(k+1) = t1+t2
                xh(kc+1) = t2-t1
            end do

            modn = mod(n,2)

            if (modn /= 0) then
                xh(ns2+2) =  4.0_wp * x(1,ns2+1)
            end if

            xh(1) = 0.0_wp
            lnxh = np1
            lnsv = np1 + int(log(real(np1, kind=wp))/log(2.0_wp)) + 4
            lnwk = np1

            call rfft1f(np1,1,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sintb1',-5)
                return
            end if

            if (mod(np1,2) == 0) then
                xh(np1) = xh(np1)+xh(np1)
            end if

            fnp1s4 = real(np1, kind=wp)/4
            x(1,1) = fnp1s4*xh(1)
            dsum = x(1,1)

            do i=3,n,2
                x(1,i-1) = fnp1s4*xh(i)
                dsum = dsum+fnp1s4*xh(i-1)
                x(1,i) = dsum
            end do

            if (modn == 0) then
                x(1,n) = fnp1s4*xh(n+1)
            end if
        end if

    end subroutine sintb1




    subroutine sintf1(n, inc, x, wsave, xh, work, ier)

        integer (ip) inc
        real (wp) dsum
        integer (ip) i
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) kc
        integer (ip) lnsv
        integer (ip) lnwk
        integer (ip) lnxh
        integer (ip) modn
        integer (ip) n
        integer (ip) np1
        integer (ip) ns2
        real (wp) sfnp1
        real (wp) t1
        real (wp) t2
        real (wp) work(*)
        real (wp) wsave(*)
        real (wp) x(inc,*)
        real (wp) xh(*)
        real (wp) xhold

        ier = 0

        if (n < 2) then
            return
        else if (n == 2) then
            xhold = (x(1,1)+x(1,2))/sqrt(3.0_wp)
            x(1,2) = (x(1,1)-x(1,2))/sqrt(3.0_wp)
            x(1,1) = xhold
        else
            np1 = n+1
            ns2 = n/2
            do k=1,ns2
                kc = np1-k
                t1 = x(1,k)-x(1,kc)
                t2 = wsave(k)*(x(1,k)+x(1,kc))
                xh(k+1) = t1+t2
                xh(kc+1) = t2-t1
            end do

            modn = mod(n,2)

            if (modn /= 0) then
                xh(ns2+2) =  4.0_wp * x(1,ns2+1)
            end if

            xh(1) = 0.0_wp
            lnxh = np1
            lnsv = np1 + int(log(real(np1, kind=wp))/log(2.0_wp)) + 4
            lnwk = np1

            call rfft1f(np1,1,xh,lnxh,wsave(ns2+1),lnsv,work, lnwk,local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sintf1',-5)
                return
            end if

            if (mod(np1,2) == 0) then
                xh(np1) = xh(np1)+xh(np1)
            end if

            sfnp1 = 1.0_wp/np1
            x(1,1) = 0.5_wp * xh(1)
            dsum = x(1,1)

            do i=3,n,2
                x(1,i-1) = 0.5_wp * xh(i)
                dsum = dsum + 0.5_wp * xh(i-1)
                x(1,i) = dsum
            end do

            if (modn == 0) then
                x(1,n) = 0.5_wp * xh(n+1)
            end if
        end if

    end subroutine sintf1


    subroutine sintmb(lot, jump, n, inc, x, lenx, wsave, lensav, &
        work, lenwrk, ier)
        !
        ! SINTMB: 64-bit float precision backward sine transform, multiple vectors.
        !
        !  Purpose:
        !
        !  SINTMB computes the one-dimensional Fourier transform of multiple
        !  odd sequences within a real array.  This transform is referred to as
        !  the backward transform or Fourier synthesis, transforming the
        !  sequences from spectral to physical space.
        !
        !  This transform is normalized since a call to SINTMB followed
        !  by a call to SINTMF (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  Parameters:
        !
        !  integer LOT, the number of sequences to be transformed
        !  within the array R.
        !
        !  integer JUMP, the increment between the locations, in
        !  array R, of the first elements of two consecutive sequences.
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N+1 is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations, in
        !  array R, of two consecutive elements within the same sequence.
        !
        !  Input/real R(LENR), containing LOT sequences, each
        !  having length N.  R can have any number of dimensions, but the total
        !  number of locations must be at least LENR.  On input, R contains the data
        !  to be transformed, and on output, the transformed data.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least (LOT-1)*JUMP + INC*(N-1)+ 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to SINTMI before the first call to routine SINTMF
        !  or SINTMB for a given transform length N.  wsave's contents may be re-used
        !  for subsequent calls to SINTMF and SINTMB with the same N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least N/2 + N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least LOT*(2*N+4).
        !
        !  integer IER, error flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  4, input parameters INC,JUMP,N,LOT are not consistent;
        !  20, input error returned by lower level routine.
        !


        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) iw1
        integer (ip) iw2
        integer (ip) jump
        integer (ip) lenx
        integer (ip) lot
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        if (lenx < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('sintmb', 6)
            return
        else if ( lensav < n / 2 + n + int(log(real(n, kind=wp) ) &
            /log(2.0_wp)) +4) then
            ier = 2
            call fft_error_handler('sintmb', 8)
            return
        else if (lenwrk < lot*(2*n+4)) then
            ier = 3
            call fft_error_handler('sintmb', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('sintmb', -1)
            return
        else
            ier = 0
        end if

        iw1 = 2*lot+1
        iw2 = iw1+lot*(n+1)

        call msntb1(lot,jump,n,inc,x,wsave,work,work(iw1),work(iw2),local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('sintmb',-5)
            return
        end if

    end subroutine sintmb


    subroutine sintmf(lot, jump, n, inc, x, lenx, wsave, lensav, &
        work, lenwrk, ier)
        !
        !  sintmf: forward sine transform, multiple vectors.
        !
        !  Purpose:
        !
        !  Computes the 1-dimensional fourier transform of multiple
        !  odd sequences within a real array. This transform is referred to as
        !  the forward transform or fourier analysis, transforming the sequences
        !  from physical to spectral space.
        !
        !  This transform is normalized since a call to sintmf followed
        !  by a call to sintmb (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  parameters:
        !
        !  integer lot, the number of sequences to be
        !  transformed within.
        !
        !  integer jump, the increment between the locations,
        !  in array r, of the first elements of two consecutive sequences.
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n+1 is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array r, of two consecutive elements within the same sequence.
        !
        !  input/real r(lenr), containing lot sequences, each
        !  having length n.  r can have any number of dimensions, but the total
        !  number of locations must be at least lenr.  on input, r contains the data
        !  to be transformed, and on output, the transformed data.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least (lot-1)*jump + inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to sintmi before the first call to routine sintmf
        !  or sintmb for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to sintmf and sintmb with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n/2 + n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least lot*(2*n+4).
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  4, input parameters inc,jump,n,lot are not consistent;
        !  20, input error returned by lower level routine.
        !
        !--------------------------------------------------------------
        ! Dictionary: Calling arguments
        !--------------------------------------------------------------
        integer (ip) inc
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) iw1
        integer (ip) iw2
        integer (ip) jump
        integer (ip) lenx
        integer (ip) lot
        integer (ip) n
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) x(inc,*)

        !
        !==> Check validity of input arguments
        !
        if ( lenx < ( lot - 1) * jump + inc * ( n - 1 ) + 1 ) then
            ier = 1
            call fft_error_handler('sintmf', 6)
            return
        else if ( lensav < n / 2 + n + int(log(real(n, kind=wp) ) &
            / log(2.0_wp ) ) + 4 ) then
            ier = 2
            call fft_error_handler('sintmf', 8)
            return
        else if ( lenwrk < lot * ( 2 * n + 4 ) ) then
            ier = 3
            call fft_error_handler('sintmf', 10)
            return
        else if (.not. fft_consistent(inc, jump, n, lot)) then
            ier = 4
            call fft_error_handler('sintmf', -1)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        iw1 = 2 * lot + 1
        iw2 = iw1 + lot * (n + 1)

        call msntf1(lot, jump, n, inc, x, wsave, work, work(iw1), work(iw2), local_error_flag)

        ! Check error flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('sintmf', -5)
        end if

    end subroutine sintmf


    subroutine sintmi(n, wsave, lensav, ier)
        !
        ! sintmi: initialization for sintmb and sintmf.
        !
        !  Purpose:
        !
        !  Initializes array wsave for use in its companion routines
        !  sintmf and sintmb.  The prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  Separate wsave arrays are required for different
        !  values of n.
        !
        !  Parameters:
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n/2 + n + int(log(real(n))) + 4.
        !
        !  real wsave(lensav), containing the prime factors
        !  of n and also containing certain trigonometric values which will be used
        !  in routines sintmb or sintmf.
        !
        !  integer ier, error flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  20, input error returned by lower level routine.
        !
        !--------------------------------------------------------------
        ! Dictionary: Calling arguments
        !--------------------------------------------------------------
        integer (ip) lensav
        real (wp) dt
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) k
        integer (ip) lnsv
        integer (ip) n
        integer (ip) np1
        integer (ip) ns2
        real (wp), parameter :: PI = acos(-1.0_wp)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lensav < n / 2 + n + int(log(real(n, kind=wp)) &
            / log(2.0_wp ) ) + 4 ) then
            ier = 2
            call fft_error_handler('sintmi', 3 )
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n > 1) then

            ns2 = n / 2
            np1 = n + 1
            dt = PI/np1

            do k = 1, ns2
                wsave(k) = 2.0_wp * sin(real(k, kind=wp) * dt)
            end do

            associate( nsv => np1+int(log(real(np1, kind=wp) )/log(2.0_wp))+4)

                call rfftmi(np1, wsave(ns2+1), lnsv, local_error_flag)

            end associate

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sintmi', -5)
            end if
        end if

    end subroutine sintmi


    subroutine w2r(ldr, ldw, l, m, r, w)
        !
        ! Purpose:
        !
        ! Copies a 2D array, allowing for different leading dimensions.
        !
        !--------------------------------------------------------------
        ! Dictionary: Calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)     :: ldr
        integer (ip), intent (in)     :: ldw
        integer (ip), intent (in)     :: m
        integer (ip), intent (in)     :: l
        real (wp),    intent (in out) :: r(ldr,m)
        real (wp),    intent (in)     :: w(ldw,m)
        !--------------------------------------------------------------

        r(1:l,:) = w(1:l,:)

    end subroutine w2r


    pure function fft_consistent(inc, jump, n, lot) result (return_value)
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
        ! Dictionary: Calling arguments
        !--------------------------------------------------------------
        integer (ip), intent (in) :: inc
        integer (ip), intent (in) :: jump
        integer (ip), intent (in) :: n
        integer (ip), intent (in) :: lot
        logical                   :: return_value
        !--------------------------------------------------------------
        ! Dictionary: Calling arguments
        !--------------------------------------------------------------
        integer (ip) :: lcm
        !--------------------------------------------------------------

        !
        !==> Set least common multiple of inc and jump.
        !
        lcm = get_least_common_multiple(inc, jump)

        !
        !==> Check consistency
        !
        if ( lcm <= (n - 1) * inc .and. lcm <= (lot - 1) * jump ) then
            return_value = .false.
        else
            return_value = .true.
        end if

    contains

        pure function get_least_common_multiple(a, b) result (return_value)
            !--------------------------------------------------------------
            ! Dictionary: Calling arguments
            !--------------------------------------------------------------
            integer (ip), intent (in) :: a
            integer (ip), intent (in) :: b
            integer (ip)              :: return_value
            !--------------------------------------------------------------
            ! Dictionary: local variables
            !--------------------------------------------------------------
            integer (ip) :: i, j, jnew !! Counters
            !--------------------------------------------------------------

            i = a
            j = b

            do while (j /= 0)
                jnew = mod(i, j)
                i = j
                j = jnew
            end do

            ! Return least common multiple of a and b
            return_value = (a * b) / i

        end function get_least_common_multiple

    end function fft_consistent


    subroutine fft_error_handler(subroutine_name, info)
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
        ! Dictionary: calling arguments
        !--------------------------------------------------------------
        integer (ip),      intent (in) :: info
        character (len=*), intent (in) :: subroutine_name
        !--------------------------------------------------------------

        write ( stderr, '(A)' ) ''
        write ( stderr, '(A)' ) ' Object of class (FFTpack): '
        write ( stderr, '(A)' ) ' FATAL ERROR '

        select case (info)
            case (-1)
                write( stderr, '(4(A))') '  On entry to ', &
                    trim(subroutine_name), &
                    ' parameters lot, jump, n and inc are inconsistent.'
            case (-2)
                write( stderr, '(4(A))')  '  On ntry to ', &
                    trim(subroutine_name), &
                    ' parameter l is greater than ldim.'
            case (-3)
                write( stderr, '(4(A))')  '  On entry to ', &
                    trim(subroutine_name), &
                    ' parameter m is greater than mdim.'
            case (-5)
                write( stderr, '(4(A))')  '  Within ', &
                    trim(subroutine_name), &
                    ' input error returned by lower level routine.'
            case (-6)
                write( stderr, '(4(A))')  '  On entry to ', &
                    trim(subroutine_name), &
                    ' parameter ldim is less than 2*(l/2+1).'
            case default
                write ( stderr, '(3(A),I3,A)') '  On entry to ',&
                    trim(subroutine_name),&
                    ' parameter number ', info, ' had an illegal value.'
        end select

    end subroutine fft_error_handler

end module type_FFTpack
