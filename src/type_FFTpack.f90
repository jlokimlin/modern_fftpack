module type_FFTpack

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use auxiliary_routines

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use complex_transform_routines

    use real_transform_routines

    use cosine_transform_routines

    use sine_transform_routines

    use quarter_cosine_transform_routines

    use quarter_sine_transform_routines

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: FFTpack


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
        procedure, private :: real_1d_forward
        procedure, private :: real_1d_backward
        procedure, private :: real_2d_forward
        procedure, private :: real_2d_backward
        procedure, private :: real_nd_forward
        procedure, private :: real_nd_backward
        procedure, private :: complex_1d_forward
        procedure, private :: complex_1d_backward
        procedure, private :: complex_2d_forward
        procedure, private :: complex_2d_backward
        !procedure, private :: complex_nd_forward
        !procedure, private :: complex_nd_backward
        procedure, private :: cost_1d_forward
        procedure, private :: cost_1d_backward

        generic,   public  :: fft => real_1d_forward, complex_1d_forward
        generic,   public  :: ifft => real_1d_backward, complex_1d_backward
        generic,   public  :: dct => cost_1d_forward
        generic,   public  :: idct => cost_1d_backward
        generic,   public  :: fft2 => real_2d_forward, complex_2d_forward
        generic,   public  :: ifft2 => real_2d_backward, complex_2d_backward
        generic,   public  :: fftn => real_nd_forward!, complex_nd_forward
        generic,   public  :: ifftn => real_nd_backward!, complex_nd_backward

        procedure, public  :: destroy => destroy_fftpack
        procedure, public  :: create => create_fftpack


        procedure, nopass, public :: get_1d_saved_workspace_length
        procedure, nopass, public :: get_1d_saved_workspace
        procedure, nopass, public :: get_1d_workspace_length
        procedure, nopass, public :: get_1d_workspace
        procedure, nopass, public :: get_1d_sin_workspace_length
        procedure, nopass, public :: get_1d_sin_workspace
        procedure, nopass, public :: get_2d_saved_workspace_length
        procedure, nopass, public :: get_2d_saved_workspace
        procedure, nopass, public :: get_2d_workspace_length
        procedure, nopass, public :: get_2d_workspace

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


    subroutine real_1d_forward(this, real_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        real(wp),        intent(in out)  :: real_data(:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: n, error_flag
        !-----------------------------------------------------------------

        n = size(real_data)

        associate( &
            r => real_data, &
            lenwrk => get_real_1d_workspace_length(n), &
            lensav => get_real_1d_saved_workspace_length(n) , &
            lenr => size(real_data), &
            inc => 1, &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize real transform
                !
                call rfft1i(n,wsave,lensav,ier)

                ! Check error_flag
                if(ier == 2) then
                    error stop "rfft1i: lensave not big enough"
                end if

                !
                !==> Perform real forward transform
                !
                call rfft1f(n,inc,r,lenr,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error_flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop "rfft1f: lenr not big enough"
            case (2)
                error stop "rfft1f: lensav not big enough"
            case (3)
                error stop "rfft1f: lenwrk not big enough"
            case (20)
                error stop "rfft1f: input error returned by lower level routine"
        end select

    end subroutine real_1d_forward



    subroutine real_1d_backward(this, real_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        real (wp),       intent(in out)  :: real_data(:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: n, error_flag
        !-----------------------------------------------------------------

        n = size(real_data)

        associate( &
            r => real_data, &
            lenwrk => get_real_1d_workspace_length(n), &
            lensav => get_real_1d_saved_workspace_length(n) , &
            lenr => size(real_data), &
            inc => 1, &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize real transform
                !
                call rfft1i(n,wsave,lensav,ier)

                ! Check error_flag
                if(ier == 2) then
                    error stop "rfft1i: lensave not big enough"
                end if

                call rfft1b(n,inc,r,lenr,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error_flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop "rfft1b: lenr not big enough"
            case (2)
                error stop "rfft1b: lensav not big enough"
            case (3)
                error stop "rfft1b: lenwrk not big enough"
            case (20)
                error stop "rfft1b: input error returned by lower level routine"
        end select

    end subroutine real_1d_backward



    subroutine complex_1d_forward(this, complex_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        complex (wp),    intent (in out) :: complex_data(:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: n, error_flag
        !-----------------------------------------------------------------

        n = size(complex_data)

        associate( &
            c => complex_data, &
            lenwrk => get_complex_1d_workspace_length(n), &
            lensav => get_complex_1d_saved_workspace_length(n), &
            lenc => n, &
            inc => 1, &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize complex transform
                !
                call cfft1i(n,wsave,lensav,ier)

                ! Check error_flag
                if(ier == 2) then
                    error stop "cfft1i: lensav not big enough"
                end if

                !
                !==> Perform complex forward transform
                !
                call cfft1f(n,inc,c,lenc,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        select case (error_flag)
            case (0)
                return
            case (1)
                error stop "cfft1f: lenc not big enough"
            case (2)
                error stop "cfft1f: lensav not big enough"
            case (3)
                error stop "cfft1f: lenwrk not big enough"
            case (20)
                error stop "cfft1f: input error returned by lower level routine"
        end select

    end subroutine complex_1d_forward



    subroutine complex_1d_backward(this, complex_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        complex (wp),    intent (in out) :: complex_data(:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: n, error_flag
        !-----------------------------------------------------------------

        n = size(complex_data)

        associate( &
            c => complex_data, &
            lenwrk => get_complex_1d_workspace_length(n), &
            lensav => get_complex_1d_saved_workspace_length(n), &
            lenc => n, &
            inc => 1, &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize complex transform
                !
                call cfft1i(n,wsave,lensav,ier)

                ! Check error_flag
                if(ier == 2) then
                    error stop "cfft1i: lensav not big enough"
                end if

                !
                !==> Perform complex backward transform
                !
                call cfft1b(n,inc,c,lenc,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        select case (error_flag)
            case (0)
                return
            case (1)
                error stop "cfft1b: lenc not big enough"
            case (2)
                error stop "cfft1b: lensav not big enough"
            case (3)
                error stop "cfft1b: lenwrk not big enough"
            case (20)
                error stop "cfft1b: input error returned by lower level routine"
        end select

    end subroutine complex_1d_backward


    subroutine cost_1d_forward(this, real_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        real (wp),       intent(in out)  :: real_data(:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: n, error_flag
        !-----------------------------------------------------------------

        n = size(real_data)

        associate( &
            r => real_data, &
            lenwrk => get_cost_1d_workspace_length(n), &
            lensav => get_cost_1d_saved_workspace_length(n) , &
            lenr => n, &
            inc => 1, &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize transform
                !
                call cost1i(n, wsave, lensav, ier)

                ! Check error_flag
                if (ier == 2) then
                    error stop "cost_1d_backward: lensave not big enough"
                else if (ier == 20) then
                    error stop "cost_1d_backward: input error returned by lower level routine"
                end if

                !
                !==> Perform transform
                !
                call cost1f(n, inc, r, lenr, wsave, lensav, work, lenwrk, ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error_flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop "cost1f: lenr not big enough"
            case (2)
                error stop "cost1f: lensav not big enough"
            case (3)
                error stop "cost1f: lenwrk not big enough"
            case (20)
                error stop "cost1f: input error returned by lower level routine"
        end select

    end subroutine cost_1d_forward


    subroutine cost_1d_backward(this, real_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        real (wp),       intent(in out)  :: real_data(:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: n, error_flag
        !-----------------------------------------------------------------

        n = size(real_data)

        associate( &
            r => real_data, &
            lenwrk => get_cost_1d_workspace_length(n), &
            lensav => get_cost_1d_saved_workspace_length(n) , &
            lenr => n, &
            inc => 1, &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize transform
                !
                call cost1i(n, wsave, lensav, ier)

                ! Check error_flag
                if (ier == 2) then
                    error stop "cost_1d_backward: lensave not big enough"
                else if (ier == 20) then
                    error stop "cost_1d_backward: input error returned by lower level routine"
                end if

                !
                !==> Perform transform
                !
                call cost1b(n, inc, r, lenr, wsave, lensav, work, lenwrk, ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error_flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop "cost_1d_backward: lenr not big enough"
            case (2)
                error stop "cost_1d_backward: lensav not big enough"
            case (3)
                error stop "cost_1d_backward: lenwrk not big enough"
            case (20)
                error stop "cost_1d_backward: input error returned by lower level routine"
        end select

    end subroutine cost_1d_backward


    subroutine real_2d_forward(this, real_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        real (wp),       intent (in out) :: real_data(:,:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: l, m, error_flag
        !-----------------------------------------------------------------

        l = size(real_data, dim=1)
        m = size(real_data, dim=2)

        associate( &
            ldim => l, &
            r => real_data, &
            lensav => get_real_2d_saved_workspace_length(l, m), &
            lenwrk => get_real_2d_workspace_length(l, m), &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize transform
                !
                call rfft2i(l,m,wsave,lensav,ier)

                ! Check error flag
                select case (ier)
                    case (2)
                        error stop "rfft2i: lensav not big enough"
                    case (20)
                        error stop "rfft2i: input error returned by lower level routine"
                end select

                !
                !==> Perform transform
                !
                call rfft2f(ldim,l,m,real_data,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error flag
        select case (error_flag)
            case (0)
                return
            case (6)
                error stop "rfft2f: ldim is less than 2*((l+1)/2)"
            case (2)
                error stop "rfft2f: lensav not big enough"
            case (3)
                error stop "rfft2f: lenwrk not big enough"
            case (20)
                error stop "rfft2f: input error returned by lower level routine"
        end select

    end subroutine real_2d_forward


    subroutine real_2d_backward(this, real_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        real (wp),       intent (in out) :: real_data(:,:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: l, m, error_flag
        !-----------------------------------------------------------------

        l = size(real_data, dim=1)
        m = size(real_data, dim=2)

        associate( &
            ldim => l, &
            r => real_data, &
            lensav => get_real_2d_saved_workspace_length(l, m), &
            lenwrk => get_real_2d_workspace_length(l, m), &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize transform
                !
                call rfft2i(l,m,wsave,lensav,ier)

                ! Check error flag
                select case (ier)
                    case (2)
                        error stop "rfft2i: lensav not big enough"
                    case (20)
                        error stop "rfft2i: input error returned by lower level routine"
                end select

                !
                !==> perform transform
                call rfft2b(ldim,l,m,real_data,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error flag
        select case(error_flag)
            case (0)
                return
            case (6)
                error stop "rfft2b: ldim is less than 2*((l+1)/2)"
            case (2)
                error stop "rfft2b: lensav not big enough"
            case (3)
                error stop "rfft2b: lenwrk not big enough"
            case (20)
                error stop "rfft2b: input error returned by lower level routine"
        end select

    end subroutine real_2d_backward


    subroutine complex_2d_forward(this, complex_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        complex (wp),    intent (in out) :: complex_data(:,:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: l, m, error_flag
        !-----------------------------------------------------------------

        l = size(complex_data, dim=1)
        m = size(complex_data, dim=2)

        associate( &
            c => complex_data, &
            ldim => l, &
            lensav => get_complex_2d_saved_workspace_length(l, m), &
            lenwrk => get_complex_2d_workspace_length(l, m), &
            ier => error_flag &
            )

                        !
                        !==> Allocate memory
                        !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize transform
                !
                call cfft2i(l,m,wsave,lensav,ier)

                select case (ier)
                    case (2)
                        error stop "cfft2i: lensav not big enough"
                    case (20)
                        error stop "cfft2i: input error returned by lower level routine"
                end select

                !
                !==> Perform transform
                !
                call cfft2f(ldim,l,m,c,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate
        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error flag
        select case(error_flag)
            case (0)
                return
            case (5)
                error stop "cfft2f: l > ldim"
            case (2)
                error stop "cfft2f: lensav not big enough"
            case (3)
                error stop "cfft2f: lenwrk not big enough"
            case (20)
                error stop "cfft2f: input error returned by lower level routine"
        end select

    end subroutine complex_2d_forward

    subroutine complex_2d_backward(this, complex_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        complex (wp),    intent (in out) :: complex_data(:,:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: l, m, error_flag
        !-----------------------------------------------------------------

        l = size(complex_data, dim=1)
        m = size(complex_data, dim=2)

        associate( &
            c => complex_data, &
            ldim => l, &
            lensav => get_complex_2d_saved_workspace_length(l, m), &
            lenwrk => get_complex_2d_workspace_length(l, m), &
            ier => error_flag &
            )

             !
             !==> Allocate memory
             !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize transform
                !
                call cfft2i(l,m,wsave,lensav,ier)

                select case (ier)
                    case (2)
                        error stop "cfft2i: lensav not big enough"
                    case (20)
                        error stop "cfft2i: input error returned by lower level routine"
                end select

                !
                !==> Perform transform
                !
                call cfft2b(ldim, l, m, c, wsave, lensav, work, lenwrk, ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error flag
        select case (error_flag)
            case (0)
                return
            case (5)
                error stop "cfft2b: l > ldim"
            case (2)
                error stop "cfft2b: lensav not big enough"
            case (3)
                error stop "cfft2b: lenwrk not big enough"
            case (20)
                error stop "cfft2b: input error returned by lower level routine"
        end select

    end subroutine complex_2d_backward


    subroutine real_nd_forward(this, n, lot, real_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        integer (ip),    intent (in)     :: n   ! length of each sequence
        integer (ip),    intent (in)     :: lot ! number of sequences
        real(wp),        intent(in out)  :: real_data(:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: l, error_flag
        !-----------------------------------------------------------------

        l = size(real_data)

        if (mod(n*lot,l)/=0) then
            error stop "real_nd_forward: "&
                //"incommensurate input arguments "&
                //"mod(n*lot,size(real_data)) /= 0."
        end if


        associate( &
            r => real_data, &
            lenwrk => get_real_nd_workspace_length(n, lot), &
            lensav => get_real_nd_saved_workspace_length(n) , &
            lenr => n*lot, &
            jump => n, &
            inc => 1, &
            ier => error_flag &
            )

             !
             !==> Allocate memory
             !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize transform
                !
                call rfftmi(n, wsave, lensav, ier)

                ! Check error_flag
                if(ier == 2) then
                    error stop "rfftmi: lensave not big enough"
                end if

                !
                !==> Perform forward transform
                !
                call rfftmf(lot,jump,n,inc,r,lenr,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error_flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop "rfftmf: lenr not big enough"
            case (2)
                error stop "rfftmf: lensav not big enough"
            case (3)
                error stop "rfftmf: lenwrk not big enough"
            case (4)
                error stop "rfftmf: inc, jump, n, lot are not consistent"
        end select

    end subroutine real_nd_forward

    subroutine real_nd_backward(this, n, lot, real_data)
        !-----------------------------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------------------------
        class (FFTpack), intent (in out) :: this
        integer (ip),    intent (in)     :: n   ! length of each sequence
        integer (ip),    intent (in)     :: lot ! number of sequences
        real (wp),       intent(in out)  :: real_data(:)
        !-----------------------------------------------------------------
        ! Local variables
        !-----------------------------------------------------------------
        integer (ip) :: l, error_flag
        !-----------------------------------------------------------------

        l = size(real_data)

        if (mod(n*lot,l)/=0) then
            error stop "real_nd_backward: "&
                //"incommensurate input arguments "&
                //"mod(n*lot,size(real_data)) /= 0."
        end if


        associate( &
            r => real_data, &
            lenwrk => get_real_nd_workspace_length(n, lot), &
            lensav => get_real_nd_saved_workspace_length(n) , &
            lenr => n*lot, &
            jump => n, &
            inc => 1, &
            ier => error_flag &
            )

            !
            !==> Allocate memory
            !
            call this%create(lensav, lenwrk)

            associate( &
                wsave => this%saved_workspace, &
                work => this%workspace &
                )

                !
                !==> Initialize transform
                !
                call rfftmi(n, wsave, lensav, ier)

                ! Check error_flag
                if(ier == 2) then
                    error stop "rfftmi: lensave not big enough"
                end if

                !
                !==> Perform forward transform
                !
                call rfftmb(lot,jump,n,inc,r,lenr,wsave,lensav,work,lenwrk,ier)

            end associate
        end associate

        !
        !==> Release memory
        !
        call this%destroy()

        ! Check error_flag
        select case (error_flag)
            case (0)
                return
            case (1)
                error stop "rfftmb: lenr not big enough"
            case (2)
                error stop "rfftmb: lensav not big enough"
            case (3)
                error stop "rfftmb: lenwrk not big enough"
            case (4)
                error stop "rfftmb: inc, jump, n, lot are not consistent"
        end select

    end subroutine real_nd_backward


    !    subroutine complex_nd_forward(this, n, lot, complex_data)
    !        !-----------------------------------------------------------------
    !        ! Dummy arguments
    !        !-----------------------------------------------------------------
    !        class (FFTpack), intent (in out) :: this
    !        integer (ip),    intent (in)     :: n   ! length of each sequence
    !        integer (ip),    intent (in)     :: lot ! number of sequences
    !        complex (wp),    intent(in out)  :: complex_data(:)
    !        !-----------------------------------------------------------------
    !        ! Local variables
    !        !-----------------------------------------------------------------
    !        integer (ip) :: l, error_flag
    !        !-----------------------------------------------------------------
    !
    !        l = size(complex_data)
    !
    !        if (mod(n*lot,l)/=0) then
    !            error stop "complex_nd_forward: "&
    !                //"incommensurate input arguments "&
    !                //"mod(n*lot,size(complex_data)) /= 0."
    !        end if
    !
    !
    !        associate( &
    !            c => complex_data, &
    !            lenwrk => get_complex_nd_workspace_length(n, lot), &
    !            lensav => get_complex_nd_saved_workspace_length(n) , &
    !            lenc => n*lot, &
    !            jump => n, &
    !            inc => 1, &
    !            ier => error_flag &
    !            )
    !
    !            !
    !            !==> Allocate memory
    !            !
    !            allocate( this%saved_workspace(lensav) )
    !            allocate( this%workspace(lenwrk) )
    !
    !            associate( &
    !                wsave => this%saved_workspace, &
    !                work => this%workspace &
    !                )
    !
    !                !
    !                !==> Initialize transform
    !                !
    !                call cfftmi(n, wsave, lensav, ier)
    !
    !                ! Check error_flag
    !                if(ier == 2) then
    !                    error stop "cfftmi: lensave not big enough"
    !                end if
    !
    !                 !
    !                 !==> Perform transform
    !                 !
    !                call cfftmf(lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
    !
    !            end associate
    !        end associate
    !
    !        !
    !        !==> Release memory
    !        !
    !        call this%destroy()
    !
    !        ! Check error_flag
    !        select case (error_flag)
    !            case (0)
    !                return
    !            case (1)
    !                error stop "cfftmf: lenc not big enough"
    !            case (2)
    !                error stop "cfftmf: lensav not big enough"
    !            case (3)
    !                error stop "cfftmf: lenwrk not big enough"
    !            case (4)
    !                error stop "cfftmf: inc, jump, n, lot are not consistent"
    !        end select
    !
    !    end subroutine complex_nd_forward
    !
    !    subroutine complex_nd_backward(this, n, lot, complex_data)
    !        !-----------------------------------------------------------------
    !        ! Dummy arguments
    !        !-----------------------------------------------------------------
    !        class (FFTpack), intent (in out) :: this
    !        integer (ip),    intent (in)     :: n   ! length of each sequence
    !        integer (ip),    intent (in)     :: lot ! number of sequences
    !        complex (wp),    intent(in out)  :: complex_data(:)
    !        !-----------------------------------------------------------------
    !        ! Local variables
    !        !-----------------------------------------------------------------
    !        integer (ip) :: l, error_flag
    !        !-----------------------------------------------------------------
    !
    !        l = size(complex_data)
    !
    !        if (mod(n*lot,l)/=0) then
    !            error stop "complex_nd_backward: "&
    !                //"incommensurate input arguments "&
    !                //"mod(n*lot,size(complex_data)) /= 0."
    !        end if
    !
    !
    !        associate( &
    !            c => complex_data, &
    !            lenwrk => get_complex_nd_workspace_length(n, lot), &
    !            lensav => get_complex_nd_saved_workspace_length(n) , &
    !            lenc => n*lot, &
    !            jump => n, &
    !            inc => 1, &
    !            ier => error_flag &
    !            )
    !
    !            !
    !            !==> Allocate memory
    !            !
    !            allocate( this%saved_workspace(lensav) )
    !            allocate( this%workspace(lenwrk) )
    !
    !            associate( &
    !                wsave => this%saved_workspace, &
    !                work => this%workspace &
    !                )
    !
    !                !
    !                !==> Initialize transform
    !                !
    !                call cfftmi(n, wsave, lensav, ier)
    !
    !                ! Check error_flag
    !                if(ier == 2) then
    !                    error stop "cfftmi: lensave not big enough"
    !                end if
    !
    !                !
    !                !==> Perform transform
    !                !
    !                call cfftmb(lot,jump,n,inc,c,lenc,wsave,lensav,work,lenwrk,ier)
    !
    !            end associate
    !        end associate
    !
    !        !
    !        !==> Release memory
    !        !
    !        call this%destroy()
    !
    !        ! Check error_flag
    !        select case (error_flag)
    !            case (0)
    !                return
    !            case (1)
    !                error stop "cfftmb: lenc not big enough"
    !            case (2)
    !                error stop "cfftmb: lensav not big enough"
    !            case (3)
    !                error stop "cfftmb: lenwrk not big enough"
    !            case (4)
    !                error stop "cfftmb: inc, jump, n, lot are not consistent"
    !        end select
    !
    !    end subroutine complex_nd_backward

    pure function fftpack_2d_constructor(l, m) result (return_value)
        !------------------------------------------------------------------
        ! Dummy arguments
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
        ! Dummy arguments
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



    subroutine create_fftpack(this, lensav, lenwrk)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        class (FFTpack) , intent (in out) :: this
        integer (ip),     intent (in)     :: lensav
        integer (ip),     intent (in)     :: lenwrk
        !------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy()

        !
        !==> Allocate memory
        !
        allocate( this%saved_workspace(lensav) )
        allocate( this%workspace(lenwrk) )

    end subroutine create_fftpack



    subroutine destroy_fftpack(this)
        !------------------------------------------------------------------
        ! Dummy arguments
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




end module type_FFTpack
