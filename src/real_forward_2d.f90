submodule (real_transform_routines:multiple_real_forward) real_forward_2d

contains

    module subroutine rfft2f(ldim, l, m, r, wsave, lensav, work, lenwrk, ierror)
        !
        ! rfft2f: real forward fast fourier transform, 2d.
        !
        ! purpose:
        !
        !  rfft2f computes the two-dimensional discrete fourier transform of a
        !  real periodic array.  this transform is known as the forward transform
        !  or fourier analysis, transforming from physical to spectral space.
        !  routine rfft2f is normalized: a call to rfft2f followed by a call to
        !  rfft2b (or vice-versa) reproduces the original array within roundoff
        !  error.
        !
        !  parameters:
        !
        !  integer ldim, the first dimension of the 2d real
        !  array r, which must be at least 2*(l/2+1).
        !
        !  integer l, the number of elements to be transformed
        !  in the first dimension of the two-dimensional real array r.  the value
        !  of l must be less than or equal to that of ldim.  the transform is most
        !  efficient when l is a product of small primes.
        !
        !  integer m, the number of elements to be transformed
        !  in the second dimension of the two-dimensional real array r.  the
        !  transform is most efficient when m is a product of small primes.
        !
        !  input/real r(ldim,m), the real array of two
        !  dimensions.  on input, containing the l-by-m physical data to be
        !  transformed.  on output, the spectral coefficients.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to rfft2i before the first call to routine rfft2f
        !  or rfft2b with lengths l and m.  wsave's contents may be re-used for
        !  subsequent calls to rfft2f and rfft2b with the same transform lengths.
        !
        !  integer lensav, the number of elements in the wsave
        !  array.  lensav must be at least l + m + int(log(real(l)))
        !  + int(log(real(m))) + 8.
        !
        !  workspace, real (wp) work(lenwrk), provides workspace, and its
        !  contents need not be saved between calls to routines rfft2f and rfft2b.
        !
        !  integer lenwrk, the number of elements in the work
        !  array.  lenwrk must be at least ldim*m.
        !
        !  integer ierror, the error_flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  6, input parameter ldim < 2*(l+1);
        !  20, input error returned by lower level routine.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)     :: ldim
        integer (ip), intent (in)     :: l
        integer (ip), intent (in)     :: m
        real (wp),    intent (in out) :: r(ldim,m)
        real (wp),    intent (in)     :: wsave(lensav)
        integer (ip), intent (in)     :: lensav
        real (wp),    intent (out)    :: work(lenwrk)
        integer (ip), intent (in)     :: lenwrk
        integer (ip), intent (out)    :: ierror
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: local_error_flag
        integer (ip) :: ldh, ldw, ldx, lwsav
        integer (ip) :: mmsav, modl, modm, mwsav
        !--------------------------------------------------------------

        ! initialize error flag
        ierror = 0
        !
        !==> verify lensav
        !
        lwsav = l+int(log(real (l, kind=wp))/log(TWO))+4
        mwsav = get_1d_saved_workspace_length(m)
        mmsav = m+int(log(real(m, kind=wp))/log(TWO))+4

        if (lensav < lwsav+mwsav+mmsav) then
            ierror = 2
            call fft_error_handler('rfft2f', 6)
            return
        end if
        !
        !==>  verify lenwrk
        !
        if (lenwrk < (l+1)*m) then
            ierror = 3
            call fft_error_handler('rfft2f', 8)
            return
        end if
        !
        !==>  verify ldim is as big as l
        !
        if (ldim < l) then
            ierror = 5
            call fft_error_handler('rfft2f', -6)
            return
        end if
        !
        !==>  Transform first dimension of array
        !
        associate( &
            arg_1 => m*ldim, &
            arg_2 => l+int(log(real(l, kind=wp))/log(TWO))+4 &
            )

            call rfftmf(m,ldim,l,1,r,arg_1,wsave(1), arg_2,work,size(work),local_error_flag)

        end associate

        if (local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('rfft2f',-5)
            return
        end if

        ldx = 2*int((l+1)/2)-1

        r(2:ldx,1:m) = HALF * r(2:ldx,1:m)

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

            r(1,2:mj) = HALF * r(1,2:mj)

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
                ierror = 20
                call fft_error_handler('rfft2f',-5)
                return
            end if

            call w2r(ldim,ldw,l,m,r,work)
        end if

        if (modl == 0) then

            call rfftmf(1, 1, m, ldim, r(l,1), m*ldim, wsave(lwsav+mwsav+1), &
                mmsav, work, size(work), local_error_flag)

            associate( mj => 2*((m+1)/2)-1 )

                r(l,2:mj) = HALF * r(l,2:mj)

            end associate

            r(l,3:m:2) = -r(l,3:m:2)

        end if

        if (local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('rfft2f',-5)
        end if

    end subroutine rfft2f


end submodule real_forward_2d
