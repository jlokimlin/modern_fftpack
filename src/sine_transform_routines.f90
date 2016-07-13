module sine_transform_routines

    use auxiliary_routines, only: &
        fft_error_handler, &
        fft_consistent

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use real_transform_routines, only: &
        rfft1i, &
        rfft1f, &
        rfftmi, &
        rfftmf

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: sint1i, sint1b, sint1f, sintmi, sintmb, sintmf


    !---------------------------------------------------------------------------------
    ! Variables confined to the module
    !---------------------------------------------------------------------------------
    real (wp), parameter :: ZERO = 0.0_wp
    real (wp), parameter :: HALF = 0.5_wp
    real (wp), parameter :: ONE = 1.0_wp
    real (wp), parameter :: TWO = 2.0_wp
    real (wp), parameter :: THREE = 3.0_wp
    real (wp), parameter :: FOUR = 4.0_wp
    real (wp), parameter :: FIVE = 5.0_wp
    !---------------------------------------------------------------------------------


contains


    subroutine sint1b(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        !  sint1b: real backward sine transform, 1d.
        !
        !  purpose:
        !
        !  sint1b computes the one-dimensional fourier transform of an odd
        !  sequence within a real array.  this transform is referred to as
        !  the backward transform or fourier synthesis, transforming the
        !  sequence from spectral to physical space.
        !
        !  this transform is normalized since a call to sint1b followed
        !  by a call to sint1f (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n+1 is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr), on input, contains the sequence
        !  to be transformed, and on output, the transformed sequence.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to sint1i before the first call to routine sint1f
        !  or sint1b for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to sint1f and sint1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n/2 + n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*n+2.
        !
        !  integer ier, error_flag.
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
            call fft_error_handler('sint1b', 6)
            return
        else if ( lensav < n / 2 + n + int(log(real(n, kind=wp) ) &
            / log(TWO ) ) + 4 ) then
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

        ! Check error_flag
        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('sint1b',-5)
        end if

    end subroutine sint1b


    subroutine sint1f(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        !  sint1f: real forward sine transform, 1d.
        !
        !  purpose:
        !
        !  sint1f computes the one-dimensional fourier transform of an odd
        !  sequence within a real array.  this transform is referred to as the
        !  forward transform or fourier analysis, transforming the sequence
        !  from physical to spectral space.
        !
        !  this transform is normalized since a call to sint1f followed
        !  by a call to sint1b (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n+1 is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr), on input, contains the sequence
        !  to be transformed, and on output, the transformed sequence.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to sint1i before the first call to routine sint1f
        !  or sint1b for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to sint1f and sint1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n/2 + n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*n+2.
        !
        !  integer ier, error_flag.
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
            call fft_error_handler('sint1f', 6)
            return
        else if (lensav < n/2 + n + int(log(real(n, kind=wp) ) &
            /log(TWO)) +4) then
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
        !  sint1i: initialization for sint1b and sint1f.
        !
        !  purpose:
        !
        !  sint1i initializes array wsave for use in its companion routines
        !  sint1f and sint1b.  the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  separate wsave arrays are required for different
        !  values of n.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n+1 is a product
        !  of small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n/2 + n + int(log(real(n))) + 4.
        !
        !  real wsave(lensav), containing the prime factors
        !  of n and also containing certain trigonometric values which will be used
        !  in routines sint1b or sint1f.
        !
        !  integer ier, error_flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
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
        real (wp), parameter :: pi = acos(-ONE)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lensav < n/2 + n + int(log(real(n, kind=wp) ) &
            /log(TWO)) +4) then
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
                wsave(k) = TWO *sin(k*dt)
            end do

            lnsv = np1 + int(log(real(np1, kind=wp))/log(TWO)) +4

            call rfft1i(np1, wsave(ns2+1), lnsv, local_error_flag)

            ! Check error_flag
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
        real (wp), parameter :: HALF_SQRT3 = sqrt(THREE)/2
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
                xh(ns2+2) =  FOUR * x(1,ns2+1)
            end if

            xh(1) = ZERO
            lnxh = np1
            lnsv = np1 + int(log(real(np1, kind=wp))/log(TWO)) + 4
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
            xhold = (x(1,1)+x(1,2))/sqrt(THREE)
            x(1,2) = (x(1,1)-x(1,2))/sqrt(THREE)
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
                xh(ns2+2) =  FOUR * x(1,ns2+1)
            end if

            xh(1) = ZERO
            lnxh = np1
            lnsv = np1 + int(log(real(np1, kind=wp))/log(TWO)) + 4
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

            sfnp1 = ONE/np1
            x(1,1) = HALF * xh(1)
            dsum = x(1,1)

            do i=3,n,2
                x(1,i-1) = HALF * xh(i)
                dsum = dsum + HALF * xh(i-1)
                x(1,i) = dsum
            end do

            if (modn == 0) then
                x(1,n) = HALF * xh(n+1)
            end if
        end if

    end subroutine sintf1


    subroutine sintmb(lot, jump, n, inc, x, lenx, wsave, lensav, &
        work, lenwrk, ier)
        !
        ! sintmb: real backward sine transform, multiple vectors.
        !
        !  purpose:
        !
        !  sintmb computes the one-dimensional fourier transform of multiple
        !  odd sequences within a real array.  this transform is referred to as
        !  the backward transform or fourier synthesis, transforming the
        !  sequences from spectral to physical space.
        !
        !  this transform is normalized since a call to sintmb followed
        !  by a call to sintmf (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  parameters:
        !
        !  integer lot, the number of sequences to be transformed
        !  within the array r.
        !
        !  integer jump, the increment between the locations, in
        !  array r, of the first elements of two consecutive sequences.
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
        !  integer ier, error_flag.
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
            /log(TWO)) +4) then
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
        !  integer ier, error_flag.
        !  0, successful exit;
        !  1, input parameter lenr not big enough;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  4, input parameters inc,jump,n,lot are not consistent;
        !  20, input error returned by lower level routine.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
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
            / log(TWO ) ) + 4 ) then
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

        ! Check error_flag
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
        !  integer ier, error_flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  20, input error returned by lower level routine.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
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
        real (wp), parameter :: PI = acos(-ONE)
        real (wp) wsave(lensav)
        !--------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (lensav < n / 2 + n + int(log(real(n, kind=wp)) &
            / log(TWO ) ) + 4 ) then
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
                wsave(k) = TWO * sin(real(k, kind=wp) * dt)
            end do

            associate( nsv => np1+int(log(real(np1, kind=wp) )/log(TWO))+4)

                call rfftmi(np1, wsave(ns2+1), lnsv, local_error_flag)

            end associate

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('sintmi', -5)
            end if
        end if

    end subroutine sintmi


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
        real (wp), parameter :: HALF_SQRT3 = sqrt(THREE)/2
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
                    xh(m1,ns2+2) =  FOUR * x(m,ns2+1)
                end do
            end if

            xh(:,1) = ZERO

            lnxh = lot-1 + lot*(np1-1) + 1
            lnsv = np1 + int(log(real(np1, kind=wp))/log(TWO)) + 4
            lnwk = lot*np1

            call rfftmf(lot,1,np1,lot,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,local_error_flag)

            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('msntb1',-5)
                return
            end if

            if (mod(np1,2) == 0) then
                xh(:,np1) = TWO * xh(:,np1)
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
                associate( sqrt3 => sqrt(THREE) )
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
                    xh(m1,ns2+2) =  FOUR * x(m,ns2+1)
                end do
            end if

            xh(:,1) = ZERO

            lnxh = lot-1 + lot*(np1-1) + 1
            lnsv = np1 + int(log(real(np1, kind=wp))/log(TWO)) + 4
            lnwk = lot*np1

            call rfftmf(lot,1,np1,lot,xh,lnxh,wsave(ns2+1),lnsv,work,lnwk,local_error_flag)

            ! Check error_flag
            if (local_error_flag /= 0) then
                ier = 20
                call fft_error_handler('msntf1',-5)
                return
            end if

            if (mod(np1,2) == 0) then
                xh(:,np1) = TWO * xh(:,np1)
            end if


            sfnp1 = ONE/np1
            m1 = 0

            do m=1,lj,jump
                m1 = m1+1
                x(m,1) = HALF * xh(m1,1)
                dsum(m1) = x(m,1)
            end do

            do i=3,n,2
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    x(m,i-1) = HALF * xh(m1,i)
                    dsum(m1) = dsum(m1)+ HALF * xh(m1,i-1)
                    x(m,i) = dsum(m1)
                end do
            end do

            if (modn == 0) then
                m1 = 0
                do m=1,lj,jump
                    m1 = m1+1
                    x(m,n) = HALF * xh(m1,n+1)
                end do
            end if
        end if

    end subroutine msntf1


end module sine_transform_routines
