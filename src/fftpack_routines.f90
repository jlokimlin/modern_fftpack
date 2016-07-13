module fftpack_routines

    use, intrinsic :: iso_fortran_env, only: &
        stderr => ERROR_UNIT

    use auxiliary_routines

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use complex_transform_routines

    use real_transform_routines

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: cfft1i, cfft1b, cfft1f, cfft2i, cfft2b, cfft2f, cfftmi, cfftmb, cfftmf
    public :: rfft1i, rfft1b, rfft1f, rfft2i, rfft2b, rfft2f, rfftmi, rfftmb, rfftmf
    public :: cost1i, cost1b, cost1f, costmi, costmb, costmf
    public :: sint1i, sint1b, sint1f, sintmi, sintmb, sintmf
    public :: cosq1i, cosq1b, cosq1f, cosqmi, cosqmb, cosqmf
    public :: sinq1i, sinq1b, sinq1f, sinqmi, sinqmb, sinqmf


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
        real (wp) x1

        !
        !==> Check validity of input arguments
        !
        if (lenx < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cosq1b', 6)
            return
        else if (lensav < get_cost_1d_saved_workspace_length(n) ) then
            ier = 2
            call fft_error_handler('cosq1b', 8)
            return
        else if (lenwrk < get_cost_1d_workspace_length(n)) then
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
        end if

        select case (n)
            case (2)
                x1 = x(1,1)+x(1,2)
                x(1,2) = (x(1,1)-x(1,2))/sqrt(TWO)
                x(1,1) = x1
            case default

                call cosqb1(n,inc,x,wsave,work,local_error_flag)

                ! check error_flag
                if (local_error_flag /= 0) then
                    ier = 20
                    call fft_error_handler('cosq1b',-5)
                end if

        end select

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
        else if (lensav < get_cost_1d_saved_workspace_length(n) ) then
            ier = 2
            call fft_error_handler('cosq1f', 8)
            return
        else if (lenwrk < get_cost_1d_workspace_length(n)) then
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
            tsqx = x(1,2)/sqrt(TWO)
            x(1,2) = HALF *x(1,1)-tsqx
            x(1,1) = HALF *x(1,1)+tsqx
        else
            ! Peform cosine transform
            call cosqf1(n,inc,x,wsave,work,local_error_flag)

            ! Check error_flag
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
        !  integer ier, error_flag.
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
        real (wp), parameter :: HALF_PI = acos(-ONE)/2
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lensav < get_cost_1d_saved_workspace_length(n)) then
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
        fk = ZERO

        do k=1,n
            fk = fk + ONE
            wsave(k) = cos(fk*dt)
        end do

        associate( lnsv => n+int(log(real(n, kind=wp))/log(TWO))+4 )

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
            x(1,i) = HALF * (x(1,i-1)-x(1,i))
            x(1,i-1) = HALF * xim1
        end do

        x(1,1) = HALF * x(1,1)
        modn = mod(n,2)

        if (modn == 0) then
            x(1,n) = HALF * x(1,n)
        end if

        associate( &
            lenx => inc*(n-1)  + 1, &
            lnsv => n + int(log(real(n, kind=wp) )/log(TWO)) + 4, &
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
            x(1,ns2+1) = TWO * wsave(ns2) * x(1,ns2+1)
        end if

        do k=2,ns2
            kc = np2-k
            x(1,k) = work(k)+work(kc)
            x(1,kc) = work(k)-work(kc)
        end do

        x(1,1) = TWO * x(1,1)

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
            work(ns2+1) = TWO * x(1,ns2+1)
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
            lnsv => n + int(log(real(n, kind=wp) )/log(TWO)) + 4, &
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
            xim1 = HALF * (x(1,i-1)+x(1,i))
            x(1,i) = HALF * (x(1,i-1)-x(1,i))
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
        else if (lensav < get_1d_saved_workspace_length(n)) then
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
                x(m,2) = (x(m,1)-x(m,2))/sqrt(TWO)
                x(m,1) = x1
            end do
        else
            call mcsqb1(lot,jump,n,inc,x,wsave,work,local_error_flag)

            ! Check error_flag
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
        else if (lensav < get_1d_saved_workspace_length(n)) then
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
                tsqx = x(m,2)/sqrt(TWO)
                x(m,2) = HALF * x(m,1)-tsqx
                x(m,1) = HALF * x(m,1)+tsqx
            end do
        else
            call mcsqf1(lot,jump,n,inc,x,wsave,work,local_error_flag)

            ! Check error_flag
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
        !  integer ier, error_flag.
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
        real (wp), parameter ::  HALF_PI = acos(-ONE)/2
        real (wp) wsave(lensav)

        if (lensav < get_1d_saved_workspace_length(n)) then
            ier = 2
            call fft_error_handler('cosqmi', 3)
            return
        else
            ier = 0
        end if

        dt = HALF_PI/n
        fk = ZERO

        do k=1,n
            fk = fk + ONE
            wsave(k) = cos(fk*dt)
        end do

        ! Set workspace index pointer
        lnsv = get_1d_saved_workspace_length(n) - n

        call rfftmi(n, wsave(n+1), lnsv, local_error_flag)

        ! Check error_flag
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
            call fft_error_handler('cost1b', 6)
            return
        else if (lensav < get_1d_saved_workspace_length(n)) then
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
            call fft_error_handler('cost1f', 6)
            return
        else if (lensav < get_1d_saved_workspace_length(n)) then
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
        !  integer ier, error_flag.
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
        real (wp), parameter :: PI = acos(-ONE)
        real (wp) wsave(lensav)


        !
        !==> Check validity of input arguments
        !
        if (lensav < get_1d_saved_workspace_length(n)) then
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
            fk = ZERO

            do k=2,ns2
                kc = np1-k
                fk = fk + ONE
                wsave(k) = TWO * sin(fk*dt)
                wsave(kc) = TWO * cos(fk*dt)
            end do

            lnsv = nm1 + int(log(real(nm1, kind=wp) )/log(TWO)) +4

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
                    x(1,ns2+1) = TWO * x(1,ns2+1)
                end if

                lenx = inc*(nm1-1) + 1
                lnsv = nm1 + int(log(real(nm1, kind=wp))/log(TWO), kind=ip) + 4
                lnwk = nm1

                call rfft1f(nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

                if (local_error_flag /= 0) then
                    ier = 20
                    call fft_error_handler('costb1',-5)
                    return
                end if

                fnm1s2 = real(nm1, kind=wp)/2
                dsum = HALF * dsum
                x(1,1) = fnm1s2*x(1,1)

                if (mod(nm1,2) == 0) then
                    x(1,nm1) = TWO * x(1,nm1)
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
            x(1,2) = HALF * (x(1,1)-x(1,2))
            x(1,1) = HALF * x1h
        else
            if (3 >= n) then
                x1p3 = x(1,1)+x(1,3)
                tx2 = x(1,2)+x(1,2)
                x(1,2) = HALF * (x(1,1)-x(1,3))
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

                ! Set workspace index pointers
                lenx = inc*(nm1-1)  + 1
                lnsv = nm1 + int(log(real(nm1, kind=wp) )/log(TWO)) + 4
                lnwk = nm1

                call rfft1f(nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

                if (local_error_flag /= 0) then
                    ier = 20
                    call fft_error_handler('costf1',-5)
                    return
                end if

                snm1 = ONE /nm1
                dsum = snm1*dsum

                if (mod(nm1,2) == 0) then
                    x(1,nm1) = x(1,nm1)+x(1,nm1)
                end if

                do i=3,n,2
                    xi = HALF * x(1,i)
                    x(1,i) = HALF * x(1,i-1)
                    x(1,i-1) = dsum
                    dsum = dsum+xi
                end do

                if (modn == 0) then
                    x(1,n) = dsum
                end if

                x(1,1) = HALF * x(1,1)
                x(1,n) = HALF * x(1,n)
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
        else if (lensav < get_1d_saved_workspace_length(n)) then
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

        ! Set workspace index pointer
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
        !  integer IER, error_flag.
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
        else if (lensav < get_1d_saved_workspace_length(n)) then
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

        ! Check error_flag
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
        !  integer ier, error_flag.
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
        real (wp), parameter :: PI = acos(-ONE)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lensav < get_1d_saved_workspace_length(n)) then
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
            fk = ZERO

            do k=2,ns2
                kc = np1-k
                fk = fk + ONE
                wsave(k) = TWO * sin(fk*dt)
                wsave(kc) = TWO * cos(fk*dt)
            end do

            lnsv = nm1 + int(log(real(nm1, kind=wp) )/log(TWO)) + 4

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
                x(m,i) = HALF * (x(m,i-1)-x(m,i))
                x(m,i-1) = HALF * xim1
            end do
        end do

        do m=1,lj,jump
            x(m,1) = HALF * x(m,1)
        end do

        modn = mod(n,2)
        if (modn == 0) then
            do m=1,lj,jump
                x(m,n) = HALF * x(m,n)
            end do
        end if

        lenx = (lot-1)*jump + inc*(n-1)  + 1
        lnsv = n + int(log(real(n, kind=wp) )/log(TWO)) + 4
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
        lnsv = n + int(log(real(n, kind=wp) )/log(TWO)) + 4
        lnwk = lot*n

        call rfftmf(lot,jump,n,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

        if (local_error_flag /= 0) then
            ier = 20
            call fft_error_handler('mcsqf1',-5)
            return
        end if

        do i=3,n,2
            do m=1,lj,jump
                xim1 = HALF * (x(m,i-1)+x(m,i))
                x(m,i) = HALF * (x(m,i-1)-x(m,i))
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
                lnsv = nm1 + int(log(real(nm1, kind=wp))/log(TWO)) + 4
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
                    dsum(m1) = HALF * dsum(m1)
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
                x(m,2) = HALF * (x(m,1)-x(m,2))
                x(m,1) = HALF * x1h
            end do
        else
            if (3 >= n) then
                do m=1,lj,jump
                    x1p3 = x(m,1)+x(m,3)
                    tx2 = x(m,2)+x(m,2)
                    x(m,2) = HALF * (x(m,1)-x(m,3))
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
                lnsv = nm1 + int(log(real(nm1, kind=wp))/log(TWO)) + 4
                lnwk = lot*nm1

                call rfftmf(lot,jump,nm1,inc,x,lenx,wsave(n+1),lnsv,work,lnwk,local_error_flag)

                if (local_error_flag /= 0) then
                    ier = 20
                    call fft_error_handler('mcstf1',-5)
                    return
                end if

                snm1 = ONE/nm1

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
                        xi = HALF * x(m,i)
                        x(m,i) = HALF * x(m,i-1)
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
                    x(m,1) = HALF * x(m,1)
                    x(m,n) = HALF * x(m,n)
                end do
            end if
        end if

    end subroutine mcstf1



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
        !  integer IER, the error_flag.
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
        else if (lensav < get_1d_saved_workspace_length(n)) then
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
        !  integer IER, the error_flag.
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
        else if (lensav < get_1d_saved_workspace_length(n)) then
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

            ! check error_flag
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
        !  integer IER, error_flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  20, input error returned by lower level routine.
        !

        integer (ip) lensav
        integer (ip) ier
        integer (ip) local_error_flag
        integer (ip) n
        real (wp) wsave(lensav)

        if (lensav < get_1d_saved_workspace_length(n)) then
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
        !  integer IER, error_flag.
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
        else if (lensav < get_1d_saved_workspace_length(n)) then
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
            x(1:lj:jump,1) =  FOUR * x(1:lj:jump,1)
        else
            ns2 = n/2
            x(2:n:2,1:lj:jump) = -x(1:lj:jump,2:n:2)

            call cosqmb(lot,jump,n,inc,x,lenx,wsave,lensav,work,lenwrk,local_error_flag)

            ! Check error_flag
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
        !  integer IER, error_flag.
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
        else if (lensav < get_1d_saved_workspace_length(n)) then
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
        !  integer IER, error_flag.
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

        if (lensav < get_1d_saved_workspace_length(n)) then
            ier = 2
            call fft_error_handler('sinqmi', 3)
            return
        end if

        call cosqmi(n, wsave, lensav, local_error_flag)

        ! Check error_flag
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
        !  integer IER, error_flag.
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
        !  integer IER, error_flag.
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
        !  integer IER, error_flag.
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
        !  integer IER, error_flag.
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

end module fftpack_routines
