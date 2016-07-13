module quarter_cosine_transform_routines

    use auxiliary_routines

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use real_transform_routines

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: cosq1i, cosq1b, cosq1f, cosqmi, cosqmb, cosqmf


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
        ! cosq1b: real backward cosine quarter wave transform, 1d.
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
        ! cosq1f: real forward cosine quarter wave transform, 1d.
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



    subroutine cosqmb(lot, jump, n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! cosqmb: real backward cosine quarter wave, multiple vectors.
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



    subroutine cosqmf(lot, jump, n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! cosqmf: real forward cosine quarter wave, multiple vectors.
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



end module quarter_cosine_transform_routines
