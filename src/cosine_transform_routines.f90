module cosine_transform_routines

    use auxiliary_routines, only: &
        get_1d_saved_workspace_length, &
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
    public :: cost1i, cost1b, cost1f, costmi, costmb, costmf


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
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
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
        !--------------------------------------------------------------

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



    subroutine cost1b(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! cost1b: real backward cosine transform, 1d.
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
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
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
        !--------------------------------------------------------------

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



    subroutine costmf(lot, jump, n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        !  costmf: real cosine transform, multiple vectors.
        !
        !  purpose:
        !
        !  costmf computes the one-dimensional fourier transform of multiple
        !  even sequences within a real array.  this transform is referred to
        !  as the forward transform or fourier analysis, transforming the
        !  sequences from physical to spectral space.
        !
        !  this transform is normalized since a call to costmf followed
        !  by a call to costmb (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  parameters:
        !
        !  integer lot, the number of sequences to be transformed
        !  within array r.
        !
        !  integer jump, the increment between the locations,
        !  in array r, of the first elements of two consecutive sequences to
        !  be transformed.
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n-1 is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the same sequence.
        !
        !  input/real r(lenr), array containing lot sequences,
        !  each having length n.  on input, the data to be transformed; on output,
        !  the transormed data.  r can have any number of dimensions, but the total
        !  number of locations must be at least lenr.
        !
        !  integer lenr, the dimension of the  r array.
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

end module cosine_transform_routines
