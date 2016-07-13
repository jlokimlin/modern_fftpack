module quarter_sine_transform_routines

    use auxiliary_routines

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use real_transform_routines

    use quarter_cosine_transform_routines

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
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

    subroutine sinq1b(n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! sinq1b: real backward sine quarter wave transform, 1d.
        !
        !  purpose:
        !
        !  computes the one-dimensional fourier transform of a sequence
        !  which is a sine series with odd wave numbers.  this transform is
        !  referred to as the backward transform or fourier synthesis,
        !  transforming the sequence from spectral to physical space.
        !
        !  this transform is normalized since a call to sinq1b followed
        !  by a call to sinq1f (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr), on input, the sequence to be
        !  transformed.  on output, the transformed sequence.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to sinq1i before the first call to routine sinq1f
        !  or sinq1b for a given transform length n.  wsave's contents may be
        !  re-used for subsequent calls to sinq1f and sinq1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least n.
        !
        !  integer ier, the error_flag.
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
        ! sinq1f: real forward sine quarter wave transform, 1d.
        !
        !  purpose:
        !
        !  computes the one-dimensional fourier transform of a sequence
        !  which is a sine series of odd wave numbers.  this transform is
        !  referred to as the forward transform or fourier analysis, transforming
        !  the sequence from physical to spectral space.
        !
        !  this transform is normalized since a call to sinq1f followed
        !  by a call to sinq1b (or vice-versa) reproduces the original
        !  array within roundoff error.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the sequence.
        !
        !  input/real r(lenr), on input, the sequence to be
        !  transformed.  on output, the transformed sequence.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to sinq1i before the first call to routine sinq1f
        !  or sinq1b for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to sinq1f and sinq1b with the same n.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  workspace, real (wp) work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least n.
        !
        !  integer ier, the error_flag.
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
        !  purpose:
        !
        !  initializes array wsave for use in its companion routines
        !  sinq1f and sinq1b. the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave. separate wsave arrays are required for different
        !  values of n.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least 2*n + int(log(real(n))) + 4.
        !
        !  real wsave(lensav), containing the prime factors
        !  of n and also containing certain trigonometric values which will be used
        ! in routines sinq1b or sinq1f.
        !
        !  integer ier, error_flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
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



    subroutine sinqmb(lot, jump, n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! sinqmb: real backward sine quarter wave, multiple vectors.
        !
        !  purpose:
        !
        !  sinqmb computes the one-dimensional fourier transform of multiple
        !  sequences within a real array, where each of the sequences is a
        !  sine series with odd wave numbers.  this transform is referred to as
        !  the backward transform or fourier synthesis, transforming the
        !  sequences from spectral to physical space.
        !
        !  this transform is normalized since a call to sinqmb followed
        !  by a call to sinqmf (or vice-versa) reproduces the original
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
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations, in
        !  array r, of two consecutive elements within the same sequence.
        !
        !  input/real r(lenr), containing lot sequences, each
        !  having length n.  r can have any number of dimensions, but the total
        !  number of locations must be at least lenr.  on input, r contains the data
        !  to be transformed, and on output the transformed data.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least (lot-1)*jump + inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to sinqmi before the first call to routine sinqmf
        !  or sinqmb for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to sinqmf and sinqmb with the same n.
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



    subroutine sinqmf(lot, jump, n, inc, x, lenx, wsave, lensav, work, lenwrk, ier)
        !
        ! sinqmf: real forward sine quarter wave, multiple vectors.
        !
        !  purpose:
        !
        !  sinqmf computes the one-dimensional fourier transform of multiple
        !  sequences within a real array, where each sequence is a sine series
        !  with odd wave numbers.  this transform is referred to as the forward
        !  transform or fourier synthesis, transforming the sequences from
        !  spectral to physical space.
        !
        !  this transform is normalized since a call to sinqmf followed
        !  by a call to sinqmb (or vice-versa) reproduces the original
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
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer inc, the increment between the locations,
        !  in array r, of two consecutive elements within the same sequence.
        !
        !  input/real r(lenr), containing lot sequences, each
        !  having length n.  r can have any number of dimensions, but the total
        !  number of locations must be at least lenr.  on input, r contains the data
        !  to be transformed, and on output the transformed data.
        !
        !  integer lenr, the dimension of the r array.
        !  lenr must be at least (lot-1)*jump + inc*(n-1)+ 1.
        !
        !  input, real (wp) wsave(lensav).  wsave's contents must be
        !  initialized with a call to sinqmi before the first call to routine sinqmf
        !  or sinqmb for a given transform length n.  wsave's contents may be re-used
        !  for subsequent calls to sinqmf and sinqmb with the same n.
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
        !  sinqmi: initialization for sinqmb and sinqmf.
        !
        !  purpose:
        !
        !  sinqmi initializes array wsave for use in its companion routines
        !  sinqmf and sinqmb.  the prime factorization of n together with a
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
        !  real wsave(lensav), containing the prime factors
        !  of n and also containing certain trigonometric values which will be used
        !  in routines sinqmb or sinqmf.
        !
        !  integer ier, error_flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
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



end module quarter_sine_transform_routines
