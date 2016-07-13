submodule (real_transform_routines) real_initialization_routines

contains

    module subroutine rfft1i(n, wsave, lensav, ierror)
        !
        !  rfft1i: initialization for rfft1b and rfft1f.
        !
        !  purpose:
        !
        !  rfft1i initializes array wsave for use in its companion routines
        !  rfft1b and rfft1f.  the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  separate wsave arrays are required for different
        !  values of n.
        !
        !  parameters:
        !
        !  integer n, the length of the sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  real wsave(lensav), containing the prime factors of
        !  n and also containing certain trigonometric values which will be used in
        !  routines rfft1b or rfft1f.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n + int(log(real(n))) + 4.
        !
        !  integer ierror, error_flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: wsave(lensav)
        integer (ip), intent (in)  :: lensav
        integer (ip), intent (out) :: ierror
        !--------------------------------------------------------------

        if (lensav < n + int(log(real(n, kind=wp) )/log(TWO)) + 4) then
            ierror = 2
            call fft_error_handler('rfft1i ', 3)
        else
            ierror = 0
        end if

        if (n /= 1) call rffti1(n, wsave(1), wsave(n+1))

    end subroutine rfft1i



    module subroutine rfft2i(l, m, wsave, lensav, ierror)
        !
        ! rfft2i: initialization for rfft2b and rfft2f.
        !
        !  purpose:
        !  rfft2i initializes real array wsave for use in its companion routines
        !  rfft2f and rfft2b for computing the two-dimensional fast fourier
        !  transform of real data.  prime factorizations of l and m, together with
        !  tabulations of the trigonometric functions, are computed and stored in
        !  array wsave.  rfft2i must be called prior to the first call to rfft2f
        !  or rfft2b.  separate wsave arrays are required for different values of
        !  l or m.
        !
        !
        !  integer l, the number of elements to be transformed
        !  in the first dimension.  the transform is most efficient when l is a
        !  product of small primes.
        !
        !  integer m, the number of elements to be transformed
        !  in the second dimension.  the transform is most efficient when m is a
        !  product of small primes.
        !
        !  integer lensav, the number of elements in the wsave
        !  array.  lensav must be at least l + m + int(log(real(l)))
        !  + int(log(real(m))) + 8.
        !
        !  real wsave(lensav), containing the prime factors
        !  of l and m, and also containing certain trigonometric values which
        !  will be used in routines rfft2b or rfft2f.
        !
        !  integer ier, error_flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  20, input error returned by lower level routine.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: l
        integer (ip), intent (in)  :: m
        real (wp),    intent (out) :: wsave(lensav)
        integer (ip), intent (in)  :: lensav
        integer (ip), intent (out) :: ierror
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip) :: local_error_flag, lwsav, mmsav, mwsav
        !--------------------------------------------------------------

        ! initialize error flag
        ierror = 0

        !
        !==> verify lensav
        !
        lwsav = l+int(log(real(l, kind=wp) )/log(TWO)) + 4
        mwsav = 2*m+int(log(real(m, kind=wp) )/log(TWO)) + 4
        mmsav = m+int(log(real(m, kind=wp) )/log(TWO)) + 4

        if (lensav < lwsav+mwsav+mmsav) then
            ierror = 2
            call fft_error_handler('rfft2i', 4)
            return
        end if

        call rfftmi(l, wsave(1), lwsav, local_error_flag)

        if (local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('rfft2i',-5)
            return
        end if

        call cfftmi(m, wsave(lwsav+1), mwsav, local_error_flag)

        if (local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('rfft2i',-5)
            return
        end if

        call rfftmi(m, wsave(lwsav+mwsav+1), mmsav, local_error_flag)

        if (local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('rfft2i',-5)
            return
        end if

    end subroutine rfft2i



    module subroutine rfftmi(n, wsave, lensav, ierror)
        !
        ! rfftmi: initialization for rfftmb and rfftmf.
        !
        !  purpose:
        !
        !  rfftmi initializes array wsave for use in its companion routines
        !  rfftmb and rfftmf.  the prime factorization of n together with a
        !  tabulation of the trigonometric functions are computed and stored
        !  in array wsave.  separate wsave arrays are required for different
        !  values of n.
        !
        !  input
        !
        !  integer n, the length of each sequence to be
        !  transformed.  the transform is most efficient when n is a product of
        !  small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least n + int(log(real(n))) + 4.
        !
        !  output
        !  real wsave(lensav), work array containing the prime
        !  factors of n and also containing certain trigonometric
        !  values which will be used in routines rfftmb or rfftmf.
        !
        !  integer ierror, error_flag.
        !  0, successful exit;
        !  2, input parameter lensav not big enough.
        !
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        integer (ip), intent (in)  :: lensav
        real (wp),    intent (out) :: wsave(lensav)
        integer (ip), intent (out) :: ierror
        !--------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (lensav < n + int(log(real(n, kind=wp) )/log(TWO)) + 4) then
            ierror = 2
            call fft_error_handler('rfftmi ', 3)
            return
        else
            ierror = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) call mrfti1(n, wsave(1), wsave(n+1))

    end subroutine rfftmi



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
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: fac(15)
        real (wp),    intent (out) :: wa(n)
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip)            :: i, ib, ido, ii, iip, ipm, is
        integer (ip)            :: j, k1, l1, l2, ld
        integer (ip)            :: nf, nfm1, nl, nq, nr, ntry
        integer (ip), parameter :: NTRYH(*) = [4, 2, 3, 5]
        real (wp),    parameter :: TWO_PI = TWO * acos(-ONE)
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
                ntry = NTRYH(j)
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
                    if (nl /= 1) cycle inner_loop
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
                    fi = ZERO
                    do ii=3,ido,2
                        i = i+2
                        fi = fi + ONE
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



    subroutine mrfti1(n, wa, fac)
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
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: wa(n)
        real (wp),    intent (out) :: fac(15)
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip)            :: i, ib, ido, ii, iip, ipm, is
        integer (ip)            :: j, k1, l1, l2, ld
        integer (ip)            :: nf, nfm1, nl, nq, nr, ntry
        integer (ip), parameter :: NTRYH(*) = [4, 2, 3, 5]
        real (wp),    parameter :: TWO_PI = TWO * acos(-ONE)
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
                ntry = NTRYH(j)
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
                fi = ZERO
                do ii=3,ido,2
                    i = i+2
                    fi = fi + ONE
                    arg = fi*argld
                    wa(i-1) = cos(arg)
                    wa(i) = sin(arg)
                end do
                is = is+ido

            end do
            l1 = l2
        end do

    end subroutine mrfti1



end submodule real_initialization_routines
