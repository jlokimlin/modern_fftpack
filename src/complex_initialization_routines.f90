submodule (complex_transform_routines) complex_initialization_routines

contains

    module subroutine cfft1i(n, wsave, lensav, ierror)
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
        !  ierror, error_flag.
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
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: iw1, iw2
        !--------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if ( size(wsave) < get_complex_1d_saved_workspace_length(n) ) then
            ierror = 2
            call fft_error_handler('cfftmi ', 3)
        else
            ierror = 0
        end if

        !
        !==> Perform transform
        !
        if (n == 1) return

        ! Set workspace indices
        iw1 = 2*n+1
        iw2 = iw1 + 1

        call mcfti1(n, wsave, wsave(iw1), wsave(iw2))

    end subroutine cfft1i


    module subroutine cfft2i(l, m, wsave, lensav, ierror)
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
        !  input
        !  integer l, the number of elements to be transformed
        !  in the first dimension. the transform is most efficient when l is a
        !  product of small primes.
        !
        !  integer m, the number of elements to be transformed
        !  in the second dimension. the transform is most efficient when m is a
        !  product of small primes.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least
        !
        !  2*(l+m) + int(log(real(l)))+ int(log(real(m))) + 8.
        !
        !  output
        !  real wsave(lensav), contains the prime factors of l
        !  and m, and also certain trigonometric values which will be used in
        !  routines cfft2b or cfft2f.
        !
        !  integer  ierror, error_flag.
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
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) local_error_flag
        !--------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (size(wsave) < get_complex_2d_saved_workspace_length(l, m)) then
            ierror = 2
            call fft_error_handler('cfft2i', 4)
            return
        else
            ierror = 0
        end if

        associate( lnsv => get_complex_1d_saved_workspace_length(l) )

            call cfftmi(l, wsave(1:), lnsv, local_error_flag)

        end associate

        ! Check error_flag
        if ( local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('cfft2i',-5)
            return
        end if

        associate( &
            iw => get_1d_saved_workspace_length(l) - 1, &
            lnsv => get_1d_saved_workspace_length(m) &
            )

            call cfftmi(m, wsave(iw:), lnsv, local_error_flag)

        end associate

        ! Check error_flag
        if (local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('cfft2i',-5)
        end if

    end subroutine cfft2i



    module subroutine cfftmi(n, wsave, lensav, ierror)
        !
        !  cfftmi: initialization for cfftmb and cfftmf.
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
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: iw1, iw2
        !--------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (lensav < get_complex_nd_saved_workspace_length(n) ) then
            ierror = 2
            call fft_error_handler('cfftmi ', 3)
            return
        else
            ierror = 0
        end if

        !
        !==>  Perform transform
        !
        if (n == 1) return

        ! Set workspace index pointer
        iw1 = 2*n+1
        iw2 = iw1 + 1

        call mcfti1(n, wsave, wsave(iw1), wsave(iw2:))


    end subroutine cfftmi



    pure subroutine mcfti1(n, wa, fnf, fac)
        !
        ! Purpose:
        !
        ! Sets up factors and tables
        !
        !--------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------
        integer (ip), intent (in)  :: n
        real (wp),    intent (out) :: wa(*)
        real (wp),    intent (out) :: fnf
        real (wp),    intent (out) :: fac(*)
        !--------------------------------------------------
        ! Dummy arguments
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

    end subroutine mcfti1



    pure subroutine get_factorization(n, nf, fac)
        !
        ! Purpose:
        !
        ! Factors of an integer for wp-kind float precision computations.
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
        ! Dummy arguments
        !--------------------------------------------------
        integer (ip), intent (in)  :: n
        integer (ip), intent (out) :: nf
        real (wp),    intent (out) :: fac(*)
        !--------------------------------------------------
        ! Local variables
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
        ! Computes trigonometric tables for wp-kind float precision arithmetic.
        !
        !--------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------
        integer (ip), intent (in)  :: ido
        integer (ip), intent (in)  :: iip
        real (wp),    intent (out) :: wa(ido,iip-1,2)
        !--------------------------------------------------
        ! Local variables
        !--------------------------------------------------
        integer (ip)         :: i, j !! Counters
        real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
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



end submodule complex_initialization_routines
