!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2011 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                     FFTPACK  version 5.1                      *
!     *                                                               *
!     *                 A Fortran Package of Fast Fourier             *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *               Paul Swarztrauber and Dick Valent               *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
program tcfft1

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT, &
        stderr => ERROR_UNIT

    use type_FFTpack, only: &
        FFTpack

    ! Explicit typing only
    implicit none

    !------------------------------------------------------
    ! Dictionary
    !------------------------------------------------------
    type (FFTpack)          :: fft
    integer (ip), parameter :: n = 10**3
    integer (ip)            :: error_flag
    complex (wp)            :: c(n)
    real (wp)               :: rr(n), ri(n)
    !------------------------------------------------------

    !
    !==> Allocate memory
    !
    fft = FFTpack(n)

    associate( &
        wsave => fft%saved_workspace, &
        work => fft%workspace &
        )
        !
        !==> Identify test and initialize FFT
        !
        write( stdout, '(A)') 'program tcfft1 and related messages:'
        call fft%cfft1i(n, wsave, size(wsave), error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stderr, '(A,I3,A)') 'error ', error_flag, ' in routine cfft1i'
            stop
        end if

        !
        !==> Generate test vector for forward-backward test
        !
        call get_random_vector(rr, ri, c)

        !
        !==> Perform backward transform
        !
        call fft%cfft1b(n, 1, c, n, &
            wsave, size(wsave), work, size(work), error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stderr, '(A,I3,A)') 'error ', error_flag, ' in routine cfft1b !'
            stop
        end if

        !
        !==> Perform forward transform
        !
        call fft%cfft1f(n, 1, c, n, wsave, size(wsave), work, size(work), error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stderr, '(A,I3,A)') 'error ', error_flag, ' in routine cfft1f !'
            stop
        end if

        !
        !==> Print test results
        !
        associate( max_err => maxval(abs(c-cmplx(rr, ri, kind=wp))) )

            write( stdout, '(A,E23.15E3)' ) 'cfft1 backward-forward max error =', max_err

        end associate

        !
        !==> Identify test and initialize fft
        !
        write( stdout, '(A)') 'program tcfft1 and related messages:'

        call fft%cfft1i(n, wsave, size(wsave), error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stderr, '(A,I3,A)') 'error ', error_flag, ' in routine cfft1i'
            stop
        end if

        !
        !==> Generate test vector for forward-backward test
        !
        call get_random_vector(rr, ri, c)

        !
        !==> Perform forward transform
        !
        call fft%cfft1f(n, 1, c, n, wsave, size(wsave), work, size(work), error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stderr, '(A,I3,A)') 'error ', error_flag, ' in routine cfft1f !'
            stop
        end if

        !
        !==> Perform backward transform
        !
        call fft%cfft1b(n, 1, c, n, wsave, size(wsave), work, size(work), error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stderr, '(A,I3,A)') 'error ', error_flag, ' in routine cfft1b !'
            stop
        end if


        !
        !==> Print test results
        !
        associate( max_err => maxval(abs(c-cmplx(rr, ri, kind=wp))) )

            write( stdout, '(A,E23.15E3)' ) 'cfft1 forward-backward max error =', max_err
            write( stdout, '(A, /)') 'end program tcfft1 and related messages'

        end associate

    end associate
    !
    !==> Release memory
    !
    call fft%destroy()

contains

    subroutine get_random_vector(real_part, imaginary_part, complex_data)
        !--------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------
        real (wp),    intent (out) :: real_part(:)
        real (wp),    intent (out) :: imaginary_part(:)
        complex (wp), intent (out) :: complex_data(:)
        !--------------------------------------------------

        call random_seed()
        call random_number(real_part)
        call random_number(imaginary_part)

        complex_data = cmplx(real_part, imaginary_part, kind=wp)

    end subroutine get_random_vector

end program tcfft1


