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
    complex (wp)            :: complex_data(n)
    real (wp)               :: real_part(n), imaginary_part(n)
    !------------------------------------------------------

    !
    !==> Identify test
    !
    write( stdout, '(A)') 'program tcfft1 and related messages:'

    !
    !==> Generate test vector for forward-backward test
    !
    call get_random_vector(real_part, imaginary_part, complex_data)

    !
    !==> Perform backward transform
    !
    call fft%ifft(complex_data)

    !
    !==> Perform forward transform
    !
    call fft%fft(complex_data)

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(complex_data-cmplx(real_part, imaginary_part, kind=wp))) )

        write( stdout, '(A,E23.15E3)' ) 'cfft1 backward-forward max error =', max_err

    end associate

    !
    !==> Generate test vector for forward-backward test
    !
    call get_random_vector(real_part, imaginary_part, complex_data)

    !
    !==> Perform forward transform
    !
    call fft%fft(complex_data)


    !
    !==> Perform backward transform
    !
    call fft%ifft(complex_data)


    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(complex_data-cmplx(real_part, imaginary_part, kind=wp))) )

        write( stdout, '(A,E23.15E3)' ) 'cfft1 forward-backward max error =', max_err
        write( stdout, '(A, /)') 'end program tcfft1 and related messages'

    end associate

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


