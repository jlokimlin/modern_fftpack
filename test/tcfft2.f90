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
!     parameter(lensav=2*(l+m) + int(log(real(l))/log(2.))
!    .                         + int(log(real(m))/log(2.)) + 8)
program tcfft2

    use, intrinsic :: ISO_Fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT, &
        stderr => ERROR_UNIT

    use type_FFTpack, only: FFTpack

    ! Explicit typing only
    implicit none

    !------------------------------------------------------------------
    ! Dictionary
    !------------------------------------------------------------------
    type (FFTpack)          :: my_fft
    integer (ip), parameter :: l=100, m=100
    complex (wp)            :: complex_data(l,m)
    real (wp)               :: real_part(l,m), imaginary_part(l,m)
    !------------------------------------------------------


    !
    !==> Identify test and initialize FFT
    !
    write( stdout, '(A)') 'program fft2 (complex) and related messages:'

    !
    !==> generate test matrix for forward-backward test
    !
    call get_random_matrix(real_part, imaginary_part, complex_data)

    !
    !==> perform forward transform
    !
    call my_fft%fft2(complex_data)

    !
    !==> perform backward transform
    !
    call my_fft%ifft2(complex_data)

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(complex_data-cmplx(real_part,imaginary_part,kind=wp))) )

        write( stdout, '(A,E23.15E3)' ) 'fft2 (complex) forward-backward max error =', max_err

    end associate

    !==> generate test matrix for backward-forward test

    call get_random_matrix(real_part, imaginary_part, complex_data)

    !==> perform backward transform

    call my_fft%ifft2(complex_data)

    !==> perform forward transform
    !
    call my_fft%fft2(complex_data)

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(complex_data-cmplx(real_part,imaginary_part,kind=wp))) )

        write( stdout, '(A,E23.15E3)' ) 'fft2 (complex) backward-forward max error =', max_err
        write( stdout, '(A,/)') 'end program tcfft2 and related messages'

    end associate

contains

    subroutine get_random_matrix(real_part, imaginary_part, complex_data)
        !--------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------
        real (wp),    intent (out) :: real_part(:,:)
        real (wp),    intent (out) :: imaginary_part(:,:)
        complex (wp), intent (out) :: complex_data(:,:)
        !--------------------------------------------------

        call random_seed()
        call random_number(real_part)
        call random_number(imaginary_part)

        complex_data = cmplx(real_part, imaginary_part, kind=wp)

    end subroutine get_random_matrix

end program tcfft2
