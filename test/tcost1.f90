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
program tcost1

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT, &
        stderr => ERROR_UNIT

    use type_FFTpack, only: &
        FFTpack

    ! Explicit typing only
    implicit none

    !------------------------------------------------------------------
    ! Dictionary
    !------------------------------------------------------------------
    type (FFTpack)          :: my_fft
    integer (ip), parameter :: n=10**3
    real (wp)               :: real_data(n), data_copy(n)
    !------------------------------------------------------

    !
    !==> Identify test
    !
    write( stdout, '(A)') 'program tcost1 and related messages:'

    !
    !==> Generate test vector for forward-backward test
    !
    call get_random_vector(real_data, data_copy)

    !
    !==> Perform forward transform
    !
    call my_fft%dct(real_data)

    !
    !==> Perform backward transform
    !
    call my_fft%idct(real_data)

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(real_data-data_copy)) )

        write( stdout, '(A,E23.15E3)' ) 'cost1 forward-backward max error =', max_err

    end associate

    !
    !==> Generate test vector for backward-forward test
    !
    call get_random_vector(real_data, data_copy)

    !
    !==> Perform backward transform
    !
    call my_fft%idct(real_data)

    !
    !==> Perform forward transform
    !
    call my_fft%dct(real_data)
    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(real_data-data_copy)) )

        write( stdout, '(A,E23.15E3)' ) 'cost1 backward-forward max error =', max_err
        write( stdout, '(A,/)') 'end program tcost1 and related messages'

    end associate


contains


    subroutine get_random_vector(real_data, real_copy)
        !--------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------
        real (wp), intent (out) :: real_data(:)
        real (wp), intent (out) :: real_copy(:)
        !--------------------------------------------------

        call random_seed()
        call random_number(real_data)
        real_copy = real_data

    end subroutine get_random_vector

end program tcost1
