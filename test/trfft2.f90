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
! parameter(lensav= l + 3*m + int(log(real(l))/log(2.))
!    .      + 2*int(log(real(m))/log(2.)) + 12)
!
program trfft2

    use, intrinsic :: iso_fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT, &
        stderr => ERROR_UNIT, &
        compiler_version, &
        compiler_options

    use type_FFTpack, only: &
        FFTpack

    ! Explicit typing only
    implicit none

    !------------------------------------------------------------------
    ! Dictionary
    !------------------------------------------------------------------
    type (FFTpack)          :: fft
    integer (ip)            :: error_flag
    integer (ip), parameter :: l=100, m=100, ldim=l
    integer (ip), parameter :: lensav= 430
    integer (ip), parameter :: lenwrk=(l+1)*m
    real (wp)               :: real_data(l,m), data_copy(l,m)
    real (wp)               :: wsave(lensav), work(lenwrk)
    !------------------------------------------------------------------

    !
    !==> Identify test and initialize FFT
    !
    write( stdout, '(A)') 'program trfft2 and related messages:'
    call fft%rfft2i(l,m,wsave, size(wsave),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine rfft2i'
        stop
    end if

    !
    !==> Generate test matrix for forward-backward test
    !
    call random_seed()
    call random_number(real_data)
    data_copy = real_data

    !
    !==> Perform forward transform
    !
    call fft%rfft2f (ldim,l,m,real_data, wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine rfft2f !'
        stop
    end if

    !
    !==> Perform backward transform
    !
    call fft%rfft2b (ldim,l,m,real_data, wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine rfft2b !'
        stop
    end if


    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(real_data-data_copy)) )

        write( stdout, '(A,E23.15E3)' ) 'rfft2 forward-backward max error =', max_err

    end associate

    !==> Generate test matrix for backward-forward test

    call random_seed()
    call random_number(real_data)
    data_copy = real_data

    !
    !==> Perform backward transform
    !
    call fft%rfft2b (ldim,l,m,real_data, &
        wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine rfft2b !'
        stop
    end if

    !
    !==> Perform forward transform
    !
    call fft%rfft2f (ldim,l,m,real_data, wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine rfft2f !'
        stop
    end if

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(real_data-data_copy)) )

        write( stdout, '(A,E23.15E3)' ) 'rfft2 backward-forward max error =', max_err
        write( stdout, '(A,/)') 'end program trfft2 and related messages'

    end associate

    !
    !==> Print compiler info
    !
    write( stdout, '(A)' ) ' '
    write( stdout, '(4A)' ) 'This result was compiled by ', &
        compiler_version(), ' using the options ', &
        compiler_options()
    write( stdout, '(A)' ) ' '

end program trfft2
