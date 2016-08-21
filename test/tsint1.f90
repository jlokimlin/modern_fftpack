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
! parameter(lensav= n/2 + n + int(log(real(n))/log(2.)) + 4)
!
program tsint1

    use, intrinsic :: ISO_Fortran_env, only: &
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
    type (FFTpack)          :: fft
    integer (ip)            :: error_flag
    integer (ip), parameter :: n = 10**3
    real (wp)               :: real_data(n), data_copy(n)
    real (wp), allocatable  :: wsave(:), work(:)
    !------------------------------------------------------

    !
    !==> Allocate memory
    !
    wsave = fft%get_1d_saved_workspace(n)
    work = fft%get_1d_sin_workspace(n)

    !
    !==> Identify program and initialize FFT
    !
    write( stdout, '(A)') 'program tsint1 and related messages:'
    call fft%sint1i(n,wsave, size(wsave),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine sint1i'
        stop
    end if

    !
    !==> Generate test vector for forward-backward test
    !
    call random_seed()
    call random_number(real_data)
    data_copy = real_data

    !
    !==> Perform forward transform
    !
    call fft%sint1f(n,1,real_data,n, wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine sint1f !'
        stop
    end if

    !
    !==> Perform backward transform
    !
    call fft%sint1b(n,1,real_data,n, &
        wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine sint1b !'
        stop
    end if

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(real_data-data_copy)) )

        write( stdout, '(A,E23.15E3)' ) 'sint1 forward-backward max error =', max_err

    end associate

    !
    !==> Generate test vector for backward-forward test
    !
    call random_seed()
    call random_number(real_data)
    data_copy = real_data

    !
    !==> Perform backward transform
    !
    call fft%sint1b(n,1,real_data,n, &
        wsave, size(wsave),work, size(work),error_flag)
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine sint1b !'
        stop
    end if

    !
    !==> Perform forward transform
    !
    call fft%sint1f(n,1,real_data,n, wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine sint1f !'
        stop
    end if

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(real_data-data_copy)) )

        write( stdout, '(A,E23.15E3)' ) 'sint1 backward-forward max error =', max_err
        write( stdout, '(A,/)') 'end program tsint1 and related messages'

    end associate

    !
    !==> Release memory
    !
    deallocate( wsave )
    deallocate( work )

end program tsint1
