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
program tcosq1

    use, intrinsic :: ISO_Fortran_env, only: &
        wp => REAL64, &
        ip => INT32, &
        stdout => OUTPUT_UNIT, &
        stderr => ERROR_UNIT

    use type_FFTpack, only: FFTpack

    ! Explicit typing only
    implicit none

    !------------------------------------------------------
    ! Dictionary
    !------------------------------------------------------
    type (FFTpack)          :: fft
    integer (ip)            :: error_flag
    integer (ip), parameter :: n=10**3
    real (wp)               :: real_data(n), data_copy(n)
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
        write( stdout, '(A)') 'program tcosq1 and related messages:'
        call fft%cosq1i(n,wsave,size(wsave),error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine cosq1i'
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
        call fft%cosq1f(n,1,real_data,n, wsave,size(wsave),work, size(work),error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine cosq1f !'
            stop
        end if

        !==> Perform backward transform

        call fft%cosq1b(n,1,real_data,n, wsave,size(wsave),work, size(work),error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine cosq1b !'
            stop
        end if

        !
        !==> Print test results
        !
        associate( max_err => maxval(abs(real_data-data_copy)) )

            write( stdout, '(A,E23.15E3)' ) 'cosq1 forward-backward max error =', max_err

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
        call fft%cosq1b(n,1,real_data,n, wsave,size(wsave),work, size(work),error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine cosq1b !'
            stop
        end if

        !
        !==> Perform forward transform
        !
        call fft%cosq1f(n,1,real_data,n, wsave,size(wsave),work, size(work),error_flag)

        ! Check error flag
        if (error_flag /= 0) then
            write( stdout, '(A,I3,A)') 'error ',error_flag,' in routine cosq1f !'
            stop
        end if

        !
        !==> Print test results
        !
        associate( max_err => maxval(abs(real_data-data_copy)) )

            write( stdout, '(A,E23.15E3)' ) 'cosq1 backward-forward max error =', max_err
            write( stdout, '(A,/)') 'end program tcosq1 and related messages'

        end associate

    end associate
    !
    !==> Release memory
    !
    call fft%destroy()

end program tcosq1
