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

    use, intrinsic :: iso_fortran_env, only: &
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
    type (FFTpack)          :: fft
    integer (ip)            :: error_flag
    integer (ip), parameter :: l=100, m=100, ldim=l
    complex (wp)            :: c(l,m)
    real (wp)               :: rr(l,m), ri(l,m)
    real (wp), allocatable  :: wsave(:), work(:)
    !------------------------------------------------------

    !
    !==> Allocate memory
    !
    wsave = fft%get_2d_saved_workspace(l,m)
    work = fft%get_2d_workspace(l,m)

    !
    !==> Identify test and initialize FFT
    !
    write( stdout, '(A)') 'program tcfft2 and related messages:'
    call fft%cfft2i(l,m,wsave, size(wsave),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stderr, '(A,I3,A)') 'error ',error_flag,' in routine cfft2i'
        stop
    end if

    !
    !==> generate test matrix for forward-backward test
    !
    call get_random_matrix(rr, ri, c)

    !
    !==> perform forward transform
    !
    call fft%cfft2f(ldim,l,m,c, &
        wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stderr, '(A,I3,A)') 'error ',error_flag,' in routine cfft2f !'
        stop
    end if

    !
    !==> perform backward transform
    !
    call fft%cfft2b(ldim,l,m,c, wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stderr, '(A,I3,A)') 'error ',error_flag,' in routine cfft2b !'
        stop
    end if

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(c-cmplx(rr,ri,kind=wp))) )

        write( stdout, '(A,E23.15E3)' ) 'cfft2 forward-backward max error =', max_err

    end associate

    !==> generate test matrix for backward-forward test

    call get_random_matrix(rr, ri, c)

    !==> perform backward transform

    call fft%cfft2b(ldim,l,m,c, wsave, size(wsave),work, size(work),error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stderr, '(A,I3,A)') 'error ',error_flag,' in routine cfft2b !'
        stop
    end if

    !==> perform forward transform
    !
    call fft%cfft2f(ldim,l,m,c, wsave, size(wsave), work, size(work), error_flag)

    ! Check error flag
    if (error_flag /= 0) then
        write( stderr, '(A,I3,A)') 'error ',error_flag,' in routine cfft2f !'
        stop
    end if

    !
    !==> Print test results
    !
    associate( max_err => maxval(abs(c-cmplx(rr,ri,kind=wp))) )

        write( stdout, '(A,E23.15E3)' ) 'cfft2 backward-forward max error =', max_err
        write( stdout, '(A,/)') 'end program tcfft2 and related messages'

    end associate

    !
    !==> Release memory
    !
    deallocate( wsave )
    deallocate( work )

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
