module fftpack_precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: wp, ip
    public :: pimach, epmach

    !-----------------------------------------------
    ! Dictionary: precision constants
    !-----------------------------------------------
    integer,   parameter :: ip = selected_int_kind(r=9)
    integer,   parameter :: sp = selected_real_kind(p=6, r=37)
    integer,   parameter :: dp = selected_real_kind(p=15, r=307)
    integer,   parameter :: qp = selected_real_kind(p=33, r=4931)
    integer,   parameter :: wp = dp
    !-----------------------------------------------

contains


    pure function pimach() result (return_value)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        real (wp) :: return_value
        !-----------------------------------------------

        return_value = 3.141592653589793238462643383279502884197169399375105820974_wp

    end function pimach


    pure function epmach() result (return_value)
        !
        ! Purpose:
        !
        ! Computes an approximate machine epsilon (accuracy), i.e.,
        !
        ! the smallest number x of the kind wp such that 1 + x > 1.
        !
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        real (wp) :: return_value
        !-----------------------------------------------

        return_value = epsilon(1.0_wp)

    end function epmach


end module fftpack_precision
