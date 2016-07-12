submodule (auxiliary_routines) error_handlers

contains

    module procedure fft_consistent
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: lcm
        !--------------------------------------------------------------

        !
        !==> Set least common multiple of inc and jump.
        !
        lcm = get_least_common_multiple(inc, jump)

        !
        !==> Check consistency
        !
        if ( lcm <= (n - 1) * inc .and. lcm <= (lot - 1) * jump ) then
            return_value = .false.
        else
            return_value = .true.
        end if

    end procedure fft_consistent


    pure function get_least_common_multiple(a, b) result (return_value)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in) :: a
        integer (ip), intent (in) :: b
        integer (ip)              :: return_value
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: i, j, jnew !! Counters
        !--------------------------------------------------------------

        i = a
        j = b

        do while (j /= 0)
            jnew = mod(i, j)
            i = j
            j = jnew
        end do

        ! Return least common multiple of a and b
        return_value = (a * b) / i

    end function get_least_common_multiple



    module procedure fft_error_handler

        use, intrinsic :: iso_fortran_env, only: stderr => ERROR_UNIT

        write (stderr, '(/a)') ' fftpack routine:: '
        write (stderr, '(a)') ' FATAL ERROR '

        select case (info)
            case (-1)
                write (stderr, '(4(a))') '  On entry to ', &
                    trim(subroutine_name), &
                    ' parameters lot, jump, n and inc are inconsistent.'
            case (-2)
                write (stderr, '(4(a))')  '  On ntry to ', &
                    trim(subroutine_name), &
                    ' parameter l is greater than ldim.'
            case (-3)
                write (stderr, '(4(a))')  '  On entry to ', &
                    trim(subroutine_name), &
                    ' parameter m is greater than mdim.'
            case (-5)
                write (stderr, '(4(a))')  '  Within ', &
                    trim(subroutine_name), &
                    ' input error returned by lower level routine.'
            case (-6)
                write (stderr, '(4(a))')  '  On entry to ', &
                    trim(subroutine_name), &
                    ' parameter ldim is less than 2*(l/2+1).'
            case default
                write (stderr, '(3(a),i3,a)') '  On entry to ',&
                    trim(subroutine_name),&
                    ' parameter number ', info, ' had an illegal value.'
        end select

    end procedure fft_error_handler

end submodule error_handlers
