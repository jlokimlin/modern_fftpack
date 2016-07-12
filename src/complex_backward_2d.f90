submodule (complex_transform_routines:multiple_complex_backward) complex_backward_2d

contains

    module subroutine cfft2b(ldim, l, m, c, wsave, lensav, work, lenwrk, ierror)
        !
        ! cfft2b: complex fast ror transform, 2d.
        !
        !  purpose:
        !
        !  cfft2b computes the two-dimensional discrete ror transform of a
        !  complex periodic array.  this transform is known as the backward
        !  transform or synthesis, transforming from spectral to
        !  physical space.  routine cfft2b is normalized, in that a call to
        !  cfft2b followed by a call to cfft2f (or vice-versa) reproduces the
        !  original array within roundoff error.
        !
        !  bug fix
        !  on 10 may 2010, this code was modified by changing the value
        !  of an index into the wsave array.
        !
        !  parameters:
        !
        !  input
        !  integer ldim, the first dimension of c.
        !
        !  integer l, the number of elements to be transformed
        !  in the first dimension of the two-dimensional complex array c.  the value
        !  of l must be less than or equal to that of ldim.  the transform is
        !  most efficient when l is a product of small primes.
        !
        !  integer m, the number of elements to be transformed in
        !  the second dimension of the two-dimensional complex array c.  the transform
        !  is most efficient when m is a product of small primes.
        !
        !  real wsave(lensav). wsave's contents must be initialized with
        !  a call to cfft2i before the first call to routine cfft2f
        !  or cfft2b with transform lengths l and m.  wsave's contents may be
        !  re-used for subsequent calls to cfft2f and cfft2b with the same
        !  transform lengths l and m.
        !
        !  integer lensav, the dimension of the wsave array.
        !  lensav must be at least
        !
        !  2*(l+m) + int(log(real(l))) + int(log(real(m))) + 8.
        !
        !  real work(lenwrk).
        !
        !  integer lenwrk, the dimension of the work array.
        !  lenwrk must be at least 2*l*m.
        !
        !  input/output
        !  complex c(ldim,m), on intput, the array of two dimensions
        !  containing the (l,m) subarray to be transformed.  on
        !  output, the transformed data.
        !
        !
        !  output
        !
        !  integer ierror, the error_flag.
        !
        !  0, successful exit;
        !  2, input parameter lensav not big enough;
        !  3, input parameter lenwrk not big enough;
        !  5, input parameter ldim < l;
        !  20, input error returned by lower level routine.
        !
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: ldim
        integer (ip), intent (in)     :: l
        integer (ip), intent (in)     :: m
        complex (wp), intent (in out) :: c(ldim,m)
        real (wp),    intent (in)     :: wsave(lensav)
        integer (ip), intent (in)     :: lensav
        real (wp),    intent (in out) :: work(lenwrk)
        integer (ip), intent (in)     :: lenwrk
        integer (ip), intent (out)    :: ierror
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip)           :: local_error_flag
        real (wp), allocatable :: real_copy(:,:,:)
        !------------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if ( ldim < l ) then
            ierror = 5
            call fft_error_handler('cfft2b', -2)
            return
        else if (size(wsave) < get_complex_2d_saved_workspace_length(l, m)) then
            ierror = 2
            call fft_error_handler('cfft2b', 6)
            return
        else if (size(work) < get_complex_2d_workspace_length(l, m)) then
            ierror = 3
            call fft_error_handler('cfft2b', 8)
            return
        else
            ierror = 0
        end if


        !
        !==> Allocate memory
        !
        allocate( real_copy(2,ldim,m) )

        ! Make a copy: complex to real
        real_copy(1,:,:) = real(c)
        real_copy(2,:,:) = aimag(c)

        !
        !==> transform x lines of real_copy
        !
        associate( &
            iw =>  get_1d_saved_workspace_length(m)-1, &
            iw1 => (l-1) + ldim*(m-1)+1, &
            iw2 => get_1d_saved_workspace_length(m), &
            iw3 => get_2d_workspace_length(l, m), &
            local_error_flag => local_error_flag &
            )

            ! Perform transform
            call cfftmb(l, 1, m, ldim, real_copy, iw1 , wsave(iw), iw2, work, iw3, local_error_flag)

        end associate

        ! Check error_flag
        if (local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('cfft2b',-5)
            return
        end if
        !
        !==>  transform y lines of real_copy
        !
        associate( &
            iw => 1, &
            iw1 => (m-1)*ldim + l, &
            iw2 => get_1d_saved_workspace_length(l), &
            iw3 => get_2d_workspace_length(l, m), &
            local_error_flag => local_error_flag &
            )

            ! Perform transform
            call cfftmb(m, ldim, l, 1, real_copy, iw1, wsave(iw), iw2, work, iw3, local_error_flag)

        end associate

        ! Check error_flag
        if (local_error_flag /= 0) then
            ierror = 20
            call fft_error_handler('cfft2b',-5)
        end if

        ! Make copy: real to complex
        c =  cmplx(real_copy(1,:,:), real_copy(2,:,:), kind=wp)

        !
        !==> Release memory
        !
        deallocate( real_copy )

    end subroutine cfft2b

end submodule complex_backward_2d
