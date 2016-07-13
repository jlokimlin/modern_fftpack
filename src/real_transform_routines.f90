module real_transform_routines

    use auxiliary_routines

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use complex_transform_routines, only: &
        cfftmi, cfftmf, cfftmb

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: rfft1i, rfft1b, rfft1f, rfft2i, rfft2b, rfft2f, rfftmi, rfftmb, rfftmf

    interface
        !
        !==> 1D real initialization
        !
        module subroutine rfft1i(n, wsave, lensav, ierror)
            !
            !  rfft1i: initialization for rfft1b and rfft1f.
            !
            !  purpose:
            !
            !  rfft1i initializes array wsave for use in its companion routines
            !  rfft1b and rfft1f.  the prime factorization of n together with a
            !  tabulation of the trigonometric functions are computed and stored
            !  in array wsave.  separate wsave arrays are required for different
            !  values of n.
            !
            !  parameters:
            !
            !  integer n, the length of the sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  real wsave(lensav), containing the prime factors of
            !  n and also containing certain trigonometric values which will be used in
            !  routines rfft1b or rfft1f.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least n + int(log(real(n))) + 4.
            !
            !  integer ierror, error_flag.
            !  0, successful exit;
            !  2, input parameter lensav not big enough.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)  :: n
            real (wp),    intent (out) :: wsave(lensav)
            integer (ip), intent (in)  :: lensav
            integer (ip), intent (out) :: ierror
            !--------------------------------------------------------------
        end subroutine rfft1i
        !
        !==> 1D real backward
        !
        module subroutine rfft1b(n, inc, r, lenr, wsave, lensav, work, lenwrk, ierror)
            !
            ! purpose:
            !
            !  computes the one-dimensional fourier transform of a periodic
            !  sequence within a real array. this is referred to as the backward
            !  transform or fourier synthesis, transforming the sequence from
            !  spectral to physical space.  this transform is normalized since a
            !  call to rfft1b followed by a call to rfft1f (or vice-versa) reproduces
            !  the original array within roundoff error.
            !
            !  parameters:
            !
            !  integer n, the length of the sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer inc, the increment between the locations,
            !  in array r, of two consecutive elements within the sequence.
            !
            !  input/real r(lenr), on input, the data to be
            !  transformed, and on output, the transformed data.
            !
            !  integer lenr, the dimension of the r array.
            !  lenr must be at least inc*(n-1) + 1.
            !
            !  input, real (wp) wsave(lensav).  wsave's contents must be
            !  initialized with a call to rfft1i before the first call to routine
            !  rfft1f or rfft1b for a given transform length n.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least n + int(log(real(n))) + 4.
            !
            !  workspace, real (wp) work(lenwrk).
            !
            !  integer lenwrk, the dimension of the work array.
            !  lenwrk must be at least n.
            !
            !  integer ierror, error_flag.
            !  0, successful exit;
            !  1, input parameter lenr not big enough;
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip) lenr
            integer (ip) lensav
            integer (ip) lenwrk
            integer (ip) ierror
            integer (ip) inc
            integer (ip) n
            real (wp) r(lenr)
            real (wp) work(lenwrk)
            real (wp) wsave(lensav)
            !--------------------------------------------------------------

        end subroutine rfft1b
        !
        !==> 1D real forward
        !
        module subroutine rfft1f(n, inc, r, lenr, wsave, lensav, work, lenwrk, ierror)
            !
            ! rfft1f: 1d forward fast fourier transform
            !
            !  Purpose:
            !
            !  rfft1f computes the one-dimensional fourier transform of a periodic
            !  sequence within a real array. This is referred to as the forward
            !  transform or fourier analysis, transforming the sequence from physical
            !  to spectral space. This transform is normalized since a call to
            !  rfft1f followed by a call to rfft1b (or vice-versa) reproduces the
            !  original array within roundoff error.
            !
            !  Parameters:
            !
            !  integer n, the length of the sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer inc, the increment between the locations, in
            !  array r, of two consecutive elements within the sequence.
            !
            !  input/real r(lenr), on input, contains the sequence
            !  to be transformed, and on output, the transformed data.
            !
            !  integer lenr, the dimension of the r array.
            !  lenr must be at least inc*(n-1) + 1.
            !
            !  input, real (wp) wsave(lensav).  wsave's contents must be
            !  initialized with a call to rfft1i before the first call to routine rfft1f
            !  or rfft1b for a given transform length n.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least n + int(log(real(n))) + 4.
            !
            !  workspace, real (wp) work(lenwrk).
            !
            !  integer lenwrk, the dimension of the work array.
            !  lenwrk must be at least n.
            !
            !  integer ierror, error_flag.
            !  0, successful exit;
            !  1, input parameter lenr not big enough:
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: inc
            real (wp),    intent (in out) :: r(lenr)
            integer (ip), intent (in)     :: lenr
            real (wp),    intent (out)    :: work(lenwrk)
            integer (ip), intent (in)     :: lenwrk
            real (wp),    intent (in)     :: wsave(lensav)
            integer (ip), intent (in)     :: lensav
            integer (ip), intent (out)    :: ierror
            !--------------------------------------------------------------
        end subroutine rfft1f
        !
        !==> 2D real initialization
        !
        module subroutine rfft2i(l, m, wsave, lensav, ierror)
            !
            ! rfft2i: initialization for rfft2b and rfft2f.
            !
            !  purpose:
            !  rfft2i initializes real array wsave for use in its companion routines
            !  rfft2f and rfft2b for computing the two-dimensional fast fourier
            !  transform of real data.  prime factorizations of l and m, together with
            !  tabulations of the trigonometric functions, are computed and stored in
            !  array wsave.  rfft2i must be called prior to the first call to rfft2f
            !  or rfft2b.  separate wsave arrays are required for different values of
            !  l or m.
            !
            !
            !  integer l, the number of elements to be transformed
            !  in the first dimension.  the transform is most efficient when l is a
            !  product of small primes.
            !
            !  integer m, the number of elements to be transformed
            !  in the second dimension.  the transform is most efficient when m is a
            !  product of small primes.
            !
            !  integer lensav, the number of elements in the wsave
            !  array.  lensav must be at least l + m + int(log(real(l)))
            !  + int(log(real(m))) + 8.
            !
            !  real wsave(lensav), containing the prime factors
            !  of l and m, and also containing certain trigonometric values which
            !  will be used in routines rfft2b or rfft2f.
            !
            !  integer ier, error_flag.
            !  0, successful exit;
            !  2, input parameter lensav not big enough;
            !  20, input error returned by lower level routine.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)  :: l
            integer (ip), intent (in)  :: m
            real (wp),    intent (out) :: wsave(lensav)
            integer (ip), intent (in)  :: lensav
            integer (ip), intent (out) :: ierror
            !--------------------------------------------------------------
        end subroutine rfft2i
        !
        !==>  2D real backward
        !
        module subroutine rfft2b(ldim, l, m, r, wsave, lensav, work, lenwrk, ierror)
            !
            ! rfft2b: real backward fast fourier transform, 2d.
            !
            ! purpose:
            !
            !  computes the two-dimensional discrete fourier transform of the
            !  complex fourier coefficients a real periodic array.  this transform is
            !  known as the backward transform or fourier synthesis, transforming from
            !  spectral to physical space.  routine rfft2b is normalized: a call to
            !  rfft2b followed by a call to rfft2f (or vice-versa) reproduces the
            !  original array within roundoff error.
            !
            !  parameters:
            !
            !  integer ldim, the first dimension of the 2d real
            !  array r, which must be at least 2*(l/2+1).
            !
            !  integer l, the number of elements to be transformed
            !  in the first dimension of the two-dimensional real array r.  the value of
            !  l must be less than or equal to that of ldim.  the transform is most
            !  efficient when l is a product of small primes.
            !
            !  integer m, the number of elements to be transformed
            !  in the second dimension of the two-dimensional real array r.  the transform
            !  is most efficient when m is a product of small primes.
            !
            !  input/real r(ldim,m), the real array of two
            !  dimensions.  on input, r contains the l/2+1-by-m complex subarray of
            !  spectral coefficients, on output, the physical coefficients.
            !
            !  input, real (wp) wsave(lensav).  wsave's contents must be
            !  initialized with a call to rfft2i before the first call to routine rfft2f
            !  or rfft2b with lengths l and m.  wsave's contents may be re-used for
            !  subsequent calls to rfft2f and rfft2b with the same transform lengths
            !  l and m.
            !
            !  integer lensav, the number of elements in the wsave
            !  array.  lensav must be at least l + m + int(log(real(l)))
            !  + int(log(real(m))) + 8.
            !
            !  workspace, real (wp) work(lenwrk).  work provides workspace, and
            !  its contents need not be saved between calls to routines rfft2b and rfft2f.
            !
            !  integer  lenwrk, the number of elements in the work
            !  array.  lenwrk must be at least ldim*m.
            !
            !  integer ierror, the error_flag.
            !  0, successful exit;
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough;
            !  6, input parameter ldim < 2*(l/2+1);
            !  20, input error returned by lower level routine.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)     :: ldim
            integer (ip), intent (in)     :: l
            integer (ip), intent (in)     :: m
            real (wp),    intent (in out) :: r(ldim,m)
            real (wp),    intent (in)     :: wsave(lensav)
            integer (ip), intent (in)     :: lensav
            real (wp),    intent (out)    :: work(lenwrk)
            integer (ip), intent (in)     :: lenwrk
            integer (ip), intent (out)    :: ierror
            !--------------------------------------------------------------
        end subroutine rfft2b
        !
        !==> 2D real forward
        !
        module subroutine rfft2f(ldim, l, m, r, wsave, lensav, work, lenwrk, ierror)
            !
            ! rfft2f: real forward fast fourier transform, 2d.
            !
            ! purpose:
            !
            !  rfft2f computes the two-dimensional discrete fourier transform of a
            !  real periodic array.  this transform is known as the forward transform
            !  or fourier analysis, transforming from physical to spectral space.
            !  routine rfft2f is normalized: a call to rfft2f followed by a call to
            !  rfft2b (or vice-versa) reproduces the original array within roundoff
            !  error.
            !
            !  parameters:
            !
            !  integer ldim, the first dimension of the 2d real
            !  array r, which must be at least 2*(l/2+1).
            !
            !  integer l, the number of elements to be transformed
            !  in the first dimension of the two-dimensional real array r.  the value
            !  of l must be less than or equal to that of ldim.  the transform is most
            !  efficient when l is a product of small primes.
            !
            !  integer m, the number of elements to be transformed
            !  in the second dimension of the two-dimensional real array r.  the
            !  transform is most efficient when m is a product of small primes.
            !
            !  input/real r(ldim,m), the real array of two
            !  dimensions.  on input, containing the l-by-m physical data to be
            !  transformed.  on output, the spectral coefficients.
            !
            !  input, real (wp) wsave(lensav).  wsave's contents must be
            !  initialized with a call to rfft2i before the first call to routine rfft2f
            !  or rfft2b with lengths l and m.  wsave's contents may be re-used for
            !  subsequent calls to rfft2f and rfft2b with the same transform lengths.
            !
            !  integer lensav, the number of elements in the wsave
            !  array.  lensav must be at least l + m + int(log(real(l)))
            !  + int(log(real(m))) + 8.
            !
            !  workspace, real (wp) work(lenwrk), provides workspace, and its
            !  contents need not be saved between calls to routines rfft2f and rfft2b.
            !
            !  integer lenwrk, the number of elements in the work
            !  array.  lenwrk must be at least ldim*m.
            !
            !  integer ierror, the error_flag.
            !  0, successful exit;
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough;
            !  6, input parameter ldim < 2*(l+1);
            !  20, input error returned by lower level routine.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)     :: ldim
            integer (ip), intent (in)     :: l
            integer (ip), intent (in)     :: m
            real (wp),    intent (in out) :: r(ldim,m)
            real (wp),    intent (in)     :: wsave(lensav)
            integer (ip), intent (in)     :: lensav
            real (wp),    intent (out)    :: work(lenwrk)
            integer (ip), intent (in)     :: lenwrk
            integer (ip), intent (out)    :: ierror
            !--------------------------------------------------------------
        end subroutine rfft2f
        !
        !==> multiple real initialization
        !
        module subroutine rfftmi(n, wsave, lensav, ierror)
            !
            ! rfftmi: initialization for rfftmb and rfftmf.
            !
            !  purpose:
            !
            !  rfftmi initializes array wsave for use in its companion routines
            !  rfftmb and rfftmf.  the prime factorization of n together with a
            !  tabulation of the trigonometric functions are computed and stored
            !  in array wsave.  separate wsave arrays are required for different
            !  values of n.
            !
            !  input
            !
            !  integer n, the length of each sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least n + int(log(real(n))) + 4.
            !
            !  output
            !  real wsave(lensav), work array containing the prime
            !  factors of n and also containing certain trigonometric
            !  values which will be used in routines rfftmb or rfftmf.
            !
            !  integer ierror, error_flag.
            !  0, successful exit;
            !  2, input parameter lensav not big enough.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)  :: n
            integer (ip), intent (in)  :: lensav
            real (wp),    intent (out) :: wsave(lensav)
            integer (ip), intent (out) :: ierror
            !--------------------------------------------------------------
        end subroutine rfftmi
        !
        !==> multiple real backward
        !
        module subroutine rfftmb(lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ierror)
            !
            ! rfftmb: real backward fft, 1d, multiple vectors.
            !
            ! purpose:
            !
            !  computes the one-dimensional fourier transform of multiple
            !  periodic sequences within a real array.  this transform is referred
            !  to as the backward transform or fourier synthesis, transforming the
            !  sequences from spectral to physical space.
            !
            !  this transform is normalized since a call to rfftmb followed
            !  by a call to rfftmf (or vice-versa) reproduces the original
            !  array  within roundoff error.
            !
            !  parameters:
            !
            !  integer lot, the number of sequences to be transformed
            !  within array r.
            !
            !  integer jump, the increment between the locations, in
            !  array r, of the first elements of two consecutive sequences to be
            !  transformed.
            !
            !  integer n, the length of each sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer inc, the increment between the locations, in
            !  array r, of two consecutive elements within the same sequence.
            !
            !  input/real r(lenr), real array containing lot
            !  sequences, each having length n.  r can have any number of dimensions,
            !  but the total number of locations must be at least lenr.  on input, the
            !  spectral data to be transformed, on output the physical data.
            !
            !  integer lenr, the dimension of the r array.
            !  lenr must be at least (lot-1)*jump + inc*(n-1) + 1.
            !
            !  input, real (wp) wsave(lensav).  wsave's contents must be
            !  initialized with a call to rfftmi before the first call to routine rfftmf
            !  or rfftmb for a given transform length n.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must  be at least n + int(log(real(n))) + 4.
            !
            !  workspace, real (wp) work(lenwrk).
            !
            !  integer lenwrk, the dimension of the work array.
            !  lenwrk must be at least lot*n.
            !
            !  integer ierror, error_flag.
            !  0, successful exit;
            !  1, input parameter lenr not big enough;
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough;
            !  4, input parameters inc, jump, n, lot are not consistent.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)     :: lot
            integer (ip), intent (in)     :: jump
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: inc
            real (wp),    intent (in out) :: r(lenr)
            integer (ip), intent (in)     :: lenr
            real (wp),    intent (in)     :: wsave(lensav)
            integer (ip), intent (in)     :: lensav
            real (wp),    intent (out)    :: work(lenwrk)
            integer (ip), intent (in)     :: lenwrk
            integer (ip), intent (out)    :: ierror
            !--------------------------------------------------------------
        end subroutine rfftmb
        !
        !==> multiple real forward
        !
        module subroutine rfftmf(lot, jump, n, inc, r, lenr, wsave, lensav, work, lenwrk, ierror)
            !
            ! rfftmf: real forward fft, 1d, multiple vectors.
            !
            !  purpose:
            !
            !  rfftmf computes the one-dimensional fourier transform of multiple
            !  periodic sequences within a real array.  this transform is referred
            !  to as the forward transform or fourier analysis, transforming the
            !  sequences from physical to spectral space.
            !
            !  this transform is normalized since a call to rfftmf followed
            !  by a call to rfftmb (or vice-versa) reproduces the original array
            !  within roundoff error.
            !
            !  parameters:
            !
            !  integer lot, the number of sequences to be transformed
            !  within array r.
            !
            !  integer jump, the increment between the locations, in
            !  array r, of the first elements of two consecutive sequences to be
            !  transformed.
            !
            !  integer n, the length of each sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer inc, the increment between the locations,
            !  in array r, of two consecutive elements within the same sequence.
            !
            !  input/real r(lenr), real array containing lot
            !  sequences, each having length n.  r can have any number of dimensions, but
            !  the total number of locations must be at least lenr.  on input, the
            !  physical data to be transformed, on output the spectral data.
            !
            !  integer lenr, the dimension of the r array.
            !  lenr must be at least (lot-1)*jump + inc*(n-1) + 1.
            !
            !  input, real (wp) wsave(lensav).  wsave's contents must be
            !  initialized with a call to rfftmi before the first call to routine rfftmf
            !  or rfftmb for a given transform length n.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least n + int(log(real(n))) + 4.
            !
            !  workspace, real (wp) work(lenwrk).
            !
            !  integer lenwrk, the dimension of the work array.
            !  lenwrk must be at least lot*n.
            !
            !  integer ierror, error_flag.
            !  0, successful exit;
            !  1, input parameter lenr not big enough;
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough;
            !  4, input parameters inc, jump, n, lot are not consistent.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)     :: lot
            integer (ip), intent (in)     :: jump
            integer (ip), intent (in)     :: n
            integer (ip), intent (in)     :: inc
            real (wp),    intent (in out) :: r(lenr)
            integer (ip), intent (in)     :: lenr
            real (wp),    intent (in)     :: wsave(lensav)
            integer (ip), intent (in)     :: lensav
            real (wp),    intent (out)    :: work(lenwrk)
            integer (ip), intent (in)     :: lenwrk
            integer (ip), intent (out)    :: ierror
            !--------------------------------------------------------------
        end subroutine rfftmf
    end interface


    !---------------------------------------------------------------------------------
    ! Variables confined to the module
    !---------------------------------------------------------------------------------
    real (wp), parameter :: ZERO = 0.0_wp
    real (wp), parameter :: HALF = 0.5_wp
    real (wp), parameter :: ONE = 1.0_wp
    real (wp), parameter :: TWO = 2.0_wp
    real (wp), parameter :: THREE = 3.0_wp
    real (wp), parameter :: FOUR = 4.0_wp
    real (wp), parameter :: FIVE = 5.0_wp
    !---------------------------------------------------------------------------------


end module real_transform_routines
