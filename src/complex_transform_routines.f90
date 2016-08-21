module complex_transform_routines

    use auxiliary_routines

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: cfft1i, cfft1b, cfft1f, cfft2i, cfft2b, cfft2f, cfftmi, cfftmb, cfftmf

    interface
        !
        !==> 1D complex initialization
        !
        module subroutine cfft1i(n, wsave, lensav, ierror)
            !
            ! cfft1i: initialization for cfft1b and cfft1f.
            !
            !  Purpose:
            !
            !  cfft1i initializes array wsave for use in its companion routines
            !  cfft1b and cfft1f. Routine cfft1i must be called before the first
            !  call to cfft1b or cfft1f, and after whenever the value of integer
            !  n changes.
            !
            !  Parameters:
            !
            !  input,
            !  n, the length of the sequence to be
            !  transformed.  the transform is most efficient when n is a product
            !  of small primes.
            !
            !  input
            !  lensav, the dimension of the wsave array.
            !  lensav must be at least 2*n + int(log(real(n))) + 4.
            !
            !  output,
            !  wsave(lensav), containing the prime factors
            !  of n and  also containing certain trigonometric values which will be used
            !  in routines cfft1b or cfft1f.
            !
            !  output
            !  ierror, error_flag.
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
        end subroutine cfft1i
        !
        !==> 1D complex backward
        !
        module subroutine cfft1b(n, inc, c, lenc, wsave, lensav, work, lenwrk, ierror)

            use, intrinsic :: ISO_C_binding, only: c_f_pointer, c_loc
            !
            !  input
            !  integer n, the length of the sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer inc, the increment between the locations, in
            !  array c, of two consecutive elements within the sequence to be transformed.
            !
            !  integer lenc, the dimension of the c array.
            !  lenc must be at least inc*(n-1) + 1.
            !
            !  real wsave(lensav). wsave's contents must be initialized with a call
            !  to cfft1i before the first call to routine cfft1f
            !  or cfft1b for a given transform length n.  wsave's contents may be
            !  re-used for subsequent calls to cfft1f and cfft1b with the same n.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least 2*n + int(log(real(n))) + 4.
            !
            !
            !  input lenwrk, the dimension of the work array.
            !  lenwrk must be at least 2*n.
            !
            !  input/output
            !  complex c(lenc) containing the sequence to be
            !  transformed.
            !
            !  real workspace work(lenwrk).
            !
            !  integer ier, error_flag.
            !  0, successful exit;
            !  1, input parameter lenc not big enough;
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough;
            !  20, input error returned by lower level routine.
            !
            !--------------------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------------------
            integer (ip), intent (in)             :: n
            integer (ip), intent (in)             :: inc
            complex (wp), intent (in out), target :: c(lenc)
            integer (ip), intent (in)             :: lenc
            integer (ip), intent (in)             :: lensav
            integer (ip), intent (in)             :: lenwrk
            integer (ip), intent (out)            :: ierror
            real (wp),    intent (in out)         :: work(lenwrk)
            real (wp),    intent (in out)         :: wsave(lensav)
            !--------------------------------------------------------------
        end subroutine cfft1b
            !
            !==> 1D complex forward
            !
        module subroutine cfft1f(n, inc, c, lenc, wsave, lensav, work, lenwrk, ierror)
            use, intrinsic :: ISO_C_binding, only: c_f_pointer, c_loc
            !
            ! cfft1f: complex forward fast fourier transform, 1d.
            !
            !  purpose:
            !
            !  cfft1f computes the one-dimensional fourier transform of a single
            !  periodic sequence within a complex array. this transform is referred
            !  to as the forward transform or fourier analysis, transforming the
            !  sequence from physical to spectral space.
            !
            !  this transform is normalized since a call to cfft1f followed
            !  by a call to cfft1b (or vice-versa) reproduces the original
            !  array within roundoff error.
            !
            !  input
            !
            !  integer n, the length of the sequence to be
            !  transformed. the transform is most efficient when
            !  n is a product of small primes.
            !
            !  integer inc, the increment between the locations, in
            !  array c, of two consecutive elements within the sequence to be transformed.
            !
            !  real wsave(lensav).  wsave's contents must be
            !  initialized with a call to cfft1i before the first call to routine cfft1f
            !  or cfft1b for a given transform length n.  wsave's contents may be re-used
            !  for subsequent calls to cfft1f and cfft1b with the same n.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least 2*n + int(log(real(n))) + 4.
            !
            !  integer lenwrk, the dimension of the work array.
            !  lenwrk must be at least 2*n.
            !
            !  input/output
            !  complex c(lenc) containing the sequence to be transformed.
            !
            !  real work(lenwrk), workspace array
            !  integer lenc, the dimension of the c array.
            !  lenc must be at least inc*(n-1) + 1.
            !
            !  output
            !  integer ier, error_flag.
            !  0, successful exit;
            !  1, input parameter lenc not big enough;
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough;
            !  20, input error returned by lower level routine.
            !
            !------------------------------------------------------------------
            ! Dummy arguments
            !------------------------------------------------------------------
            integer (ip), intent (in)             :: n
            integer (ip), intent (in)             :: inc
            complex (wp), intent (in out), target :: c(lenc)
            integer (ip), intent (in)             :: lenc
            real (wp),    intent (in)             :: wsave(lensav)
            integer (ip), intent (in)             :: lensav
            real (wp),    intent (out)            :: work(lenwrk)
            integer (ip), intent (in)             :: lenwrk
            integer (ip), intent (out)            :: ierror
            !------------------------------------------------------------------
        end subroutine cfft1f
        !
        !==> 2D complex initialization
        !
        module subroutine cfft2i(l, m, wsave, lensav, ierror)
            !
            ! cfft2i: initialization for cfft2b and cfft2f.
            !
            !  purpose:
            !
            !  cfft2i initializes real array wsave for use in its companion
            !  routines cfft2f and cfft2b for computing two-dimensional fast
            !  fourier transforms of complex data.  prime factorizations of l and m,
            !  together with tabulations of the trigonometric functions, are
            !  computed and stored in array wsave.
            !
            !  on 10 may 2010, this code was modified by changing the value
            !  of an index into the wsave array.
            !
            !  input
            !  integer l, the number of elements to be transformed
            !  in the first dimension. the transform is most efficient when l is a
            !  product of small primes.
            !
            !  integer m, the number of elements to be transformed
            !  in the second dimension. the transform is most efficient when m is a
            !  product of small primes.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least
            !
            !  2*(l+m) + int(log(real(l)))+ int(log(real(m))) + 8.
            !
            !  output
            !  real wsave(lensav), contains the prime factors of l
            !  and m, and also certain trigonometric values which will be used in
            !  routines cfft2b or cfft2f.
            !
            !  integer  ierror, error_flag.
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
        end subroutine cfft2i
        !
        !==> 2D complex backward
        !
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
        end subroutine cfft2b
        !
        !==> 2D complex forward
        !
        module subroutine cfft2f(ldim, l, m, c, wsave, lensav, work, lenwrk, ierror)
            !
            ! cfft2f: complex forward fast fourier transform, 2d.
            !
            !  Purpose:
            !
            !  cfft2f computes the two-dimensional discrete fourier transform of
            !  a complex periodic array. This transform is known as the forward
            !  transform or fourier analysis, transforming from physical to
            !  spectral space. routine cfft2f is normalized, in that a call to
            !  cfft2f followed by a call to cfft2b (or vice-versa) reproduces the
            !  original array within roundoff error.
            !
            !  BUG FIX
            !  On 10 May 2010, this code was modified by changing the value
            !  of an index into the wsave array.
            !
            !  INPUT
            !
            !  integer ldim, the first dimension of the array c.
            !
            !  integer l, the number of elements to be transformed
            !  in the first dimension of the two-dimensional complex array c.  the value
            !  of l must be less than or equal to that of ldim.  the transform is most
            !  efficient when l is a product of small primes.
            !
            !  integer m, the number of elements to be transformed
            !  in the second dimension of the two-dimensional complex array c.  the
            !  transform is most efficient when m is a product of small primes.
            !
            !  real wsave(lensav). wsave's contents must be
            !  initialized with a call to cfft2i before the first call to routine cfft2f
            !  or cfft2b with transform lengths l and m.  wsave's contents may be re-used
            !  for subsequent calls to cfft2f and cfft2b having those same
            !  transform lengths.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least
            !
            !  2*(l+m) + int(log(real(l))) + int(log(real(m))) + 8.
            !
            !
            !  integer lenwrk, the dimension of the work array.
            !  lenwrk must be at least 2*l*m.
            !
            !  INPUT/OUTPUT
            !  complex c(ldim,m), on input, the array of two
            !  dimensions containing the (l,m) subarray to be transformed. On output, the
            !  transformed data.
            !
            !  real work(lenwrk), workspace array
            !
            !  OUTPUT
            !  integer ierror, error_flag.
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
        end subroutine cfft2f
        !
        !==> multiple complex initialization
        !
        module subroutine cfftmi(n, wsave, lensav, ierror)
            !
            !  cfftmi: initialization for cfftmb and cfftmf.
            !
            !  purpose:
            !
            !  cfftmi initializes array wsave for use in its companion routines
            !  cfftmb and cfftmf.  cfftmi must be called before the first call
            !  to cfftmb or cfftmf, and after whenever the value of integer n changes.
            !
            !  parameters:
            !
            !  integer n, the length of each sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least 2*n + int(log(real(n))) + 4.
            !
            !  real wsave(lensav), containing the prime factors
            !  of n and also containing certain trigonometric values which will be used in
            !  routines cfftmb or cfftmf.
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
        end subroutine cfftmi
        !
        !==> multiple complex backward
        !
        module subroutine cfftmb(lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ierror)
            !
            ! cfftmb: complex backward fft, 1d, multiple vectors.
            !
            !  purpose:
            !
            !  cfftmb computes the one-dimensional fourier transform of multiple
            !  periodic sequences within a complex array.  this transform is referred
            !  to as the backward transform or fourier synthesis, transforming the
            !  sequences from spectral to physical space.  this transform is
            !  normalized since a call to cfftmf followed by a call to cfftmb (or
            !  vice-versa) reproduces the original array within roundoff error.
            !
            !  the parameters inc, jump, n and lot are consistent if equality
            !  i1*inc + j1*jump = i2*inc + j2*jump for i1,i2 < n and j1,j2 < lot
            !  implies i1=i2 and j1=j2.  for multiple ffts to execute correctly,
            !  input variables inc, jump, n and lot must be consistent, otherwise
            !  at least one array element mistakenly is transformed more than once.
            !
            !  parameters:
            !
            !  integer lot, the number of sequences to be transformed
            !  within array c.
            !
            !  integer jump, the increment between the locations, in
            !  array c, of the first elements of two consecutive sequences to be
            !  transformed.
            !
            !  integer n, the length of each sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer inc, the increment between the locations, in
            !  array c, of two consecutive elements within the same sequence to be
            !  transformed.
            !
            !  input/output, complex (wp) c(lenc), an array containing lot
            !  sequences, each having length n, to be transformed.  c can have any
            !  number of dimensions, but the total number of locations must be at least
            !  lenc.  on output, c contains the transformed sequences.
            !
            !  integer lenc, the dimension of the c array.
            !  lenc must be at least (lot-1)*jump + inc*(n-1) + 1.
            !
            !  input, real (wp) wsave(lensav).  wsave's contents must be
            !  initialized with a call to cfftmi before the first call to routine cfftmf
            !  or cfftmb for a given transform length n.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least 2*n + int(log(real(n))) + 4.
            !
            !  workspace, real (wp) work(lenwrk).
            !
            !  integer lenwrk, the dimension of the work array.
            !  lenwrk must be at least 2*lot*n.
            !
            !  integer ierror, error_flag.
            !  0, successful exit
            !  1, input parameter lenc not big enough;
            !  2, input parameter lensav not big enough;
            !  3, input parameter lenwrk not big enough;
            !  4, input parameters inc, jump, n, lot are not consistent.
            !
            !--------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------
            integer (ip) lot
            integer (ip) jump
            integer (ip) n
            integer (ip) inc
            real (wp) c(2,lenc)
            integer (ip) lenc
            real (wp) wsave(lensav)
            integer (ip) lensav
            real (wp) work(lenwrk)
            integer (ip) lenwrk
            integer (ip) ierror
            !--------------------------------------------------
        end subroutine cfftmb
        !
        !==> multiple complex forward
        !
        module subroutine cfftmf(lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ierror)
            !
            ! cfftmf: complex forward fft, 1d, multiple vectors.
            !
            !  Purpose:
            !
            !  cfftmf computes the one-dimensional fourier transform of multiple
            !  periodic sequences within a complex array. this transform is referred
            !  to as the forward transform or fourier analysis, transforming the
            !  sequences from physical to spectral space. this transform is
            !  normalized since a call to cfftmf followed by a call to cfftmb
            !  (or vice-versa) reproduces the original array within roundoff error.
            !
            !  the parameters integers inc, jump, n and lot are consistent if equality
            !  i1*inc + j1*jump = i2*inc + j2*jump for i1,i2 < n and j1,j2 < lot
            !  implies i1=i2 and j1=j2. for multiple ffts to execute correctly,
            !  input variables inc, jump, n and lot must be consistent, otherwise
            !  at least one array element mistakenly is transformed more than once.
            !
            !  parameters:
            !
            !  integer lot, the number of sequences to be
            !  transformed within array c.
            !
            !  integer jump, the increment between the locations,
            !  in array c, of the first elements of two consecutive sequences to be
            !  transformed.
            !
            !  integer n, the length of each sequence to be
            !  transformed.  the transform is most efficient when n is a product of
            !  small primes.
            !
            !  integer inc, the increment between the locations, in
            !  array c, of two consecutive elements within the same sequence to be
            !  transformed.
            !
            !  input/output, complex (wp) c(lenc), array containing lot sequences,
            !  each having length n, to be transformed.  c can have any number of
            !  dimensions, but the total number of locations must be at least lenc.
            !
            !  integer lenc, the dimension of the c array.
            !  lenc must be at least (lot-1)*jump + inc*(n-1) + 1.
            !
            !  input, real (wp) wsave(lensav).  wsave's contents must be
            !  initialized with a call to cfftmi before the first call to routine cfftmf
            !  or cfftmb for a given transform length n.
            !
            !  integer lensav, the dimension of the wsave array.
            !  lensav must be at least 2*n + int(log(real(n))) + 4.
            !
            !  workspace, real (wp) work(lenwrk).
            !
            !  integer lenwrk, the dimension of the work array.
            !  lenwrk must be at least 2*lot*n.
            !
            !  integer ierror, error_flag.
            !  0 successful exit;
            !  1 input parameter lenc not big enough;
            !  2 input parameter lensav not big enough;
            !  3 input parameter lenwrk not big enough;
            !  4 input parameters inc, jump, n, lot are not consistent.
            !
            !--------------------------------------------------
            ! Dummy arguments
            !--------------------------------------------------
            integer (ip) lot
            integer (ip) jump
            integer (ip) n
            integer (ip) inc
            real(wp) c(2,lenc)
            integer (ip) lenc
            real (wp) wsave(lensav)
            integer (ip) lensav
            real (wp) work(lenwrk)
            integer (ip) lenwrk
            integer (ip) ierror
            !--------------------------------------------------
        end subroutine cfftmf
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

end module complex_transform_routines
