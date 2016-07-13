module fftpack_library

    use fftpack_precision, only: &
        wp, & ! working precision
        ip ! integer precision

    use complex_transform_routines

    use real_transform_routines

    use cosine_transform_routines

    use sine_transform_routines

    use quarter_cosine_transform_routines

    use quarter_sine_transform_routines

    use type_FFTpack, only: &
        FFTpack

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: wp, ip
    public :: FFTpack
    public :: cfft1i, cfft1b, cfft1f, cfft2i, cfft2b, cfft2f, cfftmi, cfftmb, cfftmf
    public :: rfft1i, rfft1b, rfft1f, rfft2i, rfft2b, rfft2f, rfftmi, rfftmb, rfftmf
    public :: cost1i, cost1b, cost1f, costmi, costmb, costmf
    public :: sint1i, sint1b, sint1f, sintmi, sintmb, sintmf
    public :: cosq1i, cosq1b, cosq1f, cosqmi, cosqmb, cosqmf
    public :: sinq1i, sinq1b, sinq1f, sinqmi, sinqmb, sinqmf

end module fftpack_library
