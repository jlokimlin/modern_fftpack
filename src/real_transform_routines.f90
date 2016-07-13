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

contains



    subroutine mrftb1(m,im,n,in,c,ch,wa,fac)

        integer (ip) in
        integer (ip) m
        integer (ip) n

        real (wp) c(in,*)
        real (wp) ch(m,*)
        real (wp) fac(15)
        real (wp), parameter :: NEG_HALF = -HALF
        integer (ip) i
        integer (ip) idl1
        integer (ip) ido
        integer (ip) im
        integer (ip) iip
        integer (ip) iw
        integer (ip) ix2
        integer (ip) ix3
        integer (ip) ix4
        integer (ip) j
        integer (ip) k1
        integer (ip) l1
        integer (ip) l2
        integer (ip) m2
        integer (ip) modn
        integer (ip) na
        integer (ip) nf
        integer (ip) nl
        real (wp) wa(n)

        nf = int(fac(2), kind=ip)
        na = 0

        do k1=1,nf
            iip = int(fac(k1+2), kind=ip)
            na = 1-na
            if (iip <= 5) then
                cycle
            end if
            if (k1 == nf) then
                cycle
            end if
            na = 1-na
        end do


        modn = mod(n,2)

        if (modn /= 0) then
            nl = n-1
        else
            nl = n-2
        end if

        if (na /= 0) then
            m2 = 1-im
            do i=1,m
                m2 = m2+im
                ch(i,1) = c(m2,1)
                ch(i,n) = c(m2,n)
            end do
            do j=2,nl,2
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    ch(i,j) = HALF*c(m2,j)
                    ch(i,j+1) = NEG_HALF*c(m2,j+1)
                end do
            end do
        else
            do j=2,nl,2
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,j) = HALF*c(m2,j)
                    c(m2,j+1) = NEG_HALF*c(m2,j+1)
                end do
            end do
        end if

        l1 = 1
        iw = 1
        do k1=1,nf
            iip = int(fac(k1+2), kind=ip)
            l2 = iip*l1
            ido = n/l2
            idl1 = ido*l1

            select case (iip)
                case (2)
                    if (na == 0) then
                        call mradb2 (m,ido,l1,c,im,in,ch,1,m,wa(iw))
                    else
                        call mradb2 (m,ido,l1,ch,1,m,c,im,in,wa(iw))
                    end if
                    na = 1-na
                case (3)
                    ix2 = iw+ido
                    if (na == 0) then
                        call mradb3(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2))
                    else
                        call mradb3(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2))
                    end if
                    na = 1-na
                case(4)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    if (na == 0) then
                        call mradb4 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3))
                    else
                        call mradb4 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3))
                    end if
                    na = 1-na
                case (5)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    ix4 = ix3+ido
                    if (na == 0) then
                        call mradb5 (m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    else
                        call mradb5 (m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    end if
                    na = 1-na
                case default
                    if (na == 0) then
                        call mradbg (m,ido,iip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
                    else
                        call mradbg (m,ido,iip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
                    end if
                    if (ido == 1) then
                        na = 1-na
                    end if
            end select
            l1 = l2
            iw = iw+(iip-1)*ido
        end do

    contains

        subroutine mradb2(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,2,l1)
            real (wp) ch(in2,ido,l1,2)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,k,1) = cc(m1,1,1,k)+cc(m1,ido,2,k)
                    ch(m2,1,k,2) = cc(m1,1,1,k)-cc(m1,ido,2,k)
                end do
            end do

            if (ido < 2) then
                return
            else if (ido == 2) then
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,ido,k,1) = cc(m1,ido,1,k)+cc(m1,ido,1,k)
                        ch(m2,ido,k,2) = -(cc(m1,1,2,k)+cc(m1,1,2,k))
                    end do
                end do
            else
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+cc(m1,ic-1,2,k)
                            ch(m2,i,k,1) = cc(m1,i,1,k)-cc(m1,ic,2,k)
                            ch(m2,i-1,k,2) = wa1(i-2)*(cc(m1,i-1,1,k)-cc(m1,ic-1,2,k)) &
                                -wa1(i-1)*(cc(m1,i,1,k)+cc(m1,ic,2,k))

                            ch(m2,i,k,2) = wa1(i-2)*(cc(m1,i,1,k)+cc(m1,ic,2,k))+wa1(i-1) &
                                *(cc(m1,i-1,1,k)-cc(m1,ic-1,2,k))
                        end do
                    end do
                end do
                if (mod(ido,2) /= 1) then
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,ido,k,1) = cc(m1,ido,1,k)+cc(m1,ido,1,k)
                            ch(m2,ido,k,2) = -(cc(m1,1,2,k)+cc(m1,1,2,k))
                        end do
                    end do
                end if
            end if

        end subroutine mradb2

        subroutine mradb3 (m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1


            real (wp) cc(in1,ido,3,l1)
            real (wp) ch(in2,ido,l1,3)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
            real (wp), parameter :: ARG= TWO_PI/3
            real (wp), parameter :: TAUI = cos(ARG)
            real (wp), parameter :: TAUR = sin(ARG)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,k,1) = cc(m1,1,1,k)+ TWO *cc(m1,ido,2,k)
                    ch(m2,1,k,2) = cc(m1,1,1,k)+( TWO *TAUR)*cc(m1,ido,2,k) &
                        -( TWO *TAUI)*cc(m1,1,3,k)
                    ch(m2,1,k,3) = cc(m1,1,1,k)+( TWO *TAUR)*cc(m1,ido,2,k) &
                        + TWO *TAUI*cc(m1,1,3,k)
                end do
            end do

            if (ido /= 1) then
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))
                            ch(m2,i,k,1) = cc(m1,i,1,k)+(cc(m1,i,3,k)-cc(m1,ic,2,k))

                            ch(m2,i-1,k,2) = wa1(i-2)* &
                                ((cc(m1,i-1,1,k)+TAUR*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))- &
                                (TAUI*(cc(m1,i,3,k)+cc(m1,ic,2,k)))) &
                                -wa1(i-1)* &
                                ((cc(m1,i,1,k)+TAUR*(cc(m1,i,3,k)-cc(m1,ic,2,k)))+ &
                                (TAUI*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))))

                            ch(m2,i,k,2) = wa1(i-2)* &
                                ((cc(m1,i,1,k)+TAUR*(cc(m1,i,3,k)-cc(m1,ic,2,k)))+ &
                                (TAUI*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))) &
                                +wa1(i-1)* &
                                ((cc(m1,i-1,1,k)+TAUR*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))- &
                                (TAUI*(cc(m1,i,3,k)+cc(m1,ic,2,k))))

                            ch(m2,i-1,k,3) = wa2(i-2)* &
                                ((cc(m1,i-1,1,k)+TAUR*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))+ &
                                (TAUI*(cc(m1,i,3,k)+cc(m1,ic,2,k)))) &
                                -wa2(i-1)* &
                                ((cc(m1,i,1,k)+TAUR*(cc(m1,i,3,k)-cc(m1,ic,2,k)))- &
                                (TAUI*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))))

                            ch(m2,i,k,3) = wa2(i-2)* &
                                ((cc(m1,i,1,k)+TAUR*(cc(m1,i,3,k)-cc(m1,ic,2,k)))- &
                                (TAUI*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))) &
                                +wa2(i-1)* &
                                ((cc(m1,i-1,1,k)+TAUR*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))+ &
                                (TAUI*(cc(m1,i,3,k)+cc(m1,ic,2,k))))
                        end do
                    end do
                end do
            end if

        end subroutine mradb3

        subroutine mradb4(m, ido, l1, cc, im1, in1, ch, im2, in2, wa1, wa2, wa3)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,4,l1)
            real (wp) ch(in2,ido,l1,4)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp), parameter :: SQRT2 = sqrt(TWO)
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp) wa3(ido)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,k,3) = (cc(m1,1,1,k)+cc(m1,ido,4,k)) &
                        -(cc(m1,ido,2,k)+cc(m1,ido,2,k))
                    ch(m2,1,k,1) = (cc(m1,1,1,k)+cc(m1,ido,4,k)) &
                        +(cc(m1,ido,2,k)+cc(m1,ido,2,k))
                    ch(m2,1,k,4) = (cc(m1,1,1,k)-cc(m1,ido,4,k)) &
                        +(cc(m1,1,3,k)+cc(m1,1,3,k))
                    ch(m2,1,k,2) = (cc(m1,1,1,k)-cc(m1,ido,4,k)) &
                        -(cc(m1,1,3,k)+cc(m1,1,3,k))
                end do
            end do

            if (ido < 2) then
                return
            else if (ido == 2) then
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,ido,k,1) = (cc(m1,ido,1,k)+cc(m1,ido,3,k)) &
                            +(cc(m1,ido,1,k)+cc(m1,ido,3,k))
                        ch(m2,ido,k,2) = SQRT2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
                            -(cc(m1,1,2,k)+cc(m1,1,4,k)))
                        ch(m2,ido,k,3) = (cc(m1,1,4,k)-cc(m1,1,2,k)) &
                            +(cc(m1,1,4,k)-cc(m1,1,2,k))
                        ch(m2,ido,k,4) = -SQRT2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
                            +(cc(m1,1,2,k)+cc(m1,1,4,k)))
                    end do
                end do
            else
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,k,1) = (cc(m1,i-1,1,k)+cc(m1,ic-1,4,k)) &
                                +(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))
                            ch(m2,i,k,1) = (cc(m1,i,1,k)-cc(m1,ic,4,k)) &
                                +(cc(m1,i,3,k)-cc(m1,ic,2,k))
                            ch(m2,i-1,k,2)=wa1(i-2)*((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k)) &
                                -(cc(m1,i,3,k)+cc(m1,ic,2,k)))-wa1(i-1) &
                                *((cc(m1,i,1,k)+cc(m1,ic,4,k))+(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))
                            ch(m2,i,k,2)=wa1(i-2)*((cc(m1,i,1,k)+cc(m1,ic,4,k)) &
                                +(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))+wa1(i-1) &
                                *((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k))-(cc(m1,i,3,k)+cc(m1,ic,2,k)))
                            ch(m2,i-1,k,3)=wa2(i-2)*((cc(m1,i-1,1,k)+cc(m1,ic-1,4,k)) &
                                -(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)))-wa2(i-1) &
                                *((cc(m1,i,1,k)-cc(m1,ic,4,k))-(cc(m1,i,3,k)-cc(m1,ic,2,k)))
                            ch(m2,i,k,3)=wa2(i-2)*((cc(m1,i,1,k)-cc(m1,ic,4,k)) &
                                -(cc(m1,i,3,k)-cc(m1,ic,2,k)))+wa2(i-1) &
                                *((cc(m1,i-1,1,k)+cc(m1,ic-1,4,k))-(cc(m1,i-1,3,k) &
                                +cc(m1,ic-1,2,k)))
                            ch(m2,i-1,k,4)=wa3(i-2)*((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k)) &
                                +(cc(m1,i,3,k)+cc(m1,ic,2,k)))-wa3(i-1) &
                                *((cc(m1,i,1,k)+cc(m1,ic,4,k))-(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))
                            ch(m2,i,k,4)=wa3(i-2)*((cc(m1,i,1,k)+cc(m1,ic,4,k)) &
                                -(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k)))+wa3(i-1) &
                                *((cc(m1,i-1,1,k)-cc(m1,ic-1,4,k))+(cc(m1,i,3,k)+cc(m1,ic,2,k)))
                        end do
                    end do
                end do
                if (mod(ido,2) /= 1) then
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,ido,k,1) = (cc(m1,ido,1,k)+cc(m1,ido,3,k)) &
                                +(cc(m1,ido,1,k)+cc(m1,ido,3,k))
                            ch(m2,ido,k,2) = SQRT2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
                                -(cc(m1,1,2,k)+cc(m1,1,4,k)))
                            ch(m2,ido,k,3) = (cc(m1,1,4,k)-cc(m1,1,2,k)) &
                                +(cc(m1,1,4,k)-cc(m1,1,2,k))
                            ch(m2,ido,k,4) = -SQRT2*((cc(m1,ido,1,k)-cc(m1,ido,3,k)) &
                                +(cc(m1,1,2,k)+cc(m1,1,4,k)))
                        end do
                    end do
                end if
            end if

        end subroutine mradb4

        subroutine mradb5(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3,wa4)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,5,l1)
            real (wp) ch(in2,ido,l1,5)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp) wa3(ido)
            real (wp) wa4(ido)

            real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
            real (wp), parameter :: ARG = TWO_PI/5
            real (wp), parameter :: TR11=cos(ARG)
            real (wp), parameter :: TI11=sin(ARG)
            real (wp), parameter :: TR12=cos(TWO*ARG)
            real (wp), parameter :: TI12=sin(TWO*ARG)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,k,1) = cc(m1,1,1,k)+ TWO *cc(m1,ido,2,k)&
                        + TWO *cc(m1,ido,4,k)
                    ch(m2,1,k,2) = (cc(m1,1,1,k)+TR11* TWO *cc(m1,ido,2,k) &
                        +TR12* TWO *cc(m1,ido,4,k))-(TI11* TWO *cc(m1,1,3,k) &
                        +TI12* TWO *cc(m1,1,5,k))
                    ch(m2,1,k,3) = (cc(m1,1,1,k)+TR12* TWO *cc(m1,ido,2,k) &
                        +TR11* TWO *cc(m1,ido,4,k))-(TI12* TWO *cc(m1,1,3,k) &
                        -TI11* TWO *cc(m1,1,5,k))
                    ch(m2,1,k,4) = (cc(m1,1,1,k)+TR12* TWO *cc(m1,ido,2,k) &
                        +TR11* TWO *cc(m1,ido,4,k))+(TI12* TWO *cc(m1,1,3,k) &
                        -TI11* TWO *cc(m1,1,5,k))
                    ch(m2,1,k,5) = (cc(m1,1,1,k)+TR11* TWO *cc(m1,ido,2,k) &
                        +TR12* TWO *cc(m1,ido,4,k))+(TI11* TWO *cc(m1,1,3,k) &
                        +TI12* TWO *cc(m1,1,5,k))
                end do
            end do

            if (ido /= 1) then
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,k,1) = cc(m1,i-1,1,k)+(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k))
                            ch(m2,i,k,1) = cc(m1,i,1,k)+(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
                                +(cc(m1,i,5,k)-cc(m1,ic,4,k))
                            ch(m2,i-1,k,2) = wa1(i-2)*((cc(m1,i-1,1,k)+TR11* &
                                (cc(m1,i-1,3,k)+cc(m1,ic-1,2,k))+TR12 &
                                *(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(TI11*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))+TI12*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
                                -wa1(i-1)*((cc(m1,i,1,k)+TR11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
                                +TR12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))+(TI11*(cc(m1,i-1,3,k) &
                                -cc(m1,ic-1,2,k))+TI12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

                            ch(m2,i,k,2) = wa1(i-2)*((cc(m1,i,1,k)+TR11*(cc(m1,i,3,k) &
                                -cc(m1,ic,2,k))+TR12*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                +(TI11*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))+TI12 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))+wa1(i-1) &
                                *((cc(m1,i-1,1,k)+TR11*(cc(m1,i-1,3,k) &
                                +cc(m1,ic-1,2,k))+TR12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k))) &
                                -(TI11*(cc(m1,i,3,k)+cc(m1,ic,2,k))+TI12 &
                                *(cc(m1,i,5,k)+cc(m1,ic,4,k))))
                            ch(m2,i-1,k,3) = wa2(i-2) &
                                *((cc(m1,i-1,1,k)+TR12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(TI12*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))-TI11*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
                                -wa2(i-1) &
                                *((cc(m1,i,1,k)+TR12*(cc(m1,i,3,k)- &
                                cc(m1,ic,2,k))+TR11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                +(TI12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-TI11 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

                            ch(m2,i,k,3) = wa2(i-2) &
                                *((cc(m1,i,1,k)+TR12*(cc(m1,i,3,k)- &
                                cc(m1,ic,2,k))+TR11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                +(TI12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-TI11 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
                                +wa2(i-1) &
                                *((cc(m1,i-1,1,k)+TR12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))-(TI12*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))-TI11*(cc(m1,i,5,k)+cc(m1,ic,4,k))))

                            ch(m2,i-1,k,4) = wa3(i-2) &
                                *((cc(m1,i-1,1,k)+TR12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(TI12*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))-TI11*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
                                -wa3(i-1) &
                                *((cc(m1,i,1,k)+TR12*(cc(m1,i,3,k)- &
                                cc(m1,ic,2,k))+TR11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                -(TI12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-TI11 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

                            ch(m2,i,k,4) = wa3(i-2) &
                                *((cc(m1,i,1,k)+TR12*(cc(m1,i,3,k)- &
                                cc(m1,ic,2,k))+TR11*(cc(m1,i,5,k)-cc(m1,ic,4,k))) &
                                -(TI12*(cc(m1,i-1,3,k)-cc(m1,ic-1,2,k))-TI11 &
                                *(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
                                +wa3(i-1) &
                                *((cc(m1,i-1,1,k)+TR12*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR11*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(TI12*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))-TI11*(cc(m1,i,5,k)+cc(m1,ic,4,k))))

                            ch(m2,i-1,k,5) = wa4(i-2) &
                                *((cc(m1,i-1,1,k)+TR11*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(TI11*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))+TI12*(cc(m1,i,5,k)+cc(m1,ic,4,k)))) &
                                -wa4(i-1) &
                                *((cc(m1,i,1,k)+TR11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
                                +TR12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))-(TI11*(cc(m1,i-1,3,k) &
                                -cc(m1,ic-1,2,k))+TI12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k))))

                            ch(m2,i,k,5) = wa4(i-2) &
                                *((cc(m1,i,1,k)+TR11*(cc(m1,i,3,k)-cc(m1,ic,2,k)) &
                                +TR12*(cc(m1,i,5,k)-cc(m1,ic,4,k)))-(TI11*(cc(m1,i-1,3,k) &
                                -cc(m1,ic-1,2,k))+TI12*(cc(m1,i-1,5,k)-cc(m1,ic-1,4,k)))) &
                                +wa4(i-1) &
                                *((cc(m1,i-1,1,k)+TR11*(cc(m1,i-1,3,k)+cc(m1,ic-1,2,k)) &
                                +TR12*(cc(m1,i-1,5,k)+cc(m1,ic-1,4,k)))+(TI11*(cc(m1,i,3,k) &
                                +cc(m1,ic,2,k))+TI12*(cc(m1,i,5,k)+cc(m1,ic,4,k))))
                        end do
                    end do
                end do
            end if

        end subroutine mradb5

        subroutine mradbg (m,ido,iip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

            integer (ip) idl1
            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) iip
            integer (ip) l1

            real (wp) ai1
            real (wp) ai2
            real (wp) ar1
            real (wp) ar1h
            real (wp) ar2
            real (wp) ar2h, arg
            real (wp) c1(in1,ido,l1,iip)
            real (wp) c2(in1,idl1,iip)
            real (wp) cc(in1,ido,iip,l1)
            real (wp) ch(in2,ido,l1,iip)
            real (wp) ch2(in2,idl1,iip)
            real (wp) dc2, dcp
            real (wp) ds2, dsp
            integer (ip) i
            integer (ip) ic
            integer (ip) idij
            integer (ip) idp2
            integer (ip) ik
            integer (ip) im1
            integer (ip) im2
            integer (ip) iipp2
            integer (ip) iipph
            integer (ip) is
            integer (ip) j
            integer (ip) j2
            integer (ip) jc
            integer (ip) k
            integer (ip) l
            integer (ip) lc
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) nbd
            real (wp) wa(ido)

            real (wp), parameter :: TWO_PI = TWO * acos(-ONE)

            arg = TWO_PI/iip
            dcp = cos(arg)
            dsp = sin(arg)

            m1d = (m - 1) * im1 + 1
            m2s = 1 - im2

            idp2 = ido + 2
            nbd = (ido-1)/2
            iipp2 = iip+2
            iipph = (iip+1)/2

            if (ido >= l1) then
                do k=1,l1
                    do i=1,ido
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i,k,1) = cc(m1,i,1,k)
                        end do
                    end do
                end do
            else
                do i=1,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i,k,1) = cc(m1,i,1,k)
                        end do
                    end do
                end do
            end if

            do j=2,iipph
                jc = iipp2-j
                j2 = j+j
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,1,k,j) = cc(m1,ido,j2-2,k)+cc(m1,ido,j2-2,k)
                        ch(m2,1,k,jc) = cc(m1,1,j2-1,k)+cc(m1,1,j2-1,k)
                    end do
                end do
            end do

            if (ido /= 1) then
                if (nbd >= l1) then
                    do j=2,iipph
                        jc = iipp2-j
                        do k=1,l1
                            do i=3,ido,2
                                ic = idp2-i
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = cc(m1,i-1,2*j-1,k)+cc(m1,ic-1,2*j-2,k)
                                    ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k)-cc(m1,ic-1,2*j-2,k)
                                    ch(m2,i,k,j) = cc(m1,i,2*j-1,k)-cc(m1,ic,2*j-2,k)
                                    ch(m2,i,k,jc) = cc(m1,i,2*j-1,k)+cc(m1,ic,2*j-2,k)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2,iipph
                        jc = iipp2-j
                        do i=3,ido,2
                            ic = idp2-i
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = cc(m1,i-1,2*j-1,k)+cc(m1,ic-1,2*j-2,k)
                                    ch(m2,i-1,k,jc) = cc(m1,i-1,2*j-1,k)-cc(m1,ic-1,2*j-2,k)
                                    ch(m2,i,k,j) = cc(m1,i,2*j-1,k)-cc(m1,ic,2*j-2,k)
                                    ch(m2,i,k,jc) = cc(m1,i,2*j-1,k)+cc(m1,ic,2*j-2,k)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

            ar1 = ONE
            ai1 = ZERO
            do l=2,iipph
                lc = iipp2-l
                ar1h = dcp*ar1-dsp*ai1
                ai1 = dcp*ai1+dsp*ar1
                ar1 = ar1h
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        c2(m1,ik,l) = ch2(m2,ik,1)+ar1*ch2(m2,ik,2)
                        c2(m1,ik,lc) = ai1*ch2(m2,ik,iip)
                    end do
                end do
                dc2 = ar1
                ds2 = ai1
                ar2 = ar1
                ai2 = ai1
                do j=3,iipph
                    jc = iipp2-j
                    ar2h = dc2*ar2-ds2*ai2
                    ai2 = dc2*ai2+ds2*ar2
                    ar2 = ar2h
                    do ik=1,idl1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            c2(m1,ik,l) = c2(m1,ik,l)+ar2*ch2(m2,ik,j)
                            c2(m1,ik,lc) = c2(m1,ik,lc)+ai2*ch2(m2,ik,jc)
                        end do
                    end do
                end do
            end do

            do j=2,iipph
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch2(m2,ik,1) = ch2(m2,ik,1)+ch2(m2,ik,j)
                    end do
                end do
            end do

            do j=2,iipph
                jc = iipp2-j
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,1,k,j) = c1(m1,1,k,j)-c1(m1,1,k,jc)
                        ch(m2,1,k,jc) = c1(m1,1,k,j)+c1(m1,1,k,jc)
                    end do
                end do
            end do

            if (ido /= 1) then
                if (nbd >= l1) then
                    do j=2,iipph
                        jc = iipp2-j
                        do k=1,l1
                            do i=3,ido,2
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = c1(m1,i-1,k,j)-c1(m1,i,k,jc)
                                    ch(m2,i-1,k,jc) = c1(m1,i-1,k,j)+c1(m1,i,k,jc)
                                    ch(m2,i,k,j) = c1(m1,i,k,j)+c1(m1,i-1,k,jc)
                                    ch(m2,i,k,jc) = c1(m1,i,k,j)-c1(m1,i-1,k,jc)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2,iipph
                        jc = iipp2-j
                        do i=3,ido,2
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = c1(m1,i-1,k,j)-c1(m1,i,k,jc)
                                    ch(m2,i-1,k,jc) = c1(m1,i-1,k,j)+c1(m1,i,k,jc)
                                    ch(m2,i,k,j) = c1(m1,i,k,j)+c1(m1,i-1,k,jc)
                                    ch(m2,i,k,jc) = c1(m1,i,k,j)-c1(m1,i-1,k,jc)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

            if (ido /= 1) then
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        c2(m1,ik,1) = ch2(m2,ik,1)
                    end do
                end do
                do j=2,iip
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            c1(m1,1,k,j) = ch(m2,1,k,j)
                        end do
                    end do
                end do
                if (l1 >= nbd ) then
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        idij = is
                        do i=3,ido,2
                            idij = idij+2
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    c1(m1,i-1,k,j) = &
                                        wa(idij-1)*ch(m2,i-1,k,j) &
                                        - wa(idij)* ch(m2,i,k,j)
                                    c1(m1,i,k,j) = &
                                        wa(idij-1)*ch(m2,i,k,j) &
                                        + wa(idij)* ch(m2,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                else
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        do k=1,l1
                            idij = is
                            do i=3,ido,2
                                idij = idij+2
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    c1(m1,i-1,k,j) = &
                                        wa(idij-1)*ch(m2,i-1,k,j)&
                                        - wa(idij)*ch(m2,i,k,j)
                                    c1(m1,i,k,j) = &
                                        wa(idij-1)*ch(m2,i,k,j)&
                                        + wa(idij)*ch(m2,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

        end subroutine mradbg

    end subroutine mrftb1



    subroutine mrftf1(m,im,n,in,c,ch,wa,fac)

        integer (ip) in
        integer (ip) m
        integer (ip) n

        real (wp) c(in,*)
        real (wp) ch(m,*)
        real (wp) fac(15)
        integer (ip) i
        integer (ip) idl1
        integer (ip) ido
        integer (ip) im
        integer (ip) iip
        integer (ip) iw
        integer (ip) ix2
        integer (ip) ix3
        integer (ip) ix4
        integer (ip) j
        integer (ip) k1
        integer (ip) kh
        integer (ip) l1
        integer (ip) l2
        integer (ip) m2
        integer (ip) modn
        integer (ip) na
        integer (ip) nf
        integer (ip) nl
        real (wp) sn
        real (wp) tsn
        real (wp) tsnm
        real (wp) wa(n)

        nf = int(fac(2), kind=ip)
        na = 1
        l2 = n
        iw = n

        do k1=1,nf
            kh = nf-k1
            iip = int(fac(kh+3), kind=ip)
            l1 = l2/iip
            ido = n/l2
            idl1 = ido*l1
            iw = iw-(iip-1)*ido
            na = 1-na
            select case (iip)
                case (2)
                    if (na == 0) then
                        call mradf2(m,ido,l1,c,im,in,ch,1,m,wa(iw))
                    else
                        call mradf2(m,ido,l1,ch,1,m,c,im,in,wa(iw))
                    end if
                case(3)
                    ix2 = iw+ido
                    if (na == 0) then
                        call mradf3(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2))
                    else
                        call mradf3(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2))
                    end if
                case (4)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    if (na == 0) then
                        call mradf4(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3))
                    else
                        call mradf4(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3))
                    end if
                case(5)
                    ix2 = iw+ido
                    ix3 = ix2+ido
                    ix4 = ix3+ido
                    if (na == 0) then
                        call mradf5(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    else
                        call mradf5(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    end if
                case default
                    if (ido == 1) then
                        na = 1-na
                    end if
                    if (na == 0) then
                        call mradfg(m,ido,iip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
                        na = 1
                    else
                        call mradfg(m,ido,iip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
                        na = 0
                    end if
            end select
            l2 = l1
        end do

        sn = ONE/n
        tsn =  TWO/n
        tsnm = -tsn
        modn = mod(n,2)

        if (modn /= 0) then
            nl = n-1
        else
            nl = n-2
        end if

        if (na == 0) then
            m2 = 1-im

            do i=1,m
                m2 = m2+im
                c(m2,1) = sn*ch(i,1)
            end do

            do j=2,nl,2
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,j) = tsn*ch(i,j)
                    c(m2,j+1) = tsnm*ch(i,j+1)
                end do
            end do

            if (modn == 0) then
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,n) = sn*ch(i,n)
                end do
            end if
        else
            m2 = 1-im
            do i=1,m
                m2 = m2+im
                c(m2,1) = sn*c(m2,1)
            end do
            do j=2,nl,2
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,j) = tsn*c(m2,j)
                    c(m2,j+1) = tsnm*c(m2,j+1)
                end do
            end do

            if (modn == 0) then
                m2 = 1-im
                do i=1,m
                    m2 = m2+im
                    c(m2,n) = sn*c(m2,n)
                end do
            end if
        end if

    contains

        subroutine mradf2(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,l1,2)
            real (wp) ch(in2,ido,2,l1)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,1,k) = cc(m1,1,k,1)+cc(m1,1,k,2)
                    ch(m2,ido,2,k) = cc(m1,1,k,1)-cc(m1,1,k,2)
                end do
            end do

            if (ido < 2) then
                return
            else if (ido == 2) then
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,1,2,k) = -cc(m1,ido,k,2)
                        ch(m2,ido,1,k) = cc(m1,ido,k,1)
                    end do
                end do
            else
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i,1,k) = cc(m1,i,k,1)+(wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))
                            ch(m2,ic,2,k) = (wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))-cc(m1,i,k,1)
                            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+(wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))
                            ch(m2,ic-1,2,k) = cc(m1,i-1,k,1)-(wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))
                        end do
                    end do
                end do
                if (mod(ido,2) /= 1) then
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,1,2,k) = -cc(m1,ido,k,2)
                            ch(m2,ido,1,k) = cc(m1,ido,k,1)
                        end do
                    end do
                end if
            end if

        end subroutine mradf2

        subroutine mradf3(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,l1,3)
            real (wp) ch(in2,ido,3,l1)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp), parameter :: TWO_PI =TWO * acos(-ONE)
            real (wp), parameter :: ARG = TWO_PI/3
            real (wp), parameter :: TAUR = cos(ARG)
            real (wp), parameter :: TAUI = sin(ARG)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,1,k) = cc(m1,1,k,1)+(cc(m1,1,k,2)+cc(m1,1,k,3))
                    ch(m2,1,3,k) = TAUI*(cc(m1,1,k,3)-cc(m1,1,k,2))
                    ch(m2,ido,2,k) = cc(m1,1,k,1)+TAUR*(cc(m1,1,k,2)+cc(m1,1,k,3))
                end do
            end do

            if (ido /= 1) then
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3)))

                            ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3)))

                            ch(m2,i-1,3,k) = (cc(m1,i-1,k,1)+TAUR*((wa1(i-2)* &
                                cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
                                cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))))+(TAUI*((wa1(i-2)* &
                                cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
                                cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))

                            ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)+TAUR*((wa1(i-2)* &
                                cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa2(i-2)* &
                                cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))))-(TAUI*((wa1(i-2)* &
                                cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa2(i-2)* &
                                cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))))

                            ch(m2,i,3,k) = (cc(m1,i,k,1)+TAUR*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))))+(TAUI*((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2))))

                            ch(m2,ic,2,k) = (TAUI*((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2))))-(cc(m1,i,k,1)+TAUR*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))))
                        end do
                    end do
                end do
            end if

        end subroutine mradf3


        subroutine mradf4(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,l1,4)
            real (wp) ch(in2,ido,4,l1)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp) wa3(ido)
            real (wp), parameter :: HALF_SQRT2 = sqrt(TWO)/2

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,1,k) = (cc(m1,1,k,2)+cc(m1,1,k,4)) &
                        +(cc(m1,1,k,1)+cc(m1,1,k,3))
                    ch(m2,ido,4,k) = (cc(m1,1,k,1)+cc(m1,1,k,3)) &
                        -(cc(m1,1,k,2)+cc(m1,1,k,4))
                    ch(m2,ido,2,k) = cc(m1,1,k,1)-cc(m1,1,k,3)
                    ch(m2,1,3,k) = cc(m1,1,k,4)-cc(m1,1,k,2)
                end do
            end do

            if (ido < 2) then
                return
            else if (ido == 2) then
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(m2,ido,1,k) = (HALF_SQRT2*(cc(m1,ido,k,2)-cc(m1,ido,k,4)))+ &
                            cc(m1,ido,k,1)
                        ch(m2,ido,3,k) = cc(m1,ido,k,1)-(HALF_SQRT2*(cc(m1,ido,k,2)- &
                            cc(m1,ido,k,4)))
                        ch(m2,1,2,k) = (-HALF_SQRT2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))- &
                            cc(m1,ido,k,3)
                        ch(m2,1,4,k) = (-HALF_SQRT2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))+ &
                            cc(m1,ido,k,3)
                    end do
                end do
            else
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,i-1,1,k) = ((wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2))+(wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4)))+(cc(m1,i-1,k,1)+(wa2(i-2)*cc(m1,i-1,k,3)+ &
                                wa2(i-1)*cc(m1,i,k,3)))

                            ch(m2,ic-1,4,k) = (cc(m1,i-1,k,1)+(wa2(i-2)*cc(m1,i-1,k,3)+ &
                                wa2(i-1)*cc(m1,i,k,3)))-((wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))+(wa3(i-2)*cc(m1,i-1,k,4)+ &
                                wa3(i-1)*cc(m1,i,k,4)))

                            ch(m2,i,1,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))+(cc(m1,i,k,1)+(wa2(i-2)*cc(m1,i,k,3)- &
                                wa2(i-1)*cc(m1,i-1,k,3)))

                            ch(m2,ic,4,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))-(cc(m1,i,k,1)+(wa2(i-2)*cc(m1,i,k,3)- &
                                wa2(i-1)*cc(m1,i-1,k,3)))

                            ch(m2,i-1,3,k) = ((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))-(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))+(cc(m1,i-1,k,1)-(wa2(i-2)*cc(m1,i-1,k,3)+ &
                                wa2(i-1)*cc(m1,i,k,3)))

                            ch(m2,ic-1,2,k) = (cc(m1,i-1,k,1)-(wa2(i-2)*cc(m1,i-1,k,3)+ &
                                wa2(i-1)*cc(m1,i,k,3)))-((wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)* &
                                cc(m1,i-1,k,2))-(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))

                            ch(m2,i,3,k) = ((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))+(cc(m1,i,k,1)-(wa2(i-2)*cc(m1,i,k,3)- &
                                wa2(i-1)*cc(m1,i-1,k,3)))

                            ch(m2,ic,2,k) = ((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))-(cc(m1,i,k,1)-(wa2(i-2)*cc(m1,i,k,3)- &
                                wa2(i-1)*cc(m1,i-1,k,3)))
                        end do
                    end do
                end do
                if (mod(ido,2) /= 1) then
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,ido,1,k) = (HALF_SQRT2*(cc(m1,ido,k,2)-cc(m1,ido,k,4)))+ &
                                cc(m1,ido,k,1)
                            ch(m2,ido,3,k) = cc(m1,ido,k,1)-(HALF_SQRT2*(cc(m1,ido,k,2)- &
                                cc(m1,ido,k,4)))
                            ch(m2,1,2,k) = (-HALF_SQRT2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))- &
                                cc(m1,ido,k,3)
                            ch(m2,1,4,k) = (-HALF_SQRT2*(cc(m1,ido,k,2)+cc(m1,ido,k,4)))+ &
                                cc(m1,ido,k,3)
                        end do
                    end do
                end if
            end if

        end subroutine mradf4

        subroutine mradf5(m,ido,l1,cc,im1,in1,ch,im2,in2,wa1,wa2,wa3,wa4)

            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) l1

            real (wp) cc(in1,ido,l1,5)
            real (wp) ch(in2,ido,5,l1)
            integer (ip) i
            integer (ip) ic
            integer (ip) idp2
            integer (ip) im1
            integer (ip) im2
            integer (ip) k
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            real (wp) wa1(ido)
            real (wp) wa2(ido)
            real (wp) wa3(ido)
            real (wp) wa4(ido)
            real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
            real (wp), parameter :: ARG = TWO_PI/5
            real (wp), parameter :: TR11 = cos(ARG)
            real (wp), parameter :: TI11 = sin(ARG)
            real (wp), parameter :: TR12 = cos(TWO*ARG)
            real (wp), parameter :: TI12 = sin(TWO*ARG)

            m1d = (m-1)*im1+1
            m2s = 1-im2

            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(m2,1,1,k) = cc(m1,1,k,1)+(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
                        (cc(m1,1,k,4)+cc(m1,1,k,3))
                    ch(m2,ido,2,k) = cc(m1,1,k,1)+TR11*(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
                        TR12*(cc(m1,1,k,4)+cc(m1,1,k,3))
                    ch(m2,1,3,k) = TI11*(cc(m1,1,k,5)-cc(m1,1,k,2))+TI12* &
                        (cc(m1,1,k,4)-cc(m1,1,k,3))
                    ch(m2,ido,4,k) = cc(m1,1,k,1)+TR12*(cc(m1,1,k,5)+cc(m1,1,k,2))+ &
                        TR11*(cc(m1,1,k,4)+cc(m1,1,k,3))
                    ch(m2,1,5,k) = TI12*(cc(m1,1,k,5)-cc(m1,1,k,2))-TI11* &
                        (cc(m1,1,k,4)-cc(m1,1,k,3))
                end do
            end do

            if (ido /= 1) then
                idp2 = ido+2
                do k=1,l1
                    do i=3,ido,2
                        ic = idp2-i
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2

                            ch(m2,i-1,1,k) = cc(m1,i-1,k,1)+((wa1(i-2)*cc(m1,i-1,k,2)+ &
                                wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
                                cc(m1,i,k,5)))+((wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))+(wa3(i-2)*cc(m1,i-1,k,4)+ &
                                wa3(i-1)*cc(m1,i,k,4)))

                            ch(m2,i,1,k) = cc(m1,i,k,1)+((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4)))

                            ch(m2,i-1,3,k) = cc(m1,i-1,k,1)+TR11* &
                                ( wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2) &
                                +wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5))+TR12* &
                                ( wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3) &
                                +wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))+TI11* &
                                ( wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2) &
                                -(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))+TI12* &
                                ( wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3) &
                                -(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4)))

                            ch(m2,ic-1,2,k) = cc(m1,i-1,k,1)+TR11* &
                                ( wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2) &
                                +wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5))+TR12* &
                                ( wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3) &
                                +wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))-(TI11* &
                                ( wa1(i-2)*cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2) &
                                -(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))+TI12* &
                                ( wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3) &
                                -(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

                            ch(m2,i,3,k) = (cc(m1,i,k,1)+TR11*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+TR12*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4))))+(TI11*((wa4(i-2)*cc(m1,i-1,k,5)+ &
                                wa4(i-1)*cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))+TI12*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))))

                            ch(m2,ic,2,k) = (TI11*((wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
                                cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))+TI12*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))))-(cc(m1,i,k,1)+TR11*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+TR12*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4))))

                            ch(m2,i-1,5,k) = (cc(m1,i-1,k,1)+TR12*((wa1(i-2)* &
                                cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)* &
                                cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5)))+TR11*((wa2(i-2)* &
                                cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))+(wa3(i-2)* &
                                cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))))+(TI12*((wa1(i-2)* &
                                cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa4(i-2)* &
                                cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))-TI11*((wa2(i-2)* &
                                cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))-(wa3(i-2)* &
                                cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

                            ch(m2,ic-1,4,k) = (cc(m1,i-1,k,1)+TR12*((wa1(i-2)* &
                                cc(m1,i-1,k,2)+wa1(i-1)*cc(m1,i,k,2))+(wa4(i-2)* &
                                cc(m1,i-1,k,5)+wa4(i-1)*cc(m1,i,k,5)))+TR11*((wa2(i-2)* &
                                cc(m1,i-1,k,3)+wa2(i-1)*cc(m1,i,k,3))+(wa3(i-2)* &
                                cc(m1,i-1,k,4)+wa3(i-1)*cc(m1,i,k,4))))-(TI12*((wa1(i-2)* &
                                cc(m1,i,k,2)-wa1(i-1)*cc(m1,i-1,k,2))-(wa4(i-2)* &
                                cc(m1,i,k,5)-wa4(i-1)*cc(m1,i-1,k,5)))-TI11*((wa2(i-2)* &
                                cc(m1,i,k,3)-wa2(i-1)*cc(m1,i-1,k,3))-(wa3(i-2)* &
                                cc(m1,i,k,4)-wa3(i-1)*cc(m1,i-1,k,4))))

                            ch(m2,i,5,k) = (cc(m1,i,k,1)+TR12*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+TR11*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4))))+(TI12*((wa4(i-2)*cc(m1,i-1,k,5)+ &
                                wa4(i-1)*cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))-TI11*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))))

                            ch(m2,ic,4,k) = (TI12*((wa4(i-2)*cc(m1,i-1,k,5)+wa4(i-1)* &
                                cc(m1,i,k,5))-(wa1(i-2)*cc(m1,i-1,k,2)+wa1(i-1)* &
                                cc(m1,i,k,2)))-TI11*((wa3(i-2)*cc(m1,i-1,k,4)+wa3(i-1)* &
                                cc(m1,i,k,4))-(wa2(i-2)*cc(m1,i-1,k,3)+wa2(i-1)* &
                                cc(m1,i,k,3))))-(cc(m1,i,k,1)+TR12*((wa1(i-2)*cc(m1,i,k,2)- &
                                wa1(i-1)*cc(m1,i-1,k,2))+(wa4(i-2)*cc(m1,i,k,5)-wa4(i-1)* &
                                cc(m1,i-1,k,5)))+TR11*((wa2(i-2)*cc(m1,i,k,3)-wa2(i-1)* &
                                cc(m1,i-1,k,3))+(wa3(i-2)*cc(m1,i,k,4)-wa3(i-1)* &
                                cc(m1,i-1,k,4))))
                        end do
                    end do
                end do
            end if

        end subroutine mradf5

        subroutine mradfg(m,ido,iip,l1,idl1,cc,c1,c2,im1,in1,ch,ch2,im2,in2,wa)

            integer (ip) idl1
            integer (ip) ido
            integer (ip) in1
            integer (ip) in2
            integer (ip) iip
            integer (ip) l1

            real (wp) ai1
            real (wp) ai2
            real (wp) ar1
            real (wp) ar1h
            real (wp) ar2
            real (wp) ar2h
            real (wp) arg
            real (wp) c1(in1,ido,l1,iip)
            real (wp) c2(in1,idl1,iip)
            real (wp) cc(in1,ido,iip,l1)
            real (wp) ch(in2,ido,l1,iip)
            real (wp) ch2(in2,idl1,iip)
            real (wp) dc2
            real (wp) dcp
            real (wp) ds2
            real (wp) dsp
            integer (ip) i
            integer (ip) ic
            integer (ip) idij
            integer (ip) idp2
            integer (ip) ik
            integer (ip) im1
            integer (ip) im2
            integer (ip) iipp2
            integer (ip) iipph
            integer (ip) is
            integer (ip) j
            integer (ip) j2
            integer (ip) jc
            integer (ip) k
            integer (ip) l
            integer (ip) lc
            integer (ip) m
            integer (ip) m1
            integer (ip) m1d
            integer (ip) m2
            integer (ip) m2s
            integer (ip) nbd
            real (wp), parameter :: TWO_PI= TWO * acos(-ONE)
            real (wp) wa(ido)

            m1d = (m-1)*im1+1
            m2s = 1-im2
            arg = TWO_PI / iip
            dcp = cos(arg)
            dsp = sin(arg)
            iipph = (iip+1)/2
            iipp2 = iip+2
            idp2 = ido+2
            nbd = (ido-1)/2

            if (ido /= 1) then

                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch2(m2,ik,1) = c2(m1,ik,1)
                    end do
                end do

                do j=2,iip
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch(m2,1,k,j) = c1(m1,1,k,j)
                        end do
                    end do
                end do

                if ( l1 >= nbd ) then
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        idij = is
                        do i=3,ido,2
                            idij = idij+2
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = wa(idij-1)*c1(m1,i-1,k,j)+wa(idij)*c1(m1,i,k,j)
                                    ch(m2,i,k,j) = wa(idij-1)*c1(m1,i,k,j)-wa(idij)*c1(m1,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                else
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        do k=1,l1
                            idij = is
                            do i=3,ido,2
                                idij = idij+2
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    ch(m2,i-1,k,j) = wa(idij-1)*c1(m1,i-1,k,j)+wa(idij)*c1(m1,i,k,j)
                                    ch(m2,i,k,j) = wa(idij-1)*c1(m1,i,k,j)-wa(idij)*c1(m1,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                end if

                if (nbd >= l1) then
                    do j=2,iipph
                        jc = iipp2-j
                        do k=1,l1
                            do i=3,ido,2
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    c1(m1,i-1,k,j) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
                                    c1(m1,i-1,k,jc) = ch(m2,i,k,j)-ch(m2,i,k,jc)
                                    c1(m1,i,k,j) = ch(m2,i,k,j)+ch(m2,i,k,jc)
                                    c1(m1,i,k,jc) = ch(m2,i-1,k,jc)-ch(m2,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2,iipph
                        jc = iipp2-j
                        do i=3,ido,2
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    c1(m1,i-1,k,j) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
                                    c1(m1,i-1,k,jc) = ch(m2,i,k,j)-ch(m2,i,k,jc)
                                    c1(m1,i,k,j) = ch(m2,i,k,j)+ch(m2,i,k,jc)
                                    c1(m1,i,k,jc) = ch(m2,i-1,k,jc)-ch(m2,i-1,k,j)
                                end do
                            end do
                        end do
                    end do
                end if
            else
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        c2(m1,ik,1) = ch2(m2,ik,1)
                    end do
                end do
            end if
            do j=2,iipph
                jc = iipp2-j
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        c1(m1,1,k,j) = ch(m2,1,k,j)+ch(m2,1,k,jc)
                        c1(m1,1,k,jc) = ch(m2,1,k,jc)-ch(m2,1,k,j)
                    end do
                end do
            end do

            ar1 = ONE
            ai1 = ZERO
            do l=2,iipph
                lc = iipp2-l
                ar1h = dcp*ar1-dsp*ai1
                ai1 = dcp*ai1+dsp*ar1
                ar1 = ar1h
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch2(m2,ik,l) = c2(m1,ik,1)+ar1*c2(m1,ik,2)
                        ch2(m2,ik,lc) = ai1*c2(m1,ik,iip)
                    end do
                end do
                dc2 = ar1
                ds2 = ai1
                ar2 = ar1
                ai2 = ai1
                do j=3,iipph
                    jc = iipp2-j
                    ar2h = dc2*ar2-ds2*ai2
                    ai2 = dc2*ai2+ds2*ar2
                    ar2 = ar2h
                    do ik=1,idl1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            ch2(m2,ik,l) = ch2(m2,ik,l)+ar2*c2(m1,ik,j)
                            ch2(m2,ik,lc) = ch2(m2,ik,lc)+ai2*c2(m1,ik,jc)
                        end do
                    end do
                end do
            end do
            do j=2,iipph
                do ik=1,idl1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch2(m2,ik,1) = ch2(m2,ik,1)+c2(m1,ik,j)
                    end do
                end do
            end do

            if (ido >= l1) then
                do k=1,l1
                    do i=1,ido
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc(m1,i,1,k) = ch(m2,i,k,1)
                        end do
                    end do
                end do
            else
                do i=1,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc(m1,i,1,k) = ch(m2,i,k,1)
                        end do
                    end do
                end do
            end if

            do j=2,iipph
                jc = iipp2-j
                j2 = j+j
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        cc(m1,ido,j2-2,k) = ch(m2,1,k,j)
                        cc(m1,1,j2-1,k) = ch(m2,1,k,jc)
                    end do
                end do
            end do
            if (ido /= 1) then
                if (nbd >= l1) then
                    do j=2,iipph
                        jc = iipp2-j
                        j2 = j+j
                        do k=1,l1
                            do i=3,ido,2
                                ic = idp2-i
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    cc(m1,i-1,j2-1,k) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
                                    cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j)-ch(m2,i-1,k,jc)
                                    cc(m1,i,j2-1,k) = ch(m2,i,k,j)+ch(m2,i,k,jc)
                                    cc(m1,ic,j2-2,k) = ch(m2,i,k,jc)-ch(m2,i,k,j)
                                end do
                            end do
                        end do
                    end do
                else
                    do j=2,iipph
                        jc = iipp2-j
                        j2 = j+j
                        do i=3,ido,2
                            ic = idp2-i
                            do k=1,l1
                                m2 = m2s
                                do m1=1,m1d,im1
                                    m2 = m2+im2
                                    cc(m1,i-1,j2-1,k) = ch(m2,i-1,k,j)+ch(m2,i-1,k,jc)
                                    cc(m1,ic-1,j2-2,k) = ch(m2,i-1,k,j)-ch(m2,i-1,k,jc)
                                    cc(m1,i,j2-1,k) = ch(m2,i,k,j)+ch(m2,i,k,jc)
                                    cc(m1,ic,j2-2,k) = ch(m2,i,k,jc)-ch(m2,i,k,j)
                                end do
                            end do
                        end do
                    end do
                end if
            end if

        end subroutine mradfg

    end subroutine mrftf1



    subroutine rfft2b ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier)
        !
        ! rfft2b: 64-bit float precision backward fast Fourier transform, 2D.
        !
        ! Purpose:
        !
        !  Computes the two-dimensional discrete Fourier transform of the
        !  complex Fourier coefficients a real periodic array.  This transform is
        !  known as the backward transform or Fourier synthesis, transforming from
        !  spectral to physical space.  Routine RFFT2B is normalized: a call to
        !  RFFT2B followed by a call to RFFT2F (or vice-versa) reproduces the
        !  original array within roundoff error.
        !
        !  Parameters:
        !
        !  integer LDIM, the first dimension of the 2D real
        !  array R, which must be at least 2*(L/2+1).
        !
        !  integer L, the number of elements to be transformed
        !  in the first dimension of the two-dimensional real array R.  The value of
        !  L must be less than or equal to that of LDIM.  The transform is most
        !  efficient when L is a product of small primes.
        !
        !  integer M, the number of elements to be transformed
        !  in the second dimension of the two-dimensional real array R.  The transform
        !  is most efficient when M is a product of small primes.
        !
        !  Input/real R(LDIM,M), the real array of two
        !  dimensions.  On input, R contains the L/2+1-by-M complex subarray of
        !  spectral coefficients, on output, the physical coefficients.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFT2I before the first call to routine RFFT2F
        !  or RFFT2B with lengths L and M.  wsave's contents may be re-used for
        !  subsequent calls to RFFT2F and RFFT2B with the same transform lengths
        !  L and M.
        !
        !  integer LENSAV, the number of elements in the wsave
        !  array.  LENSAV must be at least L + M + INT(LOG(REAL(L)))
        !  + INT(LOG(REAL(M))) + 8.
        !
        !  Workspace, real (wp) WORK(LENWRK).  WORK provides workspace, and
        !  its contents need not be saved between calls to routines RFFT2B and RFFT2F.
        !
        !  integer  LENWRK, the number of elements in the WORK
        !  array.  LENWRK must be at least LDIM*M.
        !
        !  integer IER, the error_flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  6, input parameter LDIM < 2*(L/2+1);
        !  20, input error returned by lower level routine.
        !
        integer (ip) ldim
        integer (ip) lensav
        integer (ip) lenwrk
        integer (ip) m


        integer (ip) ier
        integer (ip) local_error_flag

        integer (ip) l
        integer (ip) ldh
        integer (ip) ldw
        integer (ip) ldx
        integer (ip) lwsav
        integer (ip) mmsav
        integer (ip) modl
        integer (ip) modm
        integer (ip) mwsav
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) r(ldim,m)

        ier = 0
        !
        !==> verify lensav
        !
        lwsav = l+int(log(real(l, kind=wp))/log(TWO))+4
        mwsav = get_1d_saved_workspace_length(m)
        mmsav = m+int(log(real(m, kind=wp))/log(TWO))+4
        modl = mod(l,2)
        modm = mod(m,2)

        if (lensav < lwsav+mwsav+mmsav) then
            ier = 2
            call fft_error_handler('rfft2f', 6)
            return
        end if
        !
        ! verify lenwrk
        !
        if (lenwrk < (l+1)*m) then
            ier = 3
            call fft_error_handler('rfft2f', 8)
            return
        end if
        !
        ! verify ldim is as big as l
        !
        if (ldim < l) then
            ier = 5
            call fft_error_handler('rfft2f', -6)
            return
        end if
        !
        !==> Transform second dimension of array
        !
        associate( mj => 2*((m+1)/2)-1 )

            r(1,2:mj) = TWO * r(1,2:mj)

        end associate

        r(1,3:m:2) = -r(1,3:m:2)

        call rfftmb(1,1,m,ldim,r,m*ldim, wsave(lwsav+mwsav+1),mmsav,work,lenwrk,local_error_flag)

        ldh = int((l+1)/2)

        if ( 1 < ldh ) then
            ldw = ldh+ldh
            !
            !  r and work are switched because the the first dimension
            !  of the input to complex cfftmf must be even.
            !
            call r2w(ldim,ldw,l,m,r,work)

            call cfftmb(ldh-1,1,m,ldh,work(2),ldh*m, &
                wsave(lwsav+1),mwsav,r,l*m, local_error_flag)

            if (local_error_flag/=0) then
                ier=20
                call fft_error_handler('rfft2b',-5)
                return
            end if

            call w2r(ldim,ldw,l,m,r,work)
        end if

        if (modl == 0) then

            associate( mj => 2*((m+1)/2)-1 )

                r(l,2:mj) = TWO * r(l,2:mj)

            end associate

            r(l,3:m:2) = -r(l,3:m:2)

            call rfftmb(1,1,m,ldim,r(l,1),m*ldim, &
                wsave(lwsav+mwsav+1),mmsav,work,lenwrk,local_error_flag)
        end if
        !
        !==> Transform first dimension of array
        !
        ldx = 2*int((l+1)/2)-1

        r(2:ldx,1:m) = TWO * r(2:ldx,1:m)
        r(3:ldx:2,1:m) = -r(3:ldx:2,1:m)

        associate( &
            arg_1 => m*ldim, &
            arg_2 => l+int(log(real(l, kind=wp) )/log(TWO))+4 &
            )

            call rfftmb(m,ldim,l,1,r,arg_1,wsave(1), arg_2,work,lenwrk,local_error_flag)

        end associate

        if (local_error_flag /= 0) then
            ier=20
            call fft_error_handler('rfft2f',-5)
        end if

    end subroutine rfft2b


    subroutine rfft2f(ldim, l, m, r, wsave, lensav, work, lenwrk, ier)
        !
        ! rfft2f: 64-bit float precision forward fast Fourier transform, 2D.
        !
        ! Purpose:
        !
        !  RFFT2F computes the two-dimensional discrete Fourier transform of a
        !  real periodic array.  This transform is known as the forward transform
        !  or Fourier analysis, transforming from physical to spectral space.
        !  Routine rfft2f is normalized: a call to rfft2f followed by a call to
        !  rfft2b (or vice-versa) reproduces the original array within roundoff
        !  error.
        !
        !  Parameters:
        !
        !  integer LDIM, the first dimension of the 2D real
        !  array R, which must be at least 2*(L/2+1).
        !
        !  integer L, the number of elements to be transformed
        !  in the first dimension of the two-dimensional real array R.  The value
        !  of L must be less than or equal to that of LDIM.  The transform is most
        !  efficient when L is a product of small primes.
        !
        !  integer M, the number of elements to be transformed
        !  in the second dimension of the two-dimensional real array R.  The
        !  transform is most efficient when M is a product of small primes.
        !
        !  Input/real R(LDIM,M), the real array of two
        !  dimensions.  On input, containing the L-by-M physical data to be
        !  transformed.  On output, the spectral coefficients.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFT2I before the first call to routine RFFT2F
        !  or RFFT2B with lengths L and M.  wsave's contents may be re-used for
        !  subsequent calls to RFFT2F and RFFT2B with the same transform lengths.
        !
        !  integer LENSAV, the number of elements in the wsave
        !  array.  LENSAV must be at least L + M + INT(LOG(REAL(L)))
        !  + INT(LOG(REAL(M))) + 8.
        !
        !  Workspace, real (wp) WORK(LENWRK), provides workspace, and its
        !  contents need not be saved between calls to routines RFFT2F and RFFT2B.
        !
        !  integer LENWRK, the number of elements in the WORK
        !  array.  LENWRK must be at least LDIM*M.
        !
        !  integer IER, the error_flag.
        !  0, successful exit;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  6, input parameter LDIM < 2*(L+1);
        !  20, input error returned by lower level routine.
        !


        integer (ip) ldim
        integer (ip) lensav
        integer (ip) lenwrk
        integer (ip) m


        integer (ip) ier
        integer (ip) local_error_flag

        integer (ip) l
        integer (ip) ldh
        integer (ip) ldw
        integer (ip) ldx
        integer (ip) lwsav
        integer (ip) mmsav
        integer (ip) modl
        integer (ip) modm
        integer (ip) mwsav
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)
        real (wp) r(ldim,m)

        ier = 0
        !
        !==> verify lensav
        !
        lwsav = l+int(log(real (l, kind=wp))/log(TWO))+4
        mwsav = get_1d_saved_workspace_length(m)
        mmsav = m+int(log(real(m, kind=wp))/log(TWO))+4

        if (lensav < lwsav+mwsav+mmsav) then
            ier = 2
            call fft_error_handler('rfft2f', 6)
            return
        end if
        !
        !==>  verify lenwrk
        !
        if (lenwrk < (l+1)*m) then
            ier = 3
            call fft_error_handler('rfft2f', 8)
            return
        end if
        !
        !==>  verify ldim is as big as l
        !
        if (ldim < l) then
            ier = 5
            call fft_error_handler('rfft2f', -6)
            return
        end if
        !
        !==>  Transform first dimension of array
        !
        associate( &
            arg_1 => m*ldim, &
            arg_2 => l+int(log(real(l, kind=wp))/log(TWO))+4 &
            )

            call rfftmf(m,ldim,l,1,r,arg_1,wsave(1), arg_2,work,size(work),local_error_flag)

        end associate

        if (local_error_flag /= 0) then
            ier=20
            call fft_error_handler('rfft2f',-5)
            return
        end if

        ldx = 2*int((l+1)/2)-1

        r(2:ldx,1:m) = HALF * r(2:ldx,1:m)

        r(3:ldx:2,1:m) = -r(3:ldx:2,1:m)

        !==>  Reshuffle to add in nyquist imaginary components
        !
        modl = mod(l,2)
        modm = mod(m,2)
        !
        !==>  Transform second dimension of array
        !
        call rfftmf(1,1,m,ldim,r,m*ldim, &
            wsave(lwsav+mwsav+1),mmsav,work,size(work),local_error_flag)

        associate( mj => 2*((m+1)/2)-1 )

            r(1,2:mj) = HALF * r(1,2:mj)

        end associate

        r(1,3:m:2) = -r(1,3:m:2)

        ldh = int((l+1)/2)

        if ( 1 < ldh ) then
            ldw = 2*ldh
            !
            !==> r and work are switched because the the first dimension
            !    of the input to complex cfftmf must be even.
            !
            call r2w(ldim,ldw,l,m,r,work)
            call cfftmf(ldh-1,1,m,ldh,work(2),ldh*m, &
                wsave(lwsav+1),mwsav,r,l*m, local_error_flag)

            if (local_error_flag /= 0) then
                ier=20
                call fft_error_handler('rfft2f',-5)
                return
            end if

            call w2r(ldim,ldw,l,m,r,work)
        end if

        if (modl == 0) then

            call rfftmf(1,1,m,ldim,r(l,1),m*ldim, &
                wsave(lwsav+mwsav+1),mmsav,work,size(work),local_error_flag)

            associate( mj => 2*((m+1)/2)-1 )

                r(l,2:mj) = HALF * r(l,2:mj)

            end associate

            r(l,3:m:2) = -r(l,3:m:2)

        end if

        if (local_error_flag /= 0) then
            ier=20
            call fft_error_handler('rfft2f',-5)
        end if

    end subroutine rfft2f


    subroutine rfftmb(lot, jump, n, inc, r, lenr, wsave, &
        lensav, work, lenwrk, ier)
        !
        ! rfftmb: 64-bit float precision backward fft, 1d, multiple vectors.
        !
        ! Purpose:
        !
        !  Computes the one-dimensional Fourier transform of multiple
        !  periodic sequences within a real array.  This transform is referred
        !  to as the backward transform or Fourier synthesis, transforming the
        !  sequences from spectral to physical space.
        !
        !  This transform is normalized since a call to RFFTMB followed
        !  by a call to RFFTMF (or vice-versa) reproduces the original
        !  array  within roundoff error.
        !
        !  Parameters:
        !
        !  integer LOT, the number of sequences to be transformed
        !  within array R.
        !
        !  integer JUMP, the increment between the locations, in
        !  array R, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations, in
        !  array R, of two consecutive elements within the same sequence.
        !
        !  Input/real R(LENR), real array containing LOT
        !  sequences, each having length N.  R can have any number of dimensions,
        !  but the total number of locations must be at least LENR.  On input, the
        !  spectral data to be transformed, on output the physical data.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least (LOT-1)*JUMP + INC*(N-1) + 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFTMI before the first call to routine RFFTMF
        !  or RFFTMB for a given transform length N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must  be at least N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least LOT*N.
        !
        !  integer IER, error_flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  4, input parameters INC, JUMP, N, LOT are not consistent.
        !
        integer (ip) lenr
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) inc
        integer (ip) jump
        integer (ip) lot
        integer (ip) n
        real (wp) r(lenr)
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)

        !
        !==> Validity of calling arguments
        !
        if (lenr < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('rfftmb ', 6)
            return
        else if (lensav < n+int(log(real(n, kind=wp))/log(TWO))+4) then
            ier = 2
            call fft_error_handler('rfftmb ', 8)
            return
        else if (lenwrk < lot*n) then
            ier = 3
            call fft_error_handler('rfftmb ', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('rfftmb ', -1)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            call mrftb1(lot,jump,n,inc,r,work,wsave,wsave(n+1))
        end if

    end subroutine rfftmb


    subroutine rfftmf(lot, jump, n, inc, r, lenr, wsave, lensav, &
        work, lenwrk, ier)
        !
        ! RFFTMF: 64-bit float precision forward FFT, 1D, multiple vectors.
        !
        !  Purpose:
        !
        !  RFFTMF computes the one-dimensional Fourier transform of multiple
        !  periodic sequences within a real array.  This transform is referred
        !  to as the forward transform or Fourier analysis, transforming the
        !  sequences from physical to spectral space.
        !
        !  This transform is normalized since a call to RFFTMF followed
        !  by a call to RFFTMB (or vice-versa) reproduces the original array
        !  within roundoff error.
        !
        !  Parameters:
        !
        !  integer LOT, the number of sequences to be transformed
        !  within array R.
        !
        !  integer JUMP, the increment between the locations, in
        !  array R, of the first elements of two consecutive sequences to be
        !  transformed.
        !
        !  integer N, the length of each sequence to be
        !  transformed.  The transform is most efficient when N is a product of
        !  small primes.
        !
        !  integer INC, the increment between the locations,
        !  in array R, of two consecutive elements within the same sequence.
        !
        !  Input/real R(LENR), real array containing LOT
        !  sequences, each having length N.  R can have any number of dimensions, but
        !  the total number of locations must be at least LENR.  On input, the
        !  physical data to be transformed, on output the spectral data.
        !
        !  integer LENR, the dimension of the R array.
        !  LENR must be at least (LOT-1)*JUMP + INC*(N-1) + 1.
        !
        !  Input, real (wp) wsave(LENSAV).  wsave's contents must be
        !  initialized with a call to RFFTMI before the first call to routine RFFTMF
        !  or RFFTMB for a given transform length N.
        !
        !  integer LENSAV, the dimension of the wsave array.
        !  LENSAV must be at least N + INT(LOG(REAL(N))) + 4.
        !
        !  Workspace, real (wp) WORK(LENWRK).
        !
        !  integer LENWRK, the dimension of the WORK array.
        !  LENWRK must be at least LOT*N.
        !
        !  integer IER, error_flag.
        !  0, successful exit;
        !  1, input parameter LENR not big enough;
        !  2, input parameter LENSAV not big enough;
        !  3, input parameter LENWRK not big enough;
        !  4, input parameters INC, JUMP, N, LOT are not consistent.
        !


        integer (ip) lenr
        integer (ip) lensav
        integer (ip) lenwrk

        integer (ip) ier
        integer (ip) inc
        integer (ip) jump
        integer (ip) lot
        integer (ip) n
        real (wp) r(lenr)
        real (wp) work(lenwrk)
        real (wp) wsave(lensav)

        !
        !==> Check validity of input arguments
        !
        if (lenr < (lot-1)*jump + inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('rfftmf ', 6)
            return
        else if (lensav < n + int(log(real(n, kind=wp) )/log(TWO)) +4) then
            ier = 2
            call fft_error_handler('rfftmf ', 8)
            return
        else if (lenwrk < lot*n) then
            ier = 3
            call fft_error_handler('rfftmf ', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ier = 4
            call fft_error_handler('rfftmf ', -1)
            return
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            call mrftf1(lot,jump,n,inc,r,work,wsave,wsave(n+1))
        end if

    end subroutine rfftmf




end module real_transform_routines
