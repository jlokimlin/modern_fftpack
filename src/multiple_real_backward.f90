submodule (real_transform_routines) multiple_real_backward

contains

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

        !
        !==> Validity of calling arguments
        !
        if (lenr < (lot-1)*jump + inc*(n-1) + 1) then
            ierror = 1
            call fft_error_handler('rfftmb ', 6)
            return
        else if (lensav < n+int(log(real(n, kind=wp))/log(TWO))+4) then
            ierror = 2
            call fft_error_handler('rfftmb ', 8)
            return
        else if (lenwrk < lot*n) then
            ierror = 3
            call fft_error_handler('rfftmb ', 10)
            return
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ierror = 4
            call fft_error_handler('rfftmb ', -1)
            return
        else
            ierror = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) call mrftb1(lot, jump, n, inc, r, work, wsave, wsave(n+1))

    end subroutine rfftmb


    subroutine mrftb1(m,im,n,in,c,ch,wa,fac)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
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
                        call mradb5(m,ido,l1,c,im,in,ch,1,m,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    else
                        call mradb5(m,ido,l1,ch,1,m,c,im,in,wa(iw),wa(ix2),wa(ix3),wa(ix4))
                    end if
                    na = 1-na
                case default
                    if (na == 0) then
                        call mradbg(m,ido,iip,l1,idl1,c,c,c,im,in,ch,ch,1,m,wa(iw))
                    else
                        call mradbg(m,ido,iip,l1,idl1,ch,ch,ch,1,m,c,c,im,in,wa(iw))
                    end if
                    if (ido == 1) then
                        na = 1-na
                    end if
            end select
            l1 = l2
            iw = iw+(iip-1)*ido
        end do

    end subroutine mrftb1



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
            if (l1 <= nbd) then
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
            if (l1 <= nbd) then
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



end submodule multiple_real_backward
