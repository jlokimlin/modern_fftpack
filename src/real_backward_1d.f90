submodule (real_transform_routines) real_backward_1d

contains

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

        !
        !==> Check validity of calling arguments
        !
        if (lenr < inc*(n-1) + 1) then
            ierror = 1
            call fft_error_handler('rfft1b ', 6)
        else if (lensav < &
            n + int(log(real(n, kind=wp) )/log(TWO))+4) then
            ierror = 2
            call fft_error_handler('rfft1b ', 8)
        else if (lenwrk < n) then
            ierror = 3
            call fft_error_handler('rfft1b ', 10)
        else
            ierror = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) call rfftb1(n,inc,r,work,wsave,wsave(n+1))

    end subroutine rfft1b



    subroutine rfftb1(n, in, c, ch, wa, fac)
        !-------------------------------------------------------------
        ! Dummy arguments
        !-------------------------------------------------------------
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: in
        real (wp),    intent (in out) :: c(in,*)
        real (wp),    intent (out)    :: ch(*)
        real (wp),    intent (in)     :: wa(n)
        real (wp),    intent (in)     :: fac(15)
        !-------------------------------------------------------------
        ! Local variables
        !-------------------------------------------------------------
        integer (ip) :: idl1, ido, iip
        integer (ip) :: workspace_indices(4)
        integer (ip) :: j, k, l1, l2
        integer (ip) :: modn, na, nf, nl
        !-------------------------------------------------------------

        nf = int(fac(2), kind=ip)
        na = 0

        factorization_loop: do k=1,nf
            iip = int(fac(k+2), kind=ip)
            na = 1-na
            if (iip <= 5) then
                cycle factorization_loop
            end if
            if (k == nf) then
                cycle factorization_loop
            end if
            na = 1-na
        end do factorization_loop

        modn = mod(n,2)

        if (modn /= 0) then
            nl = n-1
        else
            nl = n-2
        end if

        if (na /= 0) then
            ch(1) = c(1,1)
            ch(n) = c(1,n)
            do j=2,nl,2
                ch(j) = HALF*c(1,j)
                ch(j+1) = -HALF*c(1,j+1)
            end do
        else
            do j=2,nl,2
                c(1,j) = HALF*c(1,j)
                c(1,j+1) = -HALF*c(1,j+1)
            end do
        end if

        associate( &
            iw1 => workspace_indices(1), &
            iw2 => workspace_indices(2), &
            iw3 => workspace_indices(3), &
            iw4 => workspace_indices(4) &
            )

            l1 = 1
            iw1 = 1
            do k=1,nf
                iip = int(fac(k+2), kind=ip)
                l2 = iip*l1
                ido = n/l2
                idl1 = ido*l1
                select case (iip)
                    case (2)
                        select case (na)
                            case (0)
                                call r1f2kb(ido,l1,c,in,ch,1,wa(iw1))
                            case default
                                call r1f2kb(ido,l1,ch,1,c,in,wa(iw1))
                        end select
                        na = 1-na
                    case (3)
                        iw2 = iw1+ido
                        select case (na)
                            case (0)
                                call r1f3kb(ido,l1,c,in,ch,1,wa(iw1),wa(iw2))
                            case default
                                call r1f3kb(ido,l1,ch,1,c,in,wa(iw1),wa(iw2))
                        end select
                        na = 1-na
                    case (4)
                        iw2 = iw1+ido
                        iw3 = iw2+ido
                        select case (na)
                            case (0)
                                call r1f4kb(ido,l1,c,in,ch,1,wa(iw1),wa(iw2),wa(iw3))
                            case default
                                call r1f4kb(ido,l1,ch,1,c,in,wa(iw1),wa(iw2),wa(iw3))
                        end select
                        na = 1-na
                    case (5)
                        iw2 = iw1+ido
                        iw3 = iw2+ido
                        iw4 = iw3+ido
                        select case (na)
                            case (0)
                                call r1f5kb(ido,l1,c,in,ch,1,wa(iw1),wa(iw2),wa(iw3),wa(iw4))
                            case default
                                call r1f5kb(ido,l1,ch,1,c,in,wa(iw1),wa(iw2),wa(iw3),wa(iw4))
                        end select
                        na = 1-na
                    case default
                        select case (na)
                            case (0)
                                call r1fgkb(ido,iip,l1,idl1,c,c,c,in,ch,ch,1,wa(iw1))
                            case default
                                call r1fgkb(ido,iip,l1,idl1,ch,ch,ch,1,c,c,in,wa(iw1))
                        end select
                        if (ido == 1) then
                            na = 1-na
                        end if
                end select
                l1 = l2
                iw1 = iw1+(iip-1)*ido
            end do
        end associate

    end subroutine rfftb1


    subroutine r1f2kb(ido,l1,cc,in1,ch,in2,wa1)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        real (wp),    intent (in out) :: cc(in1,ido,2,l1)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,l1,2)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa1(ido)
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip) :: i, ic, idp2
        !--------------------------------------------------------------

        ch(1,1,:,1) = cc(1,1,1,:)+cc(1,ido,2,:)
        ch(1,1,:,2) = cc(1,1,1,:)-cc(1,ido,2,:)

        if (ido < 2) then
            return
        else
            select case (ido)
                case (2)
                    ch(1,ido,:,1) = cc(1,ido,1,:)+cc(1,ido,1,:)
                    ch(1,ido,:,2) = -(cc(1,1,2,:)+cc(1,1,2,:))
                case default
                    idp2 = ido+2
                    do i=3,ido,2
                        ic = idp2-i

                        ch(1,i-1,:,1) = cc(1,i-1,1,:)+cc(1,ic-1,2,:)
                        ch(1,i,:,1) = cc(1,i,1,:)-cc(1,ic,2,:)

                        ch(1,i-1,:,2) = wa1(i-2)*(cc(1,i-1,1,:)-cc(1,ic-1,2,:)) &
                            -wa1(i-1)*(cc(1,i,1,:)+cc(1,ic,2,:))
                        ch(1,i,:,2) = wa1(i-2)*(cc(1,i,1,:)+cc(1,ic,2,:))+wa1(i-1) &
                            *(cc(1,i-1,1,:)-cc(1,ic-1,2,:))
                    end do
                    if (mod(ido,2) /= 1) then
                        ch(1,ido,:,1) = cc(1,ido,1,:)+cc(1,ido,1,:)
                        ch(1,ido,:,2) = -(cc(1,1,2,:)+cc(1,1,2,:))
                    end if
            end select
        end if

    end subroutine r1f2kb


    subroutine r1f3kb(ido,l1,cc,in1,ch,in2,wa1,wa2)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        real (wp),    intent (in out) :: cc(in1,ido,3,l1)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,l1,3)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa1(ido)
        real (wp),    intent (in)     :: wa2(ido)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip)         :: i, ic, idp2
        real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
        real (wp), parameter :: ARG =TWO_PI/3
        real (wp), parameter :: TAUR = cos(ARG)
        real (wp), parameter :: TAUI = sin(ARG)
        !------------------------------------------------------------------

        ch(1,1,:,1) = cc(1,1,1,:) + TWO * cc(1,ido,2,:)
        ch(1,1,:,2) = cc(1,1,1,:) + ( TWO * TAUR ) * cc(1,ido,2,:) &
            - ( TWO *TAUI)*cc(1,1,3,:)
        ch(1,1,:,3) = cc(1,1,1,:) + ( TWO *TAUR)*cc(1,ido,2,:) &
            + TWO *TAUI*cc(1,1,3,:)

        select case (ido)
            case (1)
                return
            case default
                idp2 = ido+2
                do i=3,ido,2
                    ic = idp2-i
                    ch(1,i-1,:,1) = cc(1,i-1,1,:)+(cc(1,i-1,3,:)+cc(1,ic-1,2,:))
                    ch(1,i,:,1) = cc(1,i,1,:)+(cc(1,i,3,:)-cc(1,ic,2,:))

                    ch(1,i-1,:,2) = wa1(i-2)* &
                        ((cc(1,i-1,1,:)+TAUR*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))- &
                        (TAUI*(cc(1,i,3,:)+cc(1,ic,2,:)))) &
                        -wa1(i-1)* &
                        ((cc(1,i,1,:)+TAUR*(cc(1,i,3,:)-cc(1,ic,2,:)))+ &
                        (TAUI*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))))

                    ch(1,i,:,2) = wa1(i-2)* &
                        ((cc(1,i,1,:)+TAUR*(cc(1,i,3,:)-cc(1,ic,2,:)))+ &
                        (TAUI*(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))) &
                        +wa1(i-1)* &
                        ((cc(1,i-1,1,:)+TAUR*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))- &
                        (TAUI*(cc(1,i,3,:)+cc(1,ic,2,:))))

                    ch(1,i-1,:,3) = wa2(i-2)* &
                        ((cc(1,i-1,1,:)+TAUR*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))+ &
                        (TAUI*(cc(1,i,3,:)+cc(1,ic,2,:)))) &
                        -wa2(i-1)* &
                        ((cc(1,i,1,:)+TAUR*(cc(1,i,3,:)-cc(1,ic,2,:)))- &
                        (TAUI*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))))

                    ch(1,i,:,3) = wa2(i-2)* &
                        ((cc(1,i,1,:)+TAUR*(cc(1,i,3,:)-cc(1,ic,2,:)))- &
                        (TAUI*(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))) &
                        +wa2(i-1)* &
                        ((cc(1,i-1,1,:)+TAUR*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))+ &
                        (TAUI*(cc(1,i,3,:)+cc(1,ic,2,:))))
                end do
        end select

    end subroutine r1f3kb

    subroutine r1f4kb(ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        real (wp),    intent (in out) :: cc(in1,ido,4,l1)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,l1,4)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa1(ido)
        real (wp),    intent (in)     :: wa2(ido)
        real (wp),    intent (in)     :: wa3(ido)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip)         :: i, ic, idp2
        real (wp), parameter :: SQRT2 = sqrt(TWO)
        !------------------------------------------------------------------


        ch(1,1,:,3) = (cc(1,1,1,:)+cc(1,ido,4,:)) &
            -(cc(1,ido,2,:)+cc(1,ido,2,:))
        ch(1,1,:,1) = (cc(1,1,1,:)+cc(1,ido,4,:)) &
            +(cc(1,ido,2,:)+cc(1,ido,2,:))
        ch(1,1,:,4) = (cc(1,1,1,:)-cc(1,ido,4,:)) &
            +(cc(1,1,3,:)+cc(1,1,3,:))
        ch(1,1,:,2) = (cc(1,1,1,:)-cc(1,ido,4,:)) &
            -(cc(1,1,3,:)+cc(1,1,3,:))

        if (ido < 2) then
            return
        else
            select case (ido)
                case (2)
                    ch(1,ido,:,1) = (cc(1,ido,1,:)+cc(1,ido,3,:)) &
                        +(cc(1,ido,1,:)+cc(1,ido,3,:))
                    ch(1,ido,:,2) = SQRT2*((cc(1,ido,1,:)-cc(1,ido,3,:)) &
                        -(cc(1,1,2,:)+cc(1,1,4,:)))
                    ch(1,ido,:,3) = (cc(1,1,4,:)-cc(1,1,2,:)) &
                        +(cc(1,1,4,:)-cc(1,1,2,:))
                    ch(1,ido,:,4) = -SQRT2*((cc(1,ido,1,:)-cc(1,ido,3,:)) &
                        +(cc(1,1,2,:)+cc(1,1,4,:)))
                case default
                    idp2 = ido+2
                    do i=3,ido,2
                        ic = idp2-i
                        ch(1,i-1,:,1) = (cc(1,i-1,1,:)+cc(1,ic-1,4,:)) &
                            +(cc(1,i-1,3,:)+cc(1,ic-1,2,:))
                        ch(1,i,:,1) = (cc(1,i,1,:)-cc(1,ic,4,:)) &
                            +(cc(1,i,3,:)-cc(1,ic,2,:))
                        ch(1,i-1,:,2)=wa1(i-2)*((cc(1,i-1,1,:)-cc(1,ic-1,4,:)) &
                            -(cc(1,i,3,:)+cc(1,ic,2,:)))-wa1(i-1) &
                            *((cc(1,i,1,:)+cc(1,ic,4,:))+(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))
                        ch(1,i,:,2)=wa1(i-2)*((cc(1,i,1,:)+cc(1,ic,4,:)) &
                            +(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))+wa1(i-1) &
                            *((cc(1,i-1,1,:)-cc(1,ic-1,4,:))-(cc(1,i,3,:)+cc(1,ic,2,:)))
                        ch(1,i-1,:,3)=wa2(i-2)*((cc(1,i-1,1,:)+cc(1,ic-1,4,:)) &
                            -(cc(1,i-1,3,:)+cc(1,ic-1,2,:)))-wa2(i-1) &
                            *((cc(1,i,1,:)-cc(1,ic,4,:))-(cc(1,i,3,:)-cc(1,ic,2,:)))
                        ch(1,i,:,3)=wa2(i-2)*((cc(1,i,1,:)-cc(1,ic,4,:)) &
                            -(cc(1,i,3,:)-cc(1,ic,2,:)))+wa2(i-1) &
                            *((cc(1,i-1,1,:)+cc(1,ic-1,4,:))-(cc(1,i-1,3,:) &
                            +cc(1,ic-1,2,:)))
                        ch(1,i-1,:,4)=wa3(i-2)*((cc(1,i-1,1,:)-cc(1,ic-1,4,:)) &
                            +(cc(1,i,3,:)+cc(1,ic,2,:)))-wa3(i-1) &
                            *((cc(1,i,1,:)+cc(1,ic,4,:))-(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))
                        ch(1,i,:,4)=wa3(i-2)*((cc(1,i,1,:)+cc(1,ic,4,:)) &
                            -(cc(1,i-1,3,:)-cc(1,ic-1,2,:)))+wa3(i-1) &
                            *((cc(1,i-1,1,:)-cc(1,ic-1,4,:))+(cc(1,i,3,:)+cc(1,ic,2,:)))
                    end do

                    if (mod(ido,2) /= 1) then
                        ch(1,ido,:,1) = (cc(1,ido,1,:)+cc(1,ido,3,:)) &
                            +(cc(1,ido,1,:)+cc(1,ido,3,:))
                        ch(1,ido,:,2) = SQRT2*((cc(1,ido,1,:)-cc(1,ido,3,:)) &
                            -(cc(1,1,2,:)+cc(1,1,4,:)))
                        ch(1,ido,:,3) = (cc(1,1,4,:)-cc(1,1,2,:)) &
                            +(cc(1,1,4,:)-cc(1,1,2,:))
                        ch(1,ido,:,4) = -SQRT2*((cc(1,ido,1,:)-cc(1,ido,3,:)) &
                            +(cc(1,1,2,:)+cc(1,1,4,:)))
                    end if
            end select
        end if

    end subroutine r1f4kb


    subroutine r1f5kb(ido,l1,cc,in1,ch,in2,wa1,wa2,wa3,wa4)
        !--------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        real (wp),    intent (in out) :: cc(in1,ido,5,l1)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,l1,5)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa1(ido)
        real (wp),    intent (in)     :: wa2(ido)
        real (wp),    intent (in)     :: wa3(ido)
        real (wp),    intent (in)     :: wa4(ido)
        !--------------------------------------------------
        ! Local variables
        !--------------------------------------------------
        integer (ip)         :: i, ic, idp2
        real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
        real (wp), parameter :: ARG= TWO_PI/5
        real (wp), parameter :: TR11=cos(ARG)
        real (wp), parameter :: TI11=sin(ARG)
        real (wp), parameter :: TR12=cos(TWO*ARG)
        real (wp), parameter :: TI12=sin(TWO*ARG)
        !--------------------------------------------------

        ch(1,1,:,1) = cc(1,1,1,:)+ TWO *cc(1,ido,2,:)+ TWO *cc(1,ido,4,:)
        ch(1,1,:,2) = (cc(1,1,1,:)+TR11* TWO *cc(1,ido,2,:) &
            +TR12* TWO *cc(1,ido,4,:))-(TI11* TWO *cc(1,1,3,:) &
            +TI12* TWO *cc(1,1,5,:))
        ch(1,1,:,3) = (cc(1,1,1,:)+TR12* TWO *cc(1,ido,2,:) &
            +TR11* TWO *cc(1,ido,4,:))-(TI12* TWO *cc(1,1,3,:) &
            -TI11* TWO *cc(1,1,5,:))
        ch(1,1,:,4) = (cc(1,1,1,:)+TR12* TWO *cc(1,ido,2,:) &
            +TR11* TWO *cc(1,ido,4,:))+(TI12* TWO *cc(1,1,3,:) &
            -TI11* TWO *cc(1,1,5,:))
        ch(1,1,:,5) = (cc(1,1,1,:)+TR11* TWO *cc(1,ido,2,:) &
            +TR12* TWO *cc(1,ido,4,:))+(TI11* TWO *cc(1,1,3,:) &
            +TI12* TWO *cc(1,1,5,:))

        select case (ido)
            case (1)
                return
            case default
                idp2 = ido+2
                do i=3,ido,2
                    ic = idp2-i
                    ch(1,i-1,:,1) = cc(1,i-1,1,:)+(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                        +(cc(1,i-1,5,:)+cc(1,ic-1,4,:))
                    ch(1,i,:,1) = cc(1,i,1,:)+(cc(1,i,3,:)-cc(1,ic,2,:)) &
                        +(cc(1,i,5,:)-cc(1,ic,4,:))
                    ch(1,i-1,:,2) = wa1(i-2)*((cc(1,i-1,1,:)+TR11* &
                        (cc(1,i-1,3,:)+cc(1,ic-1,2,:))+TR12 &
                        *(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))-(TI11*(cc(1,i,3,:) &
                        +cc(1,ic,2,:))+TI12*(cc(1,i,5,:)+cc(1,ic,4,:)))) &
                        -wa1(i-1)*((cc(1,i,1,:)+TR11*(cc(1,i,3,:)-cc(1,ic,2,:)) &
                        +TR12*(cc(1,i,5,:)-cc(1,ic,4,:)))+(TI11*(cc(1,i-1,3,:) &
                        -cc(1,ic-1,2,:))+TI12*(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))

                    ch(1,i,:,2) = wa1(i-2)*((cc(1,i,1,:)+TR11*(cc(1,i,3,:) &
                        -cc(1,ic,2,:))+TR12*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                        +(TI11*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))+TI12 &
                        *(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))+wa1(i-1) &
                        *((cc(1,i-1,1,:)+TR11*(cc(1,i-1,3,:) &
                        +cc(1,ic-1,2,:))+TR12*(cc(1,i-1,5,:)+cc(1,ic-1,4,:))) &
                        -(TI11*(cc(1,i,3,:)+cc(1,ic,2,:))+TI12 &
                        *(cc(1,i,5,:)+cc(1,ic,4,:))))

                    ch(1,i-1,:,3) = wa2(i-2) &
                        *((cc(1,i-1,1,:)+TR12*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                        +TR11*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))-(TI12*(cc(1,i,3,:) &
                        +cc(1,ic,2,:))-TI11*(cc(1,i,5,:)+cc(1,ic,4,:)))) &
                        -wa2(i-1) &
                        *((cc(1,i,1,:)+TR12*(cc(1,i,3,:)- &
                        cc(1,ic,2,:))+TR11*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                        +(TI12*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))-TI11 &
                        *(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))

                    ch(1,i,:,3) = wa2(i-2) &
                        *((cc(1,i,1,:)+TR12*(cc(1,i,3,:)- &
                        cc(1,ic,2,:))+TR11*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                        +(TI12*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))-TI11 &
                        *(cc(1,i-1,5,:)-cc(1,ic-1,4,:)))) &
                        +wa2(i-1) &
                        *((cc(1,i-1,1,:)+TR12*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                        +TR11*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))-(TI12*(cc(1,i,3,:) &
                        +cc(1,ic,2,:))-TI11*(cc(1,i,5,:)+cc(1,ic,4,:))))

                    ch(1,i-1,:,4) = wa3(i-2) &
                        *((cc(1,i-1,1,:)+TR12*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                        +TR11*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))+(TI12*(cc(1,i,3,:) &
                        +cc(1,ic,2,:))-TI11*(cc(1,i,5,:)+cc(1,ic,4,:)))) &
                        -wa3(i-1) &
                        *((cc(1,i,1,:)+TR12*(cc(1,i,3,:)- &
                        cc(1,ic,2,:))+TR11*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                        -(TI12*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))-TI11 &
                        *(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))

                    ch(1,i,:,4) = wa3(i-2) &
                        *((cc(1,i,1,:)+TR12*(cc(1,i,3,:)- &
                        cc(1,ic,2,:))+TR11*(cc(1,i,5,:)-cc(1,ic,4,:))) &
                        -(TI12*(cc(1,i-1,3,:)-cc(1,ic-1,2,:))-TI11 &
                        *(cc(1,i-1,5,:)-cc(1,ic-1,4,:)))) &
                        +wa3(i-1) &
                        *((cc(1,i-1,1,:)+TR12*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                        +TR11*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))+(TI12*(cc(1,i,3,:) &
                        +cc(1,ic,2,:))-TI11*(cc(1,i,5,:)+cc(1,ic,4,:))))

                    ch(1,i-1,:,5) = wa4(i-2) &
                        *((cc(1,i-1,1,:)+TR11*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                        +TR12*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))+(TI11*(cc(1,i,3,:) &
                        +cc(1,ic,2,:))+TI12*(cc(1,i,5,:)+cc(1,ic,4,:)))) &
                        -wa4(i-1) &
                        *((cc(1,i,1,:)+TR11*(cc(1,i,3,:)-cc(1,ic,2,:)) &
                        +TR12*(cc(1,i,5,:)-cc(1,ic,4,:)))-(TI11*(cc(1,i-1,3,:) &
                        -cc(1,ic-1,2,:))+TI12*(cc(1,i-1,5,:)-cc(1,ic-1,4,:))))

                    ch(1,i,:,5) = wa4(i-2) &
                        *((cc(1,i,1,:)+TR11*(cc(1,i,3,:)-cc(1,ic,2,:)) &
                        +TR12*(cc(1,i,5,:)-cc(1,ic,4,:)))-(TI11*(cc(1,i-1,3,:) &
                        -cc(1,ic-1,2,:))+TI12*(cc(1,i-1,5,:)-cc(1,ic-1,4,:)))) &
                        +wa4(i-1) &
                        *((cc(1,i-1,1,:)+TR11*(cc(1,i-1,3,:)+cc(1,ic-1,2,:)) &
                        +TR12*(cc(1,i-1,5,:)+cc(1,ic-1,4,:)))+(TI11*(cc(1,i,3,:) &
                        +cc(1,ic,2,:))+TI12*(cc(1,i,5,:)+cc(1,ic,4,:))))
                end do
        end select

    end subroutine r1f5kb


    subroutine r1fgkb(ido,iip,l1,idl1,cc,c1,c2,in1,ch,ch2,in2,wa)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: iip
        integer (ip), intent (in)     :: l1
        integer (ip), intent (in)     :: idl1
        real (wp),    intent (in out) :: c1(in1,ido,l1,iip)
        real (wp),    intent (in out) :: c2(in1,idl1,iip)
        real (wp),    intent (in out) :: cc(in1,ido,iip,l1)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,l1,iip)
        real (wp),    intent (in out) :: ch2(in2,idl1,iip)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa(ido)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        real (wp)            :: ai1, ai2, ar1, ar1h, ar2, ar2h, arg
        real (wp)            :: dc2, dcp, ds2, dsp
        integer (ip)         :: i, idij, idp2, ipp2, ipph
        integer (ip)         :: j, is, j2, jc, k, l, lc, nbd
        real (wp), parameter :: TWO_PI = acos(-ONE)
        !------------------------------------------------------------------

        arg = TWO_PI /iip
        dcp = cos(arg)
        dsp = sin(arg)
        idp2 = ido+2
        nbd = (ido-1)/2
        ipp2 = iip+2
        ipph = (iip+1)/2

        ch(1,:,:,1) = cc(1,:,1,:)

        do j=2,ipph
            jc = ipp2-j
            j2 = j*2
            ch(1,1,:,j) = TWO * cc(1,ido,j2-2,:)
            ch(1,1,:,jc) = TWO * cc(1,1,j2-1,:)
        end do

        if (ido /= 1) then
            if (nbd >= l1) then
                do j=2,ipph
                    jc = ipp2-j
                    ch(1,2:ido-1:2,:, j) = cc(1,2:ido-1:2, 2*j-1,:) + cc(1,idp2-4: &
                        idp2-1-ido:(-2), 2*j-2,:)
                    ch(1,2:ido-1:2,:, jc) = cc(1,2:ido-1:2, 2*j-1,:) - cc(1,idp2-4: &
                        idp2-1-ido:(-2), 2*j-2,:)
                    ch(1,3:ido:2,:, j) = cc(1,3:ido:2, 2*j-1,:) - cc(1,idp2-3:idp2- &
                        ido:(-2), 2*j-2,:)
                    ch(1,3:ido:2,:, jc) = cc(1,3:ido:2, 2*j-1,:) + cc(1,idp2-3:idp2- &
                        ido:(-2), 2*j-2,:)
                end do
            else
                do j=2,ipph
                    jc = ipp2-j
                    ch(1,2:ido-1:2,:, j) = cc(1,2:ido-1:2, 2*j-1,:) + cc(1,idp2-4: &
                        idp2-1-ido:(-2), 2*j-2,:)
                    ch(1,2:ido-1:2,:, jc) = cc(1,2:ido-1:2, 2*j-1,:) - cc(1,idp2-4: &
                        idp2-1-ido:(-2), 2*j-2,:)
                    ch(1,3:ido:2,:, j) = cc(1,3:ido:2, 2*j-1,:) - cc(1,idp2-3:idp2- &
                        ido:(-2), 2*j-2,:)
                    ch(1,3:ido:2,:, jc) = cc(1,3:ido:2, 2*j-1,:) + cc(1,idp2-3:idp2- &
                        ido:(-2), 2*j-2,:)
                end do
            end if
        end if

        ar1 = ONE
        ai1 = ZERO

        do l=2,ipph
            lc = ipp2-l
            ar1h = dcp*ar1-dsp*ai1
            ai1 = dcp*ai1+dsp*ar1
            ar1 = ar1h
            c2(1,:,l) = ch2(1,:,1)+ar1*ch2(1,:,2)
            c2(1,:,lc) = ai1*ch2(1,:,iip)
            dc2 = ar1
            ds2 = ai1
            ar2 = ar1
            ai2 = ai1

            do j=3,ipph
                jc = ipp2-j
                ar2h = dc2*ar2-ds2*ai2
                ai2 = dc2*ai2+ds2*ar2
                ar2 = ar2h
                c2(1,:,l) = c2(1,:,l)+ar2*ch2(1,:,j)
                c2(1,:,lc) = c2(1,:,lc)+ai2*ch2(1,:,jc)
            end do
        end do

        do j=2,ipph
            ch2(1,:,1) = ch2(1,:,1)+ch2(1,:,j)
        end do

        do j=2,ipph
            jc = ipp2-j
            ch(1,1,:,j) = c1(1,1,:,j)-c1(1,1,:,jc)
            ch(1,1,:,jc) = c1(1,1,:,j)+c1(1,1,:,jc)
        end do

        if (ido /= 1) then
            if (nbd >= l1) then
                do j=2,ipph
                    jc = ipp2-j
                    ch(1,2:ido-1:2,:, j) = c1(1,2:ido-1:2,:, j) - c1(1,3:ido:2,:, jc)
                    ch(1,2:ido-1:2,:, jc) = c1(1,2:ido-1:2,:, j) + c1(1,3:ido:2,:, jc)
                    ch(1,3:ido:2,:, j) = c1(1,3:ido:2,:, j) + c1(1,2:ido-1:2,:, jc)
                    ch(1,3:ido:2,:, jc) = c1(1,3:ido:2,:, j) - c1(1,2:ido-1:2,:, jc)
                end do
            else
                do j=2,ipph
                    jc = ipp2-j
                    ch(1,2:ido-1:2,:, j) = c1(1,2:ido-1:2,:, j) - c1(1,3:ido:2,:, jc)
                    ch(1,2:ido-1:2,:, jc) = c1(1,2:ido-1:2,:, j) + c1(1,3:ido:2,:, jc)
                    ch(1,3:ido:2,:, j) = c1(1,3:ido:2,:, j) + c1(1,2:ido-1:2,:, jc)
                    ch(1,3:ido:2,:, jc) = c1(1,3:ido:2,:, j) - c1(1,2:ido-1:2,:, jc)
                end do
            end if
        end if

        select case (ido)
            case (1)
                return
            case default

                c2(1,:,1) = ch2(1,:,1)
                c1(1,1,:,2:iip) = ch(1,1,:,2:iip)

                if ( l1 >= nbd ) then
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        idij = is
                        do i=3,ido,2
                            idij = idij+2
                            c1(1,i-1,:,j) = &
                                wa(idij-1)*ch(1,i-1,:,j) &
                                -wa(idij)*ch(1,i,:,j)
                            c1(1,i,:,j) = &
                                wa(idij-1)*ch(1,i,:,j) &
                                +wa(idij)*ch(1,i-1,:,j)
                        end do
                    end do
                else
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        do k=1,l1
                            idij = is
                            c1(1,2:ido-1:2, k, j) = &
                                wa(idij+1:ido-2+idij:2)*ch(1,2:ido-1:2, &
                                k, j) - wa(idij+2:ido-1+idij:2)*ch(1,3:ido:2, k, j)
                            c1(1,3:ido:2, k, j) = &
                                wa(idij+1:ido-2+idij:2)*ch(1,3:ido:2, k, j) &
                                + wa(idij+2:ido-1+idij:2)*ch(1,2:ido-1:2, k, j)
                        end do
                    end do
                end if
        end select

    end subroutine r1fgkb



end submodule real_backward_1d
