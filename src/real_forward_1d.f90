submodule (real_transform_routines) real_forward_1d

contains

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

        !
        !==> Check validity of input arguments
        !
        if (lenr < inc*(n-1) + 1) then
            ierror = 1
            call fft_error_handler('rfft1f ', 6)
        else if (lensav < n + int(log(real(n, kind=wp) )/log(TWO)) +4) then
            ierror = 2
            call fft_error_handler('rfft1f ', 8)
        else if (lenwrk < n) then
            ierror = 3
            call fft_error_handler('rfft1f ', 10)
        else
            ierror = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) call real_pass_forward(n,inc,r,work,wsave,wsave(n+1))


    end subroutine rfft1f



    subroutine real_pass_forward(n, in, c, ch, wa, fac)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: in
        real (wp),    intent (in out) :: c(in,*)
        real (wp),    intent (out)    :: ch(:)
        real (wp),    intent (in)     :: wa(n)
        real (wp),    intent (in)     :: fac(15)
        !--------------------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------------------
        integer (ip) :: idl1, ido, iip
        integer (ip) :: workspace_indices(4)
        integer (ip) :: j,  k1, kh, l1, l2
        integer (ip) :: modn, na, nf, nl
        real (wp)    :: sn, tsn, tsnm
        !--------------------------------------------------------------

        nf = int(fac(2), kind=ip)
        na = 1
        l2 = n

        associate( &
            iw1 => workspace_indices(1), &
            iw2 => workspace_indices(2), &
            iw3 => workspace_indices(3), &
            iw4 => workspace_indices(4) &
            )

            iw1 = n

            do k1=1,nf

                kh = nf-k1
                iip = int(fac(kh+3), kind=ip)
                l1 = l2/iip
                ido = n/l2
                idl1 = ido*l1
                iw1 = iw1-(iip-1)*ido
                na = 1-na

                select case (iip)
                    case (2)
                        select case (na)
                            case (0)
                                call real_pass_2_forward(ido,l1,c,in,ch,1,wa(iw1))
                            case default
                                call real_pass_2_forward(ido,l1,ch,1,c,in,wa(iw1))
                        end select
                    case (3)
                        iw2 = iw1+ido
                        select case (na)
                            case (0)
                                call real_pass_3_forward(ido,l1,c,in,ch,1,wa(iw1),wa(iw2))
                            case default
                                call real_pass_3_forward(ido,l1,ch,1,c,in,wa(iw1),wa(iw2))
                        end select
                    case (4)
                        iw2 = iw1+ido
                        iw3 = iw2+ido
                        select case (na)
                            case (0)
                                call real_pass_4_forward(ido,l1,c,in,ch,1,wa(iw1),wa(iw2),wa(iw3))
                            case default
                                call real_pass_4_forward(ido,l1,ch,1,c,in,wa(iw1),wa(iw2),wa(iw3))
                        end select
                    case (5)
                        iw2 = iw1+ido
                        iw3 = iw2+ido
                        iw4 = iw3+ido
                        select case (na)
                            case (0)
                                call real_pass_5_forward(ido,l1,c,in,ch,1,wa(iw1),wa(iw2),wa(iw3),wa(iw4))
                            case default
                                call real_pass_5_forward(ido,l1,ch,1,c,in,wa(iw1),wa(iw2),wa(iw3),wa(iw4))
                        end select
                    case default
                        if (ido == 1) na = 1-na
                        select case (na)
                            case (0)
                                call real_pass_n_forward(ido,iip,l1,idl1,c,c,c,in,ch,ch,1,wa(iw1))
                                na = 1
                            case default
                                call real_pass_n_forward(ido,iip,l1,idl1,ch,ch,ch,1,c,c,in,wa(iw1))
                                na = 0
                        end select
                end select
                l2 = l1
            end do

        end associate

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
            c(1,1) = sn*ch(1)
            do j=2,nl,2
                c(1,j) = tsn*ch(j)
                c(1,j+1) = tsnm*ch(j+1)
            end do
            if (modn == 0) c(1,n) = sn*ch(n)
        else
            c(1,1) = sn*c(1,1)
            do j=2,nl,2
                c(1,j) = tsn*c(1,j)
                c(1,j+1) = tsnm*c(1,j+1)
            end do
            if (modn == 0) c(1,n) = sn*c(1,n)
        end if

    end subroutine real_pass_forward


    subroutine real_pass_2_forward(ido,l1,cc,in1,ch,in2,wa1)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        real (wp),    intent (in out) :: cc(in1,ido,l1,2)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,2,l1)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa1(ido)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer (ip) :: i, ic, idp2
        !-----------------------------------------------

        ch(1,1,1,:) = cc(1,1,:,1)+cc(1,1,:,2)
        ch(1,ido,2,:) = cc(1,1,:,1)-cc(1,1,:,2)

        if (ido < 2) return

        select case (ido)
            case (2)
                ch(1,1,2,:) = -cc(1,ido,:,2)
                ch(1,ido,1,:) = cc(1,ido,:,1)
            case default
                idp2 = ido+2
                do i=3,ido,2
                    ic = idp2-i

                    ch(1,i,1,:) = &
                        cc(1,i,:,1)+(wa1(i-2)*cc(1,i,:,2) &
                        -wa1(i-1)*cc(1,i-1,:,2))

                    ch(1,ic,2,:) = &
                        (wa1(i-2)*cc(1,i,:,2) &
                        -wa1(i-1)*cc(1,i-1,:,2))-cc(1,i,:,1)

                    ch(1,i-1,1,:) = &
                        cc(1,i-1,:,1)+(wa1(i-2)*cc(1,i-1,:,2) &
                        +wa1(i-1)*cc(1,i,:,2))

                    ch(1,ic-1,2,:) = &
                        cc(1,i-1,:,1)-(wa1(i-2)*cc(1,i-1,:,2) &
                        +wa1(i-1)*cc(1,i,:,2))
                end do

                if (mod(ido,2) /= 1) then
                    ch(1,1,2,:) = -cc(1,ido,:,2)
                    ch(1,ido,1,:) = cc(1,ido,:,1)
                end if
        end select

    end subroutine real_pass_2_forward


    subroutine real_pass_3_forward(ido,l1,cc,in1,ch,in2,wa1,wa2)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        real (wp),    intent (in out) :: cc(in1,ido,l1,3)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,3,l1)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa1(ido)
        real (wp),    intent (in)     :: wa2(ido)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer (ip)         :: i, ic, idp2
        real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
        real (wp), parameter :: ARG = TWO_PI/3
        real (wp), parameter :: TAUR = cos(ARG)
        real (wp), parameter :: TAUI = sin(ARG)
        !-----------------------------------------------

        ch(1,1,1,:) = cc(1,1,:,1)+(cc(1,1,:,2)+cc(1,1,:,3))
        ch(1,1,3,:) = TAUI*(cc(1,1,:,3)-cc(1,1,:,2))
        ch(1,ido,2,:) = cc(1,1,:,1)+TAUR*(cc(1,1,:,2)+cc(1,1,:,3))

        select case (ido)
            case (1)
                return
            case default
                idp2 = ido+2
                do i=3,ido,2
                    ic = idp2-i

                    ch(1,i-1,1,:) = cc(1,i-1,:,1)+((wa1(i-2)*cc(1,i-1,:,2)+ &
                        wa1(i-1)*cc(1,i,:,2))+(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                        cc(1,i,:,3)))

                    ch(1,i,1,:) = cc(1,i,:,1)+((wa1(i-2)*cc(1,i,:,2)- &
                        wa1(i-1)*cc(1,i-1,:,2))+(wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                        cc(1,i-1,:,3)))

                    ch(1,i-1,3,:) = (cc(1,i-1,:,1)+TAUR*((wa1(i-2)* &
                        cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2))+(wa2(i-2)* &
                        cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3))))+(TAUI*((wa1(i-2)* &
                        cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2))-(wa2(i-2)* &
                        cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3))))

                    ch(1,ic-1,2,:) = (cc(1,i-1,:,1)+TAUR*((wa1(i-2)* &
                        cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2))+(wa2(i-2)* &
                        cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3))))-(TAUI*((wa1(i-2)* &
                        cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2))-(wa2(i-2)* &
                        cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3))))

                    ch(1,i,3,:) = (cc(1,i,:,1)+TAUR*((wa1(i-2)*cc(1,i,:,2)- &
                        wa1(i-1)*cc(1,i-1,:,2))+(wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                        cc(1,i-1,:,3))))+(TAUI*((wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                        cc(1,i,:,3))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2))))

                    ch(1,ic,2,:) = (TAUI*((wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                        cc(1,i,:,3))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2))))-(cc(1,i,:,1)+TAUR*((wa1(i-2)*cc(1,i,:,2)- &
                        wa1(i-1)*cc(1,i-1,:,2))+(wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                        cc(1,i-1,:,3))))
                end do
        end select

    end subroutine real_pass_3_forward



    subroutine real_pass_4_forward(ido,l1,cc,in1,ch,in2,wa1,wa2,wa3)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        real (wp),    intent (in out) :: cc(in1,ido,l1,4)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,4,l1)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa1(ido)
        real (wp),    intent (in)     :: wa2(ido)
        real (wp),    intent (in)     :: wa3(ido)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer (ip)         :: i, ic, idp2
        real (wp), parameter :: HALF_SQRT2=sqrt(TWO)/2
        !-----------------------------------------------

        ch(1,1,1,:) = (cc(1,1,:,2)+cc(1,1,:,4))+(cc(1,1,:,1)+cc(1,1,:,3))
        ch(1,ido,4,:) = (cc(1,1,:,1)+cc(1,1,:,3))-(cc(1,1,:,2)+cc(1,1,:,4))
        ch(1,ido,2,:) = cc(1,1,:,1)-cc(1,1,:,3)
        ch(1,1,3,:) = cc(1,1,:,4)-cc(1,1,:,2)

        if (ido < 2) then
            return
        end if

        select case (ido)
            case (2)
                ch(1,ido,1,:) = (HALF_SQRT2*(cc(1,ido,:,2)-cc(1,ido,:,4)))+cc(1,ido,:,1)
                ch(1,ido,3,:) = cc(1,ido,:,1)-(HALF_SQRT2*(cc(1,ido,:,2)-cc(1,ido,:,4)))
                ch(1,1,2,:) = (-HALF_SQRT2*(cc(1,ido,:,2)+cc(1,ido,:,4)))-cc(1,ido,:,3)
                ch(1,1,4,:) = (-HALF_SQRT2*(cc(1,ido,:,2)+cc(1,ido,:,4)))+cc(1,ido,:,3)
            case default
                idp2 = ido+2
                do i=3,ido,2
                    ic = idp2-i
                    ch(1,i-1,1,:) = ((wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2))+(wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                        cc(1,i,:,4)))+(cc(1,i-1,:,1)+(wa2(i-2)*cc(1,i-1,:,3)+ &
                        wa2(i-1)*cc(1,i,:,3)))
                    ch(1,ic-1,4,:) = (cc(1,i-1,:,1)+(wa2(i-2)*cc(1,i-1,:,3)+ &
                        wa2(i-1)*cc(1,i,:,3)))-((wa1(i-2)*cc(1,i-1,:,2)+ &
                        wa1(i-1)*cc(1,i,:,2))+(wa3(i-2)*cc(1,i-1,:,4)+ &
                        wa3(i-1)*cc(1,i,:,4)))
                    ch(1,i,1,:) = ((wa1(i-2)*cc(1,i,:,2)-wa1(i-1)* &
                        cc(1,i-1,:,2))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4)))+(cc(1,i,:,1)+(wa2(i-2)*cc(1,i,:,3)- &
                        wa2(i-1)*cc(1,i-1,:,3)))
                    ch(1,ic,4,:) = ((wa1(i-2)*cc(1,i,:,2)-wa1(i-1)* &
                        cc(1,i-1,:,2))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4)))-(cc(1,i,:,1)+(wa2(i-2)*cc(1,i,:,3)- &
                        wa2(i-1)*cc(1,i-1,:,3)))
                    ch(1,i-1,3,:) = ((wa1(i-2)*cc(1,i,:,2)-wa1(i-1)* &
                        cc(1,i-1,:,2))-(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4)))+(cc(1,i-1,:,1)-(wa2(i-2)*cc(1,i-1,:,3)+ &
                        wa2(i-1)*cc(1,i,:,3)))
                    ch(1,ic-1,2,:) = (cc(1,i-1,:,1)-(wa2(i-2)*cc(1,i-1,:,3)+ &
                        wa2(i-1)*cc(1,i,:,3)))-((wa1(i-2)*cc(1,i,:,2)-wa1(i-1)* &
                        cc(1,i-1,:,2))-(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4)))
                    ch(1,i,3,:) = ((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                        cc(1,i,:,4))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2)))+(cc(1,i,:,1)-(wa2(i-2)*cc(1,i,:,3)- &
                        wa2(i-1)*cc(1,i-1,:,3)))
                    ch(1,ic,2,:) = ((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                        cc(1,i,:,4))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2)))-(cc(1,i,:,1)-(wa2(i-2)*cc(1,i,:,3)- &
                        wa2(i-1)*cc(1,i-1,:,3)))
                end do

                if (mod(ido,2) /= 1) then
                    ch(1,ido,1,:) = (HALF_SQRT2*(cc(1,ido,:,2)-cc(1,ido,:,4)))+cc(1,ido,:,1)
                    ch(1,ido,3,:) = cc(1,ido,:,1)-(HALF_SQRT2*(cc(1,ido,:,2)-cc(1,ido,:,4)))
                    ch(1,1,2,:) = (-HALF_SQRT2*(cc(1,ido,:,2)+cc(1,ido,:,4)))-cc(1,ido,:,3)
                    ch(1,1,4,:) = (-HALF_SQRT2*(cc(1,ido,:,2)+cc(1,ido,:,4)))+cc(1,ido,:,3)
                end if
        end select

    end subroutine real_pass_4_forward



    subroutine real_pass_5_forward(ido,l1,cc,in1,ch,in2,wa1,wa2,wa3,wa4)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        real (wp),    intent (in out) :: cc(in1,ido,l1,5)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,5,l1)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa1(ido)
        real (wp),    intent (in)     :: wa2(ido)
        real (wp),    intent (in)     :: wa3(ido)
        real (wp),    intent (in)     :: wa4(ido)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        integer (ip)         :: i, ic, idp2
        real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
        real (wp), parameter :: ARG= TWO_PI/5
        real (wp), parameter :: TR11=cos(ARG)
        real (wp), parameter :: TI11=sin(ARG)
        real (wp), parameter :: TR12=cos(TWO *ARG)
        real (wp), parameter :: TI12=sin(TWO *ARG)
        !-----------------------------------------------

        ch(1,1,1,:) = cc(1,1,:,1)+(cc(1,1,:,5)+cc(1,1,:,2))+ &
            (cc(1,1,:,4)+cc(1,1,:,3))
        ch(1,ido,2,:) = cc(1,1,:,1)+TR11*(cc(1,1,:,5)+cc(1,1,:,2))+ &
            TR12*(cc(1,1,:,4)+cc(1,1,:,3))
        ch(1,1,3,:) = TI11*(cc(1,1,:,5)-cc(1,1,:,2))+TI12* &
            (cc(1,1,:,4)-cc(1,1,:,3))
        ch(1,ido,4,:) = cc(1,1,:,1)+TR12*(cc(1,1,:,5)+cc(1,1,:,2))+ &
            TR11*(cc(1,1,:,4)+cc(1,1,:,3))
        ch(1,1,5,:) = TI12*(cc(1,1,:,5)-cc(1,1,:,2))-TI11* &
            (cc(1,1,:,4)-cc(1,1,:,3))

        select case (ido)
            case (1)
                return
            case default
                idp2 = ido+2
                do i=3,ido,2
                    ic = idp2-i
                    ch(1,i-1,1,:) = cc(1,i-1,:,1)+((wa1(i-2)*cc(1,i-1,:,2)+ &
                        wa1(i-1)*cc(1,i,:,2))+(wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)* &
                        cc(1,i,:,5)))+((wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                        cc(1,i,:,3))+(wa3(i-2)*cc(1,i-1,:,4)+ &
                        wa3(i-1)*cc(1,i,:,4)))

                    ch(1,i,1,:) = cc(1,i,:,1)+((wa1(i-2)*cc(1,i,:,2)- &
                        wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                        cc(1,i-1,:,5)))+((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                        cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4)))

                    ch(1,i-1,3,:) = cc(1,i-1,:,1)+TR11* &
                        ( wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2) &
                        +wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)*cc(1,i,:,5))+TR12* &
                        ( wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3) &
                        +wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)*cc(1,i,:,4))+TI11* &
                        ( wa1(i-2)*cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2) &
                        -(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)*cc(1,i-1,:,5)))+TI12* &
                        ( wa2(i-2)*cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3) &
                        -(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)*cc(1,i-1,:,4)))

                    ch(1,ic-1,2,:) = cc(1,i-1,:,1)+TR11* &
                        ( wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2) &
                        +wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)*cc(1,i,:,5))+TR12* &
                        ( wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3) &
                        +wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)*cc(1,i,:,4))-(TI11* &
                        ( wa1(i-2)*cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2) &
                        -(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)*cc(1,i-1,:,5)))+TI12* &
                        ( wa2(i-2)*cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3) &
                        -(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)*cc(1,i-1,:,4))))

                    ch(1,i,3,:) = (cc(1,i,:,1)+TR11*((wa1(i-2)*cc(1,i,:,2)- &
                        wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                        cc(1,i-1,:,5)))+TR12*((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                        cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4))))+(TI11*((wa4(i-2)*cc(1,i-1,:,5)+ &
                        wa4(i-1)*cc(1,i,:,5))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2)))+TI12*((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                        cc(1,i,:,4))-(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                        cc(1,i,:,3))))

                    ch(1,ic,2,:) = (TI11*((wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)* &
                        cc(1,i,:,5))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2)))+TI12*((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                        cc(1,i,:,4))-(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                        cc(1,i,:,3))))-(cc(1,i,:,1)+TR11*((wa1(i-2)*cc(1,i,:,2)- &
                        wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                        cc(1,i-1,:,5)))+TR12*((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                        cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4))))

                    ch(1,i-1,5,:) = (cc(1,i-1,:,1)+TR12*((wa1(i-2)* &
                        cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2))+(wa4(i-2)* &
                        cc(1,i-1,:,5)+wa4(i-1)*cc(1,i,:,5)))+TR11*((wa2(i-2)* &
                        cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3))+(wa3(i-2)* &
                        cc(1,i-1,:,4)+wa3(i-1)*cc(1,i,:,4))))+(TI12*((wa1(i-2)* &
                        cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2))-(wa4(i-2)* &
                        cc(1,i,:,5)-wa4(i-1)*cc(1,i-1,:,5)))-TI11*((wa2(i-2)* &
                        cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3))-(wa3(i-2)* &
                        cc(1,i,:,4)-wa3(i-1)*cc(1,i-1,:,4))))

                    ch(1,ic-1,4,:) = (cc(1,i-1,:,1)+TR12*((wa1(i-2)* &
                        cc(1,i-1,:,2)+wa1(i-1)*cc(1,i,:,2))+(wa4(i-2)* &
                        cc(1,i-1,:,5)+wa4(i-1)*cc(1,i,:,5)))+TR11*((wa2(i-2)* &
                        cc(1,i-1,:,3)+wa2(i-1)*cc(1,i,:,3))+(wa3(i-2)* &
                        cc(1,i-1,:,4)+wa3(i-1)*cc(1,i,:,4))))-(TI12*((wa1(i-2)* &
                        cc(1,i,:,2)-wa1(i-1)*cc(1,i-1,:,2))-(wa4(i-2)* &
                        cc(1,i,:,5)-wa4(i-1)*cc(1,i-1,:,5)))-TI11*((wa2(i-2)* &
                        cc(1,i,:,3)-wa2(i-1)*cc(1,i-1,:,3))-(wa3(i-2)* &
                        cc(1,i,:,4)-wa3(i-1)*cc(1,i-1,:,4))))

                    ch(1,i,5,:) = (cc(1,i,:,1)+TR12*((wa1(i-2)*cc(1,i,:,2)- &
                        wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                        cc(1,i-1,:,5)))+TR11*((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                        cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4))))+(TI12*((wa4(i-2)*cc(1,i-1,:,5)+ &
                        wa4(i-1)*cc(1,i,:,5))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2)))-TI11*((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                        cc(1,i,:,4))-(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                        cc(1,i,:,3))))

                    ch(1,ic,4,:) = (TI12*((wa4(i-2)*cc(1,i-1,:,5)+wa4(i-1)* &
                        cc(1,i,:,5))-(wa1(i-2)*cc(1,i-1,:,2)+wa1(i-1)* &
                        cc(1,i,:,2)))-TI11*((wa3(i-2)*cc(1,i-1,:,4)+wa3(i-1)* &
                        cc(1,i,:,4))-(wa2(i-2)*cc(1,i-1,:,3)+wa2(i-1)* &
                        cc(1,i,:,3))))-(cc(1,i,:,1)+TR12*((wa1(i-2)*cc(1,i,:,2)- &
                        wa1(i-1)*cc(1,i-1,:,2))+(wa4(i-2)*cc(1,i,:,5)-wa4(i-1)* &
                        cc(1,i-1,:,5)))+TR11*((wa2(i-2)*cc(1,i,:,3)-wa2(i-1)* &
                        cc(1,i-1,:,3))+(wa3(i-2)*cc(1,i,:,4)-wa3(i-1)* &
                        cc(1,i-1,:,4))))
                end do
        end select

    end subroutine real_pass_5_forward


    subroutine real_pass_n_forward(ido,iip,l1,idl1,cc,c1,c2,in1,ch,ch2,in2,wa)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: iip
        integer (ip), intent (in)     :: l1
        integer (ip), intent (in)     :: idl1
        real (wp),    intent (out)    :: cc(in1,ido,iip,l1)
        real (wp),    intent (in out) :: c1(in1,ido,l1,iip)
        real (wp),    intent (in out) :: c2(in1,idl1,iip)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,ido,l1,iip)
        real (wp),    intent (in out) :: ch2(in2,idl1,iip)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa(ido)
        !-----------------------------------------------
        ! Local variables
        !-----------------------------------------------
        real (wp)            :: ai1, ai2, ar1, ar1h
        real (wp)            :: ar2, ar2h, arg
        real (wp)            :: dc2, dcp, ds2, dsp
        integer (ip)         :: i, idij, idp2
        integer (ip)         :: ipp2, ipph, is
        integer (ip)         :: j, jc
        integer (ip)         :: k, l, lc, nbd
        real (wp), parameter :: TWO_PI = TWO * acos(-ONE)
        !-----------------------------------------------

        arg = TWO_PI/iip
        dcp = cos(arg)
        dsp = sin(arg)
        ipph = (iip+1)/2
        ipp2 = iip+2
        idp2 = ido+2
        nbd = (ido-1)/2


        select case (ido)
            case (1)
                c2(1,:,1) = ch2(1,:,1)
            case default
                ch2(1,:,1) = c2(1,:,1)
                ch(1,1,:,2:iip) = c1(1,1,:,2:iip)
                if ( nbd <= l1 ) then
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        idij = is
                        do i=3,ido,2
                            idij = idij+2
                            ch(1,i-1,:,j) = wa(idij-1)*c1(1,i-1,:,j)+wa(idij)*c1(1,i,:,j)
                            ch(1,i,:,j) = wa(idij-1)*c1(1,i,:,j)-wa(idij)*c1(1,i-1,:,j)
                        end do
                    end do
                else
                    is = -ido
                    do j=2,iip
                        is = is+ido
                        do k=1,l1
                            idij = is
                            ch(1,2:ido-1:2,k,j) = &
                                wa(idij+1:ido-2+idij:2)*c1(1,2:ido-1:2,k,j)&
                                +wa(idij+2:ido-1+idij:2)*c1(1,3:ido:2,k,j)
                            ch(1,3:ido:2,k,j) = &
                                wa(idij+1:ido-2+idij:2)*c1(1,3:ido:2,k,j)&
                                -wa(idij+2:ido-1+idij:2)*c1(1,2:ido-1:2,k,j)
                        end do
                    end do
                end if

                if (l1 <= nbd) then
                    do j=2,ipph
                        jc = ipp2-j
                        c1(1,2:ido-1:2,:, j)=ch(1,2:ido-1:2,:, j)+ch(1,2:ido-1:2,:, jc)
                        c1(1,2:ido-1:2,:, jc) = ch(1,3:ido:2,:, j) - ch(1,3:ido:2,:, jc)
                        c1(1,3:ido:2,:, j) = ch(1,3:ido:2,:, j) + ch(1,3:ido:2,:, jc)
                        c1(1,3:ido:2,:, jc) = ch(1,2:ido-1:2,:, jc) - ch(1,2:ido-1:2,:, j)
                    end do
                else
                    do j=2,ipph
                        jc = ipp2-j
                        c1(1,2:ido-1:2,:, j) = ch(1,2:ido-1:2,:, j) + ch(1,2:ido-1:2,:, jc)
                        c1(1,2:ido-1:2,:, jc) = ch(1,3:ido:2,:, j) - ch(1,3:ido:2,:, jc)
                        c1(1,3:ido:2,:, j) = ch(1,3:ido:2,:, j) + ch(1,3:ido:2,:, jc)
                        c1(1,3:ido:2,:, jc) = ch(1,2:ido-1:2,:, jc) - ch(1,2:ido-1:2,:, j)
                    end do
                end if
        end select

        do j=2,ipph
            jc = ipp2-j
            c1(1,1,:,j) = ch(1,1,:,j)+ch(1,1,:,jc)
            c1(1,1,:,jc) = ch(1,1,:,jc)-ch(1,1,:,j)
        end do

        ar1 = ONE
        ai1 = ZERO
        do l=2,ipph
            lc = ipp2-l
            ar1h = dcp*ar1-dsp*ai1
            ai1 = dcp*ai1+dsp*ar1
            ar1 = ar1h
            ch2(1,:,l) = c2(1,:,1)+ar1*c2(1,:,2)
            ch2(1,:,lc) = ai1*c2(1,:,iip)
            dc2 = ar1
            ds2 = ai1
            ar2 = ar1
            ai2 = ai1
            do j=3,ipph
                jc = ipp2-j
                ar2h = dc2*ar2-ds2*ai2
                ai2 = dc2*ai2+ds2*ar2
                ar2 = ar2h
                ch2(1,:,l) = ch2(1,:,l)+ar2*c2(1,:,j)
                ch2(1,:,lc) = ch2(1,:,lc)+ai2*c2(1,:,jc)
            end do
        end do

        do j=2,ipph
            ch2(1,:,1) = ch2(1,:,1)+c2(1,:,j)
        end do

        cc(1,:,1,:) = ch(1,:,:,1)

        cc(1,ido,2:(ipph-1)*2:2,:) = transpose(ch(1,1,:,2:ipph))

        cc(1,1,3:ipph*2-1:2,:) = transpose(ch(1,1,:,ipp2-2:ipp2-ipph:(-1)))

        select case (ido)
            case (1)
                return
            case default
                if (l1 <= nbd) then

                    cc(1,2:ido-1:2, 3:ipph*2-1:2,:) = &
                        reshape(source = &
                        ch(1,2:ido-1:2,:,2:ipph)+ch(1,2:ido-1:2,:,ipp2-2:ipp2-ipph:(-1)), &
                        shape = [(ido -1)/2, ipph-1, l1], &
                        order = [1, 3, 2])

                    cc(1,idp2-4:idp2-1-ido:(-2), 2:(ipph-1)*2:2,:) = &
                        reshape(source = &
                        ch(1,2:ido-1:2,:, 2:ipph)-ch(1,2:ido-1:2,:,ipp2-2:ipp2-ipph:(-1)), &
                        shape = [(ido-1)/2, ipph-1, l1], &
                        order = [1, 3, 2])

                    cc(1,3:ido:2, 3:ipph*2-1:2,:) = &
                        reshape(source = &
                        ch(1,3:ido:2,:,2:ipph)+ch(1,3:ido:2,:,ipp2-2:ipp2-ipph:(-1)), &
                        shape = [(ido-1)/2, ipph-1, l1], &
                        order = [1, 3, 2])

                    cc(1,idp2-3:idp2-ido:(-2), 2:(ipph-1)*2:2,:) = &
                        reshape(source = &
                        ch(1,3:ido:2,:,ipp2-2:ipp2-ipph:(-1))-ch(1,3:ido:2,:, 2:ipph), shape &
                        = [(ido-1)/2, ipph-1, l1], order = [1, 3, 2])
                else
                    cc(1,2:ido-1:2, 3:ipph*2-1:2,:) = &
                        reshape(source = &
                        ch(1,2:ido-1:2,:, 2:ipph)+ch(1,2:ido-1:2,:,ipp2-2:ipp2-ipph:(-1)), &
                        shape = [(ido-1)/2, ipph-1, l1], &
                        order = [1, 3, 2])

                    cc(1,idp2-4:idp2-1-ido:(-2), 2:(ipph-1)*2:2,:) = &
                        reshape(source = &
                        ch(1,2:ido-1:2,:, 2:ipph)-ch(1,2:ido-1:2,:,ipp2-2:ipp2-ipph:(-1)), &
                        shape = [(ido-1)/2, ipph-1, l1], &
                        order = [1, 3, 2])

                    cc(1,3:ido:2, 3:ipph*2-1:2,:) = &
                        reshape(source = &
                        ch(1,3:ido:2,:, 2:ipph) +ch(1,3:ido:2,:, ipp2-2:ipp2-ipph:(-1)), &
                        shape = [(ido-1)/2, ipph-1, l1], &
                        order = [1, 3, 2])

                    cc(1,idp2-3:idp2-ido:(-2), 2:(ipph-1)*2:2,:) = &
                        reshape(source = &
                        ch(1,3:ido:2,:,ipp2-2:ipp2-ipph:(-1))-ch(1,3:ido:2,:,2:ipph), &
                        shape = [(ido-1)/2, ipph-1, l1], order = [1, 3, 2])
                end if
        end select

    end subroutine real_pass_n_forward



end submodule real_forward_1d
