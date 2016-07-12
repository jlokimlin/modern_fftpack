submodule (complex_transform_routines) complex_backward_1d

contains

    module subroutine cfft1b(n, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
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
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: inc
        complex (wp), intent (in out) :: c(lenc)
        integer (ip), intent (in)     :: lenc
        integer (ip), intent (in)     :: lensav
        integer (ip), intent (in)     :: lenwrk
        integer (ip), intent (out)    :: ier
        real (wp),    intent (in out) :: work(lenwrk)
        real (wp),    intent (in out) :: wsave(lensav)
        !--------------------------------------------------------------
        ! Local variables
        !--------------------------------------------------------------
        integer (ip)           :: iw1
        real (wp), allocatable :: real_copy(:,:)
        !--------------------------------------------------------------

        !
        !==> Check validity of calling arguments
        !
        if (lenc < inc * (n - 1) + 1) then
            ier = 1
            call fft_error_handler('cfft1b ', 4)
        else if (size(wsave) < get_complex_1d_saved_workspace_length(n)) then
            ier = 2
            call fft_error_handler('cfft1b ', 6)
        else if (size(work) < get_complex_1d_workspace_length(n)) then
            ier = 3
            call fft_error_handler('cfft1b ', 8)
        else
            ier = 0
        end if

        !
        !==> Perform transform
        !
        if (n /= 1) then
            !
            !==> Allocate memory
            !
            allocate( real_copy(2,size(c)) )

            !
            !==> Make copy: complex to real
            !
            real_copy(1,:) = real(c)
            real_copy(2,:) = aimag(c)

            ! Set workspace index pointer
            iw1 = 2 * n + 1

            call c1fm1b(n, inc, real_copy, work, wsave, wsave(iw1), wsave(iw1+1:))

            !
            !==> Make copy: real to complex
            !
            c =  cmplx(real_copy(1,:), real_copy(2,:), kind=wp)

            !
            !==> Release memory
            !
            deallocate( real_copy )
        end if


    end subroutine cfft1b


    subroutine c1fm1b(n, inc, c, ch, wa, fnf, fac)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: inc
        real (wp),    intent (in out) :: c(:, :)
        real (wp),    intent (out)    :: ch(:)
        real (wp),    intent (in out) :: wa(*)
        real (wp),    intent (in out) :: fnf
        real (wp),    intent (in out) :: fac(:)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip) :: ido, inc2, iip, iw
        integer (ip) :: k1, l1, l2, lid
        integer (ip) :: na, nbr, nf
        !----------------------------------------------------------------------

        inc2 = 2*inc
        nf = int(fnf, kind=ip)
        na = 0
        l1 = 1
        iw = 1

        do k1=1, nf
            iip = int(fac(k1), kind=ip)
            l2 = iip*l1
            ido = n/l2
            lid = l1*ido
            nbr = 1+na+2*min(iip-2, 4)
            select case (nbr)
                case (1)
                    call c1f2kb(ido, l1, na, c, inc2, ch, 2, wa(iw))
                case (2)
                    call c1f2kb(ido, l1, na, ch, 2, c, inc2, wa(iw))
                case (3)
                    call c1f3kb(ido, l1, na, c, inc2, ch, 2, wa(iw))
                case (4)
                    call c1f3kb(ido, l1, na, ch, 2, c, inc2, wa(iw))
                case (5)
                    call c1f4kb(ido, l1, na, c, inc2, ch, 2, wa(iw))
                case (6)
                    call c1f4kb(ido, l1, na, ch, 2, c, inc2, wa(iw))
                case (7)
                    call c1f5kb(ido, l1, na, c, inc2, ch, 2, wa(iw))
                case (8)
                    call c1f5kb(ido, l1, na, ch, 2, c, inc2, wa(iw))
                case (9)
                    call c1fgkb(ido, iip, l1, lid, na, c, c, inc2, ch, ch, 2, wa(iw))
                case (10)
                    call c1fgkb(ido, iip, l1, lid, na, ch, ch, 2, c, c, inc2, wa(iw))
            end select

            l1 = l2
            iw = iw+(iip-1)*(2*ido)

            if (iip <= 5) na = 1-na

        end do

    end subroutine c1fm1b



    subroutine c1f2kb(ido, l1, na, cc, in1, ch, in2, wa)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        integer (ip), intent (in)     :: na
        real (wp),    intent (in out) :: cc(in1,l1,ido,2)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,l1,2,ido)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa(ido,1,2)
        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer (ip)           :: i !! Counter
        real (wp), allocatable :: chold1(:), chold2(:)
        real (wp), allocatable :: ti2(:), tr2(:)
        !----------------------------------------------------------------------

        if (ido <= 1 .and. na /= 1) then
            !
            !==> Allocate memory
            !
            allocate( chold1(l1) )
            allocate( chold2(l1) )

            chold1 = cc(1,:,1,1)+cc(1,:,1,2)
            cc(1,:,1,2) = cc(1,:,1,1)-cc(1,:,1,2)
            cc(1,:,1,1) = chold1
            chold2 = cc(2,:,1,1)+cc(2,:,1,2)
            cc(2,:,1,2) = cc(2,:,1,1)-cc(2,:,1,2)
            cc(2,:,1,1) = chold2
            !
            !==> Release memory
            !
            deallocate( chold1 )
            deallocate( chold2 )
        else
            ch(1,:,1,1) = cc(1,:,1,1)+cc(1,:,1,2)
            ch(1,:,2,1) = cc(1,:,1,1)-cc(1,:,1,2)
            ch(2,:,1,1) = cc(2,:,1,1)+cc(2,:,1,2)
            ch(2,:,2,1) = cc(2,:,1,1)-cc(2,:,1,2)
            !
            !==> Allocate memory
            !
            allocate( tr2(l1) )
            allocate( ti2(l1) )

            do i=2,ido
                ch(1,:,1,i) = cc(1,:,i,1)+cc(1,:,i,2)
                tr2 = cc(1,:,i,1)-cc(1,:,i,2)
                ch(2,:,1,i) = cc(2,:,i,1)+cc(2,:,i,2)
                ti2 = cc(2,:,i,1)-cc(2,:,i,2)
                ch(2,:,2,i) = wa(i,1,1)*ti2+wa(i,1,2)*tr2
                ch(1,:,2,i) = wa(i,1,1)*tr2-wa(i,1,2)*ti2
            end do
            !
            !==> Release memory
            !
            deallocate( tr2 )
            deallocate( ti2 )
        end if

    end subroutine c1f2kb



    subroutine c1f3kb(ido, l1, na, cc, in1, ch, in2, wa)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        integer (ip), intent (in)     :: na
        real (wp),    intent (in out) :: cc(in1,l1,ido,3)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,l1,3,ido)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa(ido,2,2)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip)           :: i !! Counter
        real (wp), allocatable :: di2(:), di3(:)
        real (wp), allocatable :: dr2(:), dr3(:)
        real (wp), allocatable :: ci2(:), ci3(:)
        real (wp), allocatable :: cr2(:), cr3(:)
        real (wp), allocatable :: ti2(:), tr2(:)
        real (wp), parameter   :: TAUI = sqrt(THREE)/2!0.866025403784439_wp
        real (wp), parameter   :: TAUR = -HALF
        !----------------------------------------------------------------------

        !
        !==> Allocate memory
        !
        allocate( ci2(l1), ci3(l1) )
        allocate( cr2(l1), cr3(l1) )
        allocate( ti2(l1), tr2(l1) )

        if (1 >= ido .and. na /= 1) then
            tr2 = cc(1,:,1,2)+cc(1,:,1,3)
            cr2 = cc(1,:,1,1)+TAUR*tr2
            cc(1,:,1,1) = cc(1,:,1,1)+tr2
            ti2 = cc(2,:,1,2)+cc(2,:,1,3)
            ci2 = cc(2,:,1,1)+TAUR*ti2
            cc(2,:,1,1) = cc(2,:,1,1)+ti2
            cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
            ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
            cc(1,:,1,2) = cr2-ci3
            cc(1,:,1,3) = cr2+ci3
            cc(2,:,1,2) = ci2+cr3
            cc(2,:,1,3) = ci2-cr3
        else
            tr2 = cc(1,:,1,2)+cc(1,:,1,3)
            cr2 = cc(1,:,1,1)+TAUR*tr2
            ch(1,:,1,1) = cc(1,:,1,1)+tr2
            ti2 = cc(2,:,1,2)+cc(2,:,1,3)
            ci2 = cc(2,:,1,1)+TAUR*ti2
            ch(2,:,1,1) = cc(2,:,1,1)+ti2
            cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
            ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
            ch(1,:,2,1) = cr2-ci3
            ch(1,:,3,1) = cr2+ci3
            ch(2,:,2,1) = ci2+cr3
            ch(2,:,3,1) = ci2-cr3

            !
            !==> Allocate memory
            !
            allocate( dr2(l1), dr3(l1) )
            allocate( di2(l1), di3(l1) )

            do i=2,ido
                tr2 = cc(1,:,i,2)+cc(1,:,i,3)
                cr2 = cc(1,:,i,1)+TAUR*tr2
                ch(1,:,1,i) = cc(1,:,i,1)+tr2
                ti2 = cc(2,:,i,2)+cc(2,:,i,3)
                ci2 = cc(2,:,i,1)+TAUR*ti2
                ch(2,:,1,i) = cc(2,:,i,1)+ti2
                cr3 = TAUI*(cc(1,:,i,2)-cc(1,:,i,3))
                ci3 = TAUI*(cc(2,:,i,2)-cc(2,:,i,3))

                dr2 = cr2-ci3
                dr3 = cr2+ci3
                di2 = ci2+cr3
                di3 = ci2-cr3

                ch(2,:,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
                ch(1,:,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
                ch(2,:,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
                ch(1,:,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3

            end do
            !
            !==> Release memory
            !
            deallocate( dr2, dr3 )
            deallocate( di2, di3 )
        end if

        !
        !==> Release memory
        !
        deallocate( ci2, ci3 )
        deallocate( cr2, cr3 )
        deallocate( ti2, tr2 )

    end subroutine c1f3kb


    subroutine c1f4kb(ido, l1, na, cc, in1, ch, in2, wa)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        integer (ip), intent (in)     :: na
        real (wp),    intent (in out) :: cc(in1,l1,ido,4)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,l1,4,ido)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa(ido,3,2)
        !----------------------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------------------
        integer (ip)           :: i !! Counter
        real (wp), allocatable :: ci2(:), ci3(:), ci4(:)
        real (wp), allocatable :: cr2(:), cr3(:), cr4(:)
        real (wp), allocatable :: ti1(:), ti2(:), ti3(:), ti4(:)
        real (wp), allocatable :: tr1(:), tr2(:), tr3(:), tr4(:)
        !----------------------------------------------------------------------

        !
        !==> Allocate memory
        !
        allocate( ti1(l1), ti2(l1), ti3(l1), ti4(l1) )
        allocate( tr1(l1), tr2(l1), tr3(l1), tr4(l1) )

        if (1 >= ido .and. na /= 1) then
            ti1 = cc(2,:,1,1)-cc(2,:,1,3)
            ti2 = cc(2,:,1,1)+cc(2,:,1,3)
            tr4 = cc(2,:,1,4)-cc(2,:,1,2)
            ti3 = cc(2,:,1,2)+cc(2,:,1,4)
            tr1 = cc(1,:,1,1)-cc(1,:,1,3)
            tr2 = cc(1,:,1,1)+cc(1,:,1,3)
            ti4 = cc(1,:,1,2)-cc(1,:,1,4)
            tr3 = cc(1,:,1,2)+cc(1,:,1,4)
            cc(1,:,1,1) = tr2+tr3
            cc(1,:,1,3) = tr2-tr3
            cc(2,:,1,1) = ti2+ti3
            cc(2,:,1,3) = ti2-ti3
            cc(1,:,1,2) = tr1+tr4
            cc(1,:,1,4) = tr1-tr4
            cc(2,:,1,2) = ti1+ti4
            cc(2,:,1,4) = ti1-ti4
        else
            ti1 = cc(2,:,1,1)-cc(2,:,1,3)
            ti2 = cc(2,:,1,1)+cc(2,:,1,3)
            tr4 = cc(2,:,1,4)-cc(2,:,1,2)
            ti3 = cc(2,:,1,2)+cc(2,:,1,4)
            tr1 = cc(1,:,1,1)-cc(1,:,1,3)
            tr2 = cc(1,:,1,1)+cc(1,:,1,3)
            ti4 = cc(1,:,1,2)-cc(1,:,1,4)
            tr3 = cc(1,:,1,2)+cc(1,:,1,4)
            ch(1,:,1,1) = tr2+tr3
            ch(1,:,3,1) = tr2-tr3
            ch(2,:,1,1) = ti2+ti3
            ch(2,:,3,1) = ti2-ti3
            ch(1,:,2,1) = tr1+tr4
            ch(1,:,4,1) = tr1-tr4
            ch(2,:,2,1) = ti1+ti4
            ch(2,:,4,1) = ti1-ti4

            !
            !==> Allocate memory
            !
            allocate( ci2(l1), ci3(l1), ci4(l1) )
            allocate( cr2(l1), cr3(l1), cr4(l1) )

            do i=2,ido
                ti1 = cc(2,:,i,1)-cc(2,:,i,3)
                ti2 = cc(2,:,i,1)+cc(2,:,i,3)
                ti3 = cc(2,:,i,2)+cc(2,:,i,4)
                tr4 = cc(2,:,i,4)-cc(2,:,i,2)
                tr1 = cc(1,:,i,1)-cc(1,:,i,3)
                tr2 = cc(1,:,i,1)+cc(1,:,i,3)
                ti4 = cc(1,:,i,2)-cc(1,:,i,4)
                tr3 = cc(1,:,i,2)+cc(1,:,i,4)
                ch(1,:,1,i) = tr2+tr3
                cr3 = tr2-tr3
                ch(2,:,1,i) = ti2+ti3
                ci3 = ti2-ti3
                cr2 = tr1+tr4
                cr4 = tr1-tr4
                ci2 = ti1+ti4
                ci4 = ti1-ti4
                ch(1,:,2,i) = wa(i,1,1)*cr2-wa(i,1,2)*ci2
                ch(2,:,2,i) = wa(i,1,1)*ci2+wa(i,1,2)*cr2
                ch(1,:,3,i) = wa(i,2,1)*cr3-wa(i,2,2)*ci3
                ch(2,:,3,i) = wa(i,2,1)*ci3+wa(i,2,2)*cr3
                ch(1,:,4,i) = wa(i,3,1)*cr4-wa(i,3,2)*ci4
                ch(2,:,4,i) = wa(i,3,1)*ci4+wa(i,3,2)*cr4
            end do
            !
            !==> Release memory
            !
            deallocate( ci2, ci3, ci4 )
            deallocate( cr2, cr3, cr4 )
        end if

        !
        !==> Release memory
        !
        deallocate( ti1, ti2, ti3, ti4 )
        deallocate( tr1, tr2, tr3, tr4 )

    end subroutine c1f4kb

    subroutine c1f5kb(ido, l1, na, cc, in1, ch, in2, wa)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: l1
        integer (ip), intent (in)     :: na
        real (wp),    intent (in out) :: cc(in1,l1,ido,5)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,l1,5,ido)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa(ido,4,2)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip)           :: i !! Counter
        real (wp), allocatable :: chold1(:), chold2(:)
        real (wp), allocatable :: ci2(:), ci3(:), ci4(:), ci5(:)
        real (wp), allocatable :: cr2(:), cr3(:), cr4(:), cr5(:)
        real (wp), allocatable :: di2(:), di3(:), di4(:), di5(:)
        real (wp), allocatable :: dr2(:), dr3(:), dr4(:), dr5(:)
        real (wp), allocatable :: ti2(:), ti3(:), ti4(:), ti5(:)
        real (wp), allocatable :: tr2(:), tr3(:), tr4(:), tr5(:)
        real (wp), parameter   :: SQRT5 = sqrt(FIVE)
        real (wp), parameter   :: SQRT5_PLUS_5 = SQRT5 + FIVE
        real (wp), parameter   :: TI11 = sqrt(SQRT5_PLUS_5/2)/2             ! 0.9510565162951536_wp
        real (wp), parameter   :: TI12 = sqrt(FIVE/(TWO*SQRT5_PLUS_5)) ! 0.5877852522924731_wp
        real (wp), parameter   :: TR11 =  (SQRT5 - ONE)/4                 ! 0.3090169943749474_wp
        real (wp), parameter   :: TR12 = -(ONE + SQRT5)/4                 !-0.8090169943749474_wp
        !------------------------------------------------------------------

        !
        !==> Allocate memory
        !
        allocate( chold1(l1), chold2(l1) )
        allocate( ti2(l1), ti3(l1), ti4(l1), ti5(l1) )
        allocate( tr2(l1),  tr3(l1), tr4(l1), tr5(l1) )
        allocate( cr2(l1), cr3(l1), cr4(l1), cr5(l1) )
        allocate( ci2(l1), ci3(l1), ci4(l1), ci5(l1) )

        if (1 >= ido .and. na /= 1) then
            ti5 = cc(2,:,1,2)-cc(2,:,1,5)
            ti2 = cc(2,:,1,2)+cc(2,:,1,5)
            ti4 = cc(2,:,1,3)-cc(2,:,1,4)
            ti3 = cc(2,:,1,3)+cc(2,:,1,4)
            tr5 = cc(1,:,1,2)-cc(1,:,1,5)
            tr2 = cc(1,:,1,2)+cc(1,:,1,5)
            tr4 = cc(1,:,1,3)-cc(1,:,1,4)
            tr3 = cc(1,:,1,3)+cc(1,:,1,4)
            chold1 = cc(1,:,1,1)+tr2+tr3
            chold2 = cc(2,:,1,1)+ti2+ti3
            cr2 = cc(1,:,1,1)+TR11*tr2+TR12*tr3
            ci2 = cc(2,:,1,1)+TR11*ti2+TR12*ti3
            cr3 = cc(1,:,1,1)+TR12*tr2+TR11*tr3
            ci3 = cc(2,:,1,1)+TR12*ti2+TR11*ti3
            cc(1,:,1,1) = chold1
            cc(2,:,1,1) = chold2
            cr5 = TI11*tr5+TI12*tr4
            ci5 = TI11*ti5+TI12*ti4
            cr4 = TI12*tr5-TI11*tr4
            ci4 = TI12*ti5-TI11*ti4
            cc(1,:,1,2) = cr2-ci5
            cc(1,:,1,5) = cr2+ci5
            cc(2,:,1,2) = ci2+cr5
            cc(2,:,1,3) = ci3+cr4
            cc(1,:,1,3) = cr3-ci4
            cc(1,:,1,4) = cr3+ci4
            cc(2,:,1,4) = ci3-cr4
            cc(2,:,1,5) = ci2-cr5
        else
            ti5 = cc(2,:,1,2)-cc(2,:,1,5)
            ti2 = cc(2,:,1,2)+cc(2,:,1,5)
            ti4 = cc(2,:,1,3)-cc(2,:,1,4)
            ti3 = cc(2,:,1,3)+cc(2,:,1,4)
            tr5 = cc(1,:,1,2)-cc(1,:,1,5)
            tr2 = cc(1,:,1,2)+cc(1,:,1,5)
            tr4 = cc(1,:,1,3)-cc(1,:,1,4)
            tr3 = cc(1,:,1,3)+cc(1,:,1,4)
            ch(1,:,1,1) = cc(1,:,1,1)+tr2+tr3
            ch(2,:,1,1) = cc(2,:,1,1)+ti2+ti3
            cr2 = cc(1,:,1,1)+TR11*tr2+TR12*tr3
            ci2 = cc(2,:,1,1)+TR11*ti2+TR12*ti3
            cr3 = cc(1,:,1,1)+TR12*tr2+TR11*tr3
            ci3 = cc(2,:,1,1)+TR12*ti2+TR11*ti3
            cr5 = TI11*tr5+TI12*tr4
            ci5 = TI11*ti5+TI12*ti4
            cr4 = TI12*tr5-TI11*tr4
            ci4 = TI12*ti5-TI11*ti4
            ch(1,:,2,1) = cr2-ci5
            ch(1,:,5,1) = cr2+ci5
            ch(2,:,2,1) = ci2+cr5
            ch(2,:,3,1) = ci3+cr4
            ch(1,:,3,1) = cr3-ci4
            ch(1,:,4,1) = cr3+ci4
            ch(2,:,4,1) = ci3-cr4
            ch(2,:,5,1) = ci2-cr5

            !
            !==> Allocate memory
            !
            allocate( di2(l1), di3(l1), di4(l1), di5(l1) )
            allocate( dr2(l1), dr3(l1), dr4(l1), dr5(l1) )

            do i=2,ido
                ti5 = cc(2,:,i,2)-cc(2,:,i,5)
                ti2 = cc(2,:,i,2)+cc(2,:,i,5)
                ti4 = cc(2,:,i,3)-cc(2,:,i,4)
                ti3 = cc(2,:,i,3)+cc(2,:,i,4)
                tr5 = cc(1,:,i,2)-cc(1,:,i,5)
                tr2 = cc(1,:,i,2)+cc(1,:,i,5)
                tr4 = cc(1,:,i,3)-cc(1,:,i,4)
                tr3 = cc(1,:,i,3)+cc(1,:,i,4)
                ch(1,:,1,i) = cc(1,:,i,1)+tr2+tr3
                ch(2,:,1,i) = cc(2,:,i,1)+ti2+ti3
                cr2 = cc(1,:,i,1)+TR11*tr2+TR12*tr3
                ci2 = cc(2,:,i,1)+TR11*ti2+TR12*ti3
                cr3 = cc(1,:,i,1)+TR12*tr2+TR11*tr3
                ci3 = cc(2,:,i,1)+TR12*ti2+TR11*ti3
                cr5 = TI11*tr5+TI12*tr4
                ci5 = TI11*ti5+TI12*ti4
                cr4 = TI12*tr5-TI11*tr4
                ci4 = TI12*ti5-TI11*ti4
                dr3 = cr3-ci4
                dr4 = cr3+ci4
                di3 = ci3+cr4
                di4 = ci3-cr4
                dr5 = cr2+ci5
                dr2 = cr2-ci5
                di5 = ci2-cr5
                di2 = ci2+cr5
                ch(1,:,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
                ch(2,:,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
                ch(1,:,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
                ch(2,:,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
                ch(1,:,4,i) = wa(i,3,1)*dr4-wa(i,3,2)*di4
                ch(2,:,4,i) = wa(i,3,1)*di4+wa(i,3,2)*dr4
                ch(1,:,5,i) = wa(i,4,1)*dr5-wa(i,4,2)*di5
                ch(2,:,5,i) = wa(i,4,1)*di5+wa(i,4,2)*dr5
            end do
            !
            !==> Release memory
            !
            deallocate( di2, di3, di4, di5 )
            deallocate( dr2, dr3, dr4, dr5 )
        end if

        !
        !==> Release memory
        !
        deallocate( ti2, ti3, ti4, ti5 )
        deallocate( tr2,  tr3, tr4, tr5 )
        deallocate( chold1, chold2 )
        deallocate( cr2, cr3, cr4, cr5 )
        deallocate( ci2, ci3, ci4, ci5 )

    end subroutine c1f5kb


    subroutine c1fgkb(ido, iip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa)
        !------------------------------------------------------------------
        ! Dummy arguments
        !------------------------------------------------------------------
        integer (ip), intent (in)     :: ido
        integer (ip), intent (in)     :: iip
        integer (ip), intent (in)     :: l1
        integer (ip), intent (in)     :: lid
        integer (ip), intent (in)     :: na
        real (wp),    intent (in out) :: cc(in1,l1,iip,ido)
        real (wp),    intent (in out) :: cc1(in1,lid,iip)
        integer (ip), intent (in)     :: in1
        real (wp),    intent (in out) :: ch(in2,l1,ido,iip)
        real (wp),    intent (in out) :: ch1(in2,lid,iip)
        integer (ip), intent (in)     :: in2
        real (wp),    intent (in)     :: wa(ido,iip-1,2)
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        integer (ip)           :: i, idlj, iipp2, iipph
        integer (ip)           :: j, jc, l, lc
        real (wp)              :: wai, war
        real (wp), allocatable :: chold1(:), chold2(:)
        !------------------------------------------------------------------

        iipp2 = iip+2
        iipph = (iip+1)/2

        ch1(1,:,1) = cc1(1,:,1)
        ch1(2,:,1) = cc1(2,:,1)

        do j=2,iipph
            jc = iipp2-j
            ch1(1,:,j) =  cc1(1,:,j)+cc1(1,:,jc)
            ch1(1,:,jc) = cc1(1,:,j)-cc1(1,:,jc)
            ch1(2,:,j) =  cc1(2,:,j)+cc1(2,:,jc)
            ch1(2,:,jc) = cc1(2,:,j)-cc1(2,:,jc)
        end do

        do j=2,iipph
            cc1(1,:,1) = cc1(1,:,1)+ch1(1,:,j)
            cc1(2,:,1) = cc1(2,:,1)+ch1(2,:,j)
        end do

        do l=2,iipph

            lc = iipp2-l
            cc1(1,:,l) = ch1(1,:,1)+wa(1,l-1,1)*ch1(1,:,2)
            cc1(1,:,lc) = wa(1,l-1,2)*ch1(1,:,iip)
            cc1(2,:,l) = ch1(2,:,1)+wa(1,l-1,1)*ch1(2,:,2)
            cc1(2,:,lc) = wa(1,l-1,2)*ch1(2,:,iip)

            do j=3,iipph
                jc = iipp2-j
                idlj = mod((l-1)*(j-1),iip)
                war = wa(1,idlj,1)
                wai = wa(1,idlj,2)
                cc1(1,:,l) = cc1(1,:,l)+war*ch1(1,:,j)
                cc1(1,:,lc) = cc1(1,:,lc)+wai*ch1(1,:,jc)
                cc1(2,:,l) = cc1(2,:,l)+war*ch1(2,:,j)
                cc1(2,:,lc) = cc1(2,:,lc)+wai*ch1(2,:,jc)
            end do
        end do

        if (1 >= ido .and. na /= 1) then
            !
            !==> Allocate memory
            !
            allocate( chold1(lid), chold2(lid) )

            do j=2,iipph
                jc = iipp2-j
                chold1 = cc1(1,:,j)-cc1(2,:,jc)
                chold2 = cc1(1,:,j)+cc1(2,:,jc)
                cc1(1,:,j) = chold1
                cc1(2,:,jc) = cc1(2,:,j)-cc1(1,:,jc)
                cc1(2,:,j) = cc1(2,:,j)+cc1(1,:,jc)
                cc1(1,:,jc) = chold2
            end do
            !
            !==> Release memory
            !
            deallocate( chold1, chold2 )
        else
            ch1(1,:,1) = cc1(1,:,1)
            ch1(2,:,1) = cc1(2,:,1)

            do j=2,iipph
                jc = iipp2-j
                ch1(1,:,j) = cc1(1,:,j)-cc1(2,:,jc)
                ch1(1,:,jc) = cc1(1,:,j)+cc1(2,:,jc)
                ch1(2,:,jc) = cc1(2,:,j)-cc1(1,:,jc)
                ch1(2,:,j) = cc1(2,:,j)+cc1(1,:,jc)
            end do

            if (ido /= 1) then
                cc(1,:,1,:) = ch(1,:,:,1)
                cc(2,:,1,:) = ch(2,:,:,1)

                cc(1,:,2:iip,1) = ch(1,:,1,2:iip)
                cc(2,:,2:iip,1) = ch(2,:,1,2:iip)

                do j=2,iip
                    do i=2,ido
                        cc(1,:,j,i) = wa(i,j-1,1)*ch(1,:,i,j) &
                            -wa(i,j-1,2)*ch(2,:,i,j)
                        cc(2,:,j,i) = wa(i,j-1,1)*ch(2,:,i,j) &
                            +wa(i,j-1,2)*ch(1,:,i,j)
                    end do
                end do
            end if
        end if

    end subroutine c1fgkb


end submodule complex_backward_1d
