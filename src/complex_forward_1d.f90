submodule (complex_transform_routines) complex_forward_1d

contains

    module subroutine cfft1f(n, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
        !
        ! cfft1f: complex 64-bit precision forward fast fourier transform, 1d.
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
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: inc
        complex (wp), intent (in out) :: c(lenc)
        integer (ip), intent (in)     :: lenc
        real (wp),    intent (in)     :: wsave(lensav)
        integer (ip), intent (in)     :: lensav
        real (wp),    intent (out)    :: work(lenwrk)
        integer (ip), intent (in)     :: lenwrk
        integer (ip), intent (out)    :: ier
        !------------------------------------------------------------------
        ! Local variables
        !------------------------------------------------------------------
        real (wp), allocatable :: real_copy(:,:)
        integer (ip)           :: iw1
        !------------------------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (lenc < inc*(n-1) + 1) then
            ier = 1
            call fft_error_handler('cfft1f ', 4)
        else if (size(wsave) < get_complex_1d_saved_workspace_length(n)) then
            ier = 2
            call fft_error_handler('cfft1f ', 6)
        else if (size(work) < get_complex_1d_workspace_length(n)) then
            ier = 3
            call fft_error_handler('cfft1f ', 8)
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

            ! Make copy
            real_copy(1,:) = real(c)
            real_copy(2,:) = aimag(c)

            ! Set workspace index pointer
            iw1 = 2 * n+1

            call c1fm1f(n,inc,real_copy,work,wsave,wsave(iw1),wsave(iw1+1:))

            ! Make copy
            c =  cmplx(real_copy(1,:), real_copy(2,:), kind=wp)

            !
            !==> Release memory
            !
            deallocate( real_copy )
        end if

    end subroutine cfft1f



    subroutine c1fm1f(n, inc, c, ch, wa, fnf, fac)
        !----------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------
        integer (ip), intent (in)     :: n
        integer (ip), intent (in)     :: inc
        real (wp),    intent (in out) :: c(:,:)
        real (wp),    intent (out)    :: ch(:)
        real (wp),    intent (in)     :: wa(*)
        real (wp),    intent (in)     :: fnf
        real (wp),    intent (in)     :: fac(:)
        !----------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------
        integer (ip) :: ido, inc2, iip, iw
        integer (ip) :: k1, l1, l2, lid
        integer (ip) :: na, nbr, nf
        !----------------------------------------------------------

        inc2 = 2*inc
        nf = int(fnf, kind=ip)
        na = 0
        l1 = 1
        iw = 1

        do k1=1,nf
            iip = int(fac(k1), kind=ip)
            l2 = iip*l1
            ido = n/l2
            lid = l1*ido
            nbr = 1+na+2*min(iip-2,4)
            select case (nbr)
                case (1)
                    call c1f2kf(ido,l1,na,c,inc2,ch,2,wa(iw))
                case (2)
                    call c1f2kf(ido,l1,na,ch,2,c,inc2,wa(iw))
                case (3)
                    call c1f3kf(ido,l1,na,c,inc2,ch,2,wa(iw))
                case (4)
                    call c1f3kf(ido,l1,na,ch,2,c,inc2,wa(iw))
                case (5)
                    call c1f4kf(ido,l1,na,c,inc2,ch,2,wa(iw))
                case (6)
                    call c1f4kf(ido,l1,na,ch,2,c,inc2,wa(iw))
                case (7)
                    call c1f5kf(ido,l1,na,c,inc2,ch,2,wa(iw))
                case (8)
                    call c1f5kf(ido,l1,na,ch,2,c,inc2,wa(iw))
                case (9)
                    call c1fgkf(ido,iip,l1,lid,na,c,c,inc2,ch,ch,2,wa(iw))
                case (10)
                    call c1fgkf(ido,iip,l1,lid,na,ch,ch,2,c,c,inc2,wa(iw))
            end select

            l1 = l2
            iw = iw+(iip-1)*(2*ido)
            if (iip <= 5) then
                na = 1-na
            end if
        end do

    end subroutine c1fm1f


    subroutine c1f2kf(ido, l1, na, cc, in1, ch, in2, wa)
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
        integer (ip)           :: i !! counter
        real (wp)              :: sn
        real (wp), allocatable :: chold1(:), chold2(:)
        real (wp), allocatable :: ti2(:), tr2(:)
        !----------------------------------------------------------------------

        if (1 >= ido) then
            sn = ONE/(2 * l1)
            if (na /= 1) then
                !
                !==> Allocate memory
                !
                allocate( chold1(l1) )
                allocate( chold2(l1) )

                chold1 = sn*(cc(1,:,1,1)+cc(1,:,1,2))
                cc(1,:,1,2) = sn*(cc(1,:,1,1)-cc(1,:,1,2))
                cc(1,:,1,1) = chold1
                chold2 = sn*(cc(2,:,1,1)+cc(2,:,1,2))
                cc(2,:,1,2) = sn*(cc(2,:,1,1)-cc(2,:,1,2))
                cc(2,:,1,1) = chold2
                !
                !==> Release memory
                !
                deallocate( chold1 )
                deallocate( chold2 )
            else
                ch(1,:,1,1) = sn*(cc(1,:,1,1)+cc(1,:,1,2))
                ch(1,:,2,1) = sn*(cc(1,:,1,1)-cc(1,:,1,2))
                ch(2,:,1,1) = sn*(cc(2,:,1,1)+cc(2,:,1,2))
                ch(2,:,2,1) = sn*(cc(2,:,1,1)-cc(2,:,1,2))
            end if
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
                ch(2,:,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
                ch(1,:,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
            end do
            !
            !==> Release memory
            !
            deallocate( tr2 )
            deallocate( ti2 )
        end if

    end subroutine c1f2kf


    subroutine c1f3kf(ido, l1, na, cc, in1, ch, in2, wa)
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
        real (wp), allocatable :: ci2(:), ci3(:)
        real (wp), allocatable :: cr2(:), cr3(:)
        real (wp), allocatable :: di2(:), di3(:)
        real (wp), allocatable :: dr2(:), dr3(:)
        real (wp), allocatable :: ti2(:), tr2(:)
        real (wp), parameter   :: TAUI = -sqrt(THREE)/2!-0.866025403784439_wp
        real (wp), parameter   :: TAUR = -HALF
        real (wp)              :: sn
        !----------------------------------------------------------------------

        !
        !==> Allocate memory
        !
        allocate( ci2(l1), ci3(l1) )
        allocate( cr2(l1), cr3(l1) )
        allocate( ti2(l1), tr2(l1) )

        if (1 >= ido) then
            sn = ONE/(3 * l1)
            if (na /= 1) then
                tr2 = cc(1,:,1,2)+cc(1,:,1,3)
                cr2 = cc(1,:,1,1)+TAUR*tr2
                cc(1,:,1,1) = sn*(cc(1,:,1,1)+tr2)
                ti2 = cc(2,:,1,2)+cc(2,:,1,3)
                ci2 = cc(2,:,1,1)+TAUR*ti2
                cc(2,:,1,1) = sn*(cc(2,:,1,1)+ti2)
                cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
                ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
                cc(1,:,1,2) = sn*(cr2-ci3)
                cc(1,:,1,3) = sn*(cr2+ci3)
                cc(2,:,1,2) = sn*(ci2+cr3)
                cc(2,:,1,3) = sn*(ci2-cr3)
            else
                tr2 = cc(1,:,1,2)+cc(1,:,1,3)
                cr2 = cc(1,:,1,1)+TAUR*tr2
                ch(1,:,1,1) = sn*(cc(1,:,1,1)+tr2)
                ti2 = cc(2,:,1,2)+cc(2,:,1,3)
                ci2 = cc(2,:,1,1)+TAUR*ti2
                ch(2,:,1,1) = sn*(cc(2,:,1,1)+ti2)
                cr3 = TAUI*(cc(1,:,1,2)-cc(1,:,1,3))
                ci3 = TAUI*(cc(2,:,1,2)-cc(2,:,1,3))
                ch(1,:,2,1) = sn*(cr2-ci3)
                ch(1,:,3,1) = sn*(cr2+ci3)
                ch(2,:,2,1) = sn*(ci2+cr3)
                ch(2,:,3,1) = sn*(ci2-cr3)
            end if
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
            allocate( di2(l1), di3(l1) )
            allocate( dr2(l1), dr3(l1) )

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
                ch(2,:,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
                ch(1,:,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
                ch(2,:,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
                ch(1,:,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
            end do
            !
            !==> Release memory
            !
            deallocate( di2, di3 )
            deallocate( dr2, dr3 )
        end if

        !
        !==> Release memory
        !
        deallocate( ci2, ci3 )
        deallocate( cr2, cr3 )
        deallocate( ti2, tr2 )

    end subroutine c1f3kf


    subroutine c1f4kf(ido, l1, na, cc, in1, ch, in2, wa)
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
        integer (ip) :: i
        real (wp)    :: sn
        !----------------------------------------------------------------------

        if (1 >= ido) then
            sn = ONE/(4 * l1)
            if (na /= 1) then
                associate(&
                    ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
                    ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
                    tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
                    ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
                    tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
                    tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
                    ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
                    tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
                    )
                    cc(1,:,1,1) = sn*(tr2+tr3)
                    cc(1,:,1,3) = sn*(tr2-tr3)
                    cc(2,:,1,1) = sn*(ti2+ti3)
                    cc(2,:,1,3) = sn*(ti2-ti3)
                    cc(1,:,1,2) = sn*(tr1+tr4)
                    cc(1,:,1,4) = sn*(tr1-tr4)
                    cc(2,:,1,2) = sn*(ti1+ti4)
                    cc(2,:,1,4) = sn*(ti1-ti4)
                end associate
            else
                associate( &
                    ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
                    ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
                    tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
                    ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
                    tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
                    tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
                    ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
                    tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
                    )
                    ch(1,:,1,1) = sn*(tr2+tr3)
                    ch(1,:,3,1) = sn*(tr2-tr3)
                    ch(2,:,1,1) = sn*(ti2+ti3)
                    ch(2,:,3,1) = sn*(ti2-ti3)
                    ch(1,:,2,1) = sn*(tr1+tr4)
                    ch(1,:,4,1) = sn*(tr1-tr4)
                    ch(2,:,2,1) = sn*(ti1+ti4)
                    ch(2,:,4,1) = sn*(ti1-ti4)
                end associate
            end if
        else
            associate( &
                ti1 => cc(2,:,1,1)-cc(2,:,1,3), &
                ti2 => cc(2,:,1,1)+cc(2,:,1,3), &
                tr4 => cc(2,:,1,2)-cc(2,:,1,4), &
                ti3 => cc(2,:,1,2)+cc(2,:,1,4), &
                tr1 => cc(1,:,1,1)-cc(1,:,1,3), &
                tr2 => cc(1,:,1,1)+cc(1,:,1,3), &
                ti4 => cc(1,:,1,4)-cc(1,:,1,2), &
                tr3 => cc(1,:,1,2)+cc(1,:,1,4) &
                )
                ch(1,:,1,1) = tr2+tr3
                ch(1,:,3,1) = tr2-tr3
                ch(2,:,1,1) = ti2+ti3
                ch(2,:,3,1) = ti2-ti3
                ch(1,:,2,1) = tr1+tr4
                ch(1,:,4,1) = tr1-tr4
                ch(2,:,2,1) = ti1+ti4
                ch(2,:,4,1) = ti1-ti4
            end associate
            do i=2,ido
                associate( &
                    ti1 => cc(2,:,i,1)-cc(2,:,i,3), &
                    ti2 => cc(2,:,i,1)+cc(2,:,i,3), &
                    ti3 => cc(2,:,i,2)+cc(2,:,i,4), &
                    tr4 => cc(2,:,i,2)-cc(2,:,i,4), &
                    tr1 => cc(1,:,i,1)-cc(1,:,i,3), &
                    tr2 => cc(1,:,i,1)+cc(1,:,i,3), &
                    ti4 => cc(1,:,i,4)-cc(1,:,i,2), &
                    tr3 => cc(1,:,i,2)+cc(1,:,i,4) &
                    )
                    ch(1,:,1,i) = tr2+tr3
                    associate( cr3 => tr2-tr3 )
                        ch(2,:,1,i) = ti2+ti3
                        associate( &
                            ci3 => ti2-ti3, &
                            cr2 => tr1+tr4, &
                            cr4 => tr1-tr4, &
                            ci2 => ti1+ti4, &
                            ci4 => ti1-ti4 &
                            )
                            ch(1,:,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
                            ch(2,:,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
                            ch(1,:,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
                            ch(2,:,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
                            ch(1,:,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
                            ch(2,:,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
                        end associate
                    end associate
                end associate
            end do
        end if

    end subroutine c1f4kf


    subroutine c1f5kf(ido, l1, na, cc, in1, ch, in2, wa)
        !----------------------------------------------------------
        ! Dummy arguments
        !----------------------------------------------------------
        integer (ip) ido
        integer (ip) in1
        integer (ip) in2
        integer (ip) l1
        integer (ip) na
        real (wp) cc(in1,l1,ido,5)
        real (wp) ch(in2,l1,5,ido)
        real (wp) wa(ido,4,2)
        !----------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------
        integer (ip)         :: i, k
        real (wp)            :: sn
        real (wp)            :: chold1, chold2
        real (wp)            :: ci2, ci3, ci4, ci5
        real (wp)            :: cr2, cr3, cr4, cr5
        real (wp)            :: di2, di3, di4, di5
        real (wp)            :: dr2, dr3, dr4, dr5
        real (wp)            :: ti2, ti3, ti4, ti5
        real (wp)            :: tr2, tr3, tr4, tr5
        real (wp), parameter :: SQRT5 = sqrt(FIVE)
        real (wp), parameter :: SQRT5_PLUS_5 = SQRT5 + FIVE
        real (wp), parameter :: TI11 = -sqrt(SQRT5_PLUS_5/2)/2             !-0.9510565162951536_wp
        real (wp), parameter :: TI12 = -sqrt(FIVE/(TWO*SQRT5_PLUS_5)) !-0.5877852522924731_wp
        real (wp), parameter :: TR11 =  (SQRT5 - ONE)/4                 ! 0.3090169943749474_wp
        real (wp), parameter :: TR12 = -(ONE + SQRT5)/4                 !-0.8090169943749474_wp
        !----------------------------------------------------------

        if (1 >= ido) then
            sn = ONE/(5 * l1)
            if (na /= 1) then
                do k=1,l1
                    ti5 = cc(2,k,1,2)-cc(2,k,1,5)
                    ti2 = cc(2,k,1,2)+cc(2,k,1,5)
                    ti4 = cc(2,k,1,3)-cc(2,k,1,4)
                    ti3 = cc(2,k,1,3)+cc(2,k,1,4)
                    tr5 = cc(1,k,1,2)-cc(1,k,1,5)
                    tr2 = cc(1,k,1,2)+cc(1,k,1,5)
                    tr4 = cc(1,k,1,3)-cc(1,k,1,4)
                    tr3 = cc(1,k,1,3)+cc(1,k,1,4)
                    chold1 = sn*(cc(1,k,1,1)+tr2+tr3)
                    chold2 = sn*(cc(2,k,1,1)+ti2+ti3)
                    cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
                    ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
                    cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
                    ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
                    cc(1,k,1,1) = chold1
                    cc(2,k,1,1) = chold2
                    cr5 = TI11*tr5+TI12*tr4
                    ci5 = TI11*ti5+TI12*ti4
                    cr4 = TI12*tr5-TI11*tr4
                    ci4 = TI12*ti5-TI11*ti4
                    cc(1,k,1,2) = sn*(cr2-ci5)
                    cc(1,k,1,5) = sn*(cr2+ci5)
                    cc(2,k,1,2) = sn*(ci2+cr5)
                    cc(2,k,1,3) = sn*(ci3+cr4)
                    cc(1,k,1,3) = sn*(cr3-ci4)
                    cc(1,k,1,4) = sn*(cr3+ci4)
                    cc(2,k,1,4) = sn*(ci3-cr4)
                    cc(2,k,1,5) = sn*(ci2-cr5)
                end do
            else
                do k=1,l1
                    ti5 = cc(2,k,1,2)-cc(2,k,1,5)
                    ti2 = cc(2,k,1,2)+cc(2,k,1,5)
                    ti4 = cc(2,k,1,3)-cc(2,k,1,4)
                    ti3 = cc(2,k,1,3)+cc(2,k,1,4)
                    tr5 = cc(1,k,1,2)-cc(1,k,1,5)
                    tr2 = cc(1,k,1,2)+cc(1,k,1,5)
                    tr4 = cc(1,k,1,3)-cc(1,k,1,4)
                    tr3 = cc(1,k,1,3)+cc(1,k,1,4)
                    ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2+tr3)
                    ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2+ti3)
                    cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
                    ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
                    cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
                    ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
                    cr5 = TI11*tr5+TI12*tr4
                    ci5 = TI11*ti5+TI12*ti4
                    cr4 = TI12*tr5-TI11*tr4
                    ci4 = TI12*ti5-TI11*ti4
                    ch(1,k,2,1) = sn*(cr2-ci5)
                    ch(1,k,5,1) = sn*(cr2+ci5)
                    ch(2,k,2,1) = sn*(ci2+cr5)
                    ch(2,k,3,1) = sn*(ci3+cr4)
                    ch(1,k,3,1) = sn*(cr3-ci4)
                    ch(1,k,4,1) = sn*(cr3+ci4)
                    ch(2,k,4,1) = sn*(ci3-cr4)
                    ch(2,k,5,1) = sn*(ci2-cr5)
                end do
            end if
        else
            do k=1,l1
                ti5 = cc(2,k,1,2)-cc(2,k,1,5)
                ti2 = cc(2,k,1,2)+cc(2,k,1,5)
                ti4 = cc(2,k,1,3)-cc(2,k,1,4)
                ti3 = cc(2,k,1,3)+cc(2,k,1,4)
                tr5 = cc(1,k,1,2)-cc(1,k,1,5)
                tr2 = cc(1,k,1,2)+cc(1,k,1,5)
                tr4 = cc(1,k,1,3)-cc(1,k,1,4)
                tr3 = cc(1,k,1,3)+cc(1,k,1,4)
                ch(1,k,1,1) = cc(1,k,1,1)+tr2+tr3
                ch(2,k,1,1) = cc(2,k,1,1)+ti2+ti3
                cr2 = cc(1,k,1,1)+TR11*tr2+TR12*tr3
                ci2 = cc(2,k,1,1)+TR11*ti2+TR12*ti3
                cr3 = cc(1,k,1,1)+TR12*tr2+TR11*tr3
                ci3 = cc(2,k,1,1)+TR12*ti2+TR11*ti3
                cr5 = TI11*tr5+TI12*tr4
                ci5 = TI11*ti5+TI12*ti4
                cr4 = TI12*tr5-TI11*tr4
                ci4 = TI12*ti5-TI11*ti4
                ch(1,k,2,1) = cr2-ci5
                ch(1,k,5,1) = cr2+ci5
                ch(2,k,2,1) = ci2+cr5
                ch(2,k,3,1) = ci3+cr4
                ch(1,k,3,1) = cr3-ci4
                ch(1,k,4,1) = cr3+ci4
                ch(2,k,4,1) = ci3-cr4
                ch(2,k,5,1) = ci2-cr5
            end do
            do i=2,ido
                do k=1,l1
                    ti5 = cc(2,k,i,2)-cc(2,k,i,5)
                    ti2 = cc(2,k,i,2)+cc(2,k,i,5)
                    ti4 = cc(2,k,i,3)-cc(2,k,i,4)
                    ti3 = cc(2,k,i,3)+cc(2,k,i,4)
                    tr5 = cc(1,k,i,2)-cc(1,k,i,5)
                    tr2 = cc(1,k,i,2)+cc(1,k,i,5)
                    tr4 = cc(1,k,i,3)-cc(1,k,i,4)
                    tr3 = cc(1,k,i,3)+cc(1,k,i,4)
                    ch(1,k,1,i) = cc(1,k,i,1)+tr2+tr3
                    ch(2,k,1,i) = cc(2,k,i,1)+ti2+ti3
                    cr2 = cc(1,k,i,1)+TR11*tr2+TR12*tr3
                    ci2 = cc(2,k,i,1)+TR11*ti2+TR12*ti3
                    cr3 = cc(1,k,i,1)+TR12*tr2+TR11*tr3
                    ci3 = cc(2,k,i,1)+TR12*ti2+TR11*ti3
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
                    ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
                    ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
                    ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
                    ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
                    ch(1,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
                    ch(2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
                    ch(1,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
                    ch(2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
                end do
            end do
        end if

    end subroutine c1f5kf


    subroutine c1fgkf(ido, iip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa)

        integer (ip) ido
        integer (ip) in1
        integer (ip) in2
        integer (ip) iip
        integer (ip) l1
        integer (ip) lid

        real (wp) cc(in1,l1,iip,ido)
        real (wp) cc1(in1,lid,iip)
        real (wp) ch(in2,l1,ido,iip)
        real (wp) ch1(in2,lid,iip)
        real (wp) chold1
        real (wp) chold2
        integer (ip) i
        integer (ip) idlj
        integer (ip) iipp2
        integer (ip) iipph
        integer (ip) j
        integer (ip) jc
        integer (ip) ki
        integer (ip) l
        integer (ip) lc
        integer (ip) na
        real (wp) sn
        real (wp) wa(ido,iip-1,2)
        real (wp) wai
        real (wp) war

        iipp2 = iip+2
        iipph = (iip+1)/2
        do ki=1,lid
            ch1(1,ki,1) = cc1(1,ki,1)
            ch1(2,ki,1) = cc1(2,ki,1)
        end do

        do j=2,iipph
            jc = iipp2-j
            do ki=1,lid
                ch1(1,ki,j) =  cc1(1,ki,j)+cc1(1,ki,jc)
                ch1(1,ki,jc) = cc1(1,ki,j)-cc1(1,ki,jc)
                ch1(2,ki,j) =  cc1(2,ki,j)+cc1(2,ki,jc)
                ch1(2,ki,jc) = cc1(2,ki,j)-cc1(2,ki,jc)
            end do
        end do

        do j=2,iipph
            do ki=1,lid
                cc1(1,ki,1) = cc1(1,ki,1)+ch1(1,ki,j)
                cc1(2,ki,1) = cc1(2,ki,1)+ch1(2,ki,j)
            end do
        end do

        do l=2,iipph
            lc = iipp2-l
            do ki=1,lid
                cc1(1,ki,l) = ch1(1,ki,1)+wa(1,l-1,1)*ch1(1,ki,2)
                cc1(1,ki,lc) = -wa(1,l-1,2)*ch1(1,ki,iip)
                cc1(2,ki,l) = ch1(2,ki,1)+wa(1,l-1,1)*ch1(2,ki,2)
                cc1(2,ki,lc) = -wa(1,l-1,2)*ch1(2,ki,iip)
            end do
            do j=3,iipph
                jc = iipp2-j
                idlj = mod((l-1)*(j-1),iip)
                war = wa(1,idlj,1)
                wai = -wa(1,idlj,2)
                do ki=1,lid
                    cc1(1,ki,l) = cc1(1,ki,l)+war*ch1(1,ki,j)
                    cc1(1,ki,lc) = cc1(1,ki,lc)+wai*ch1(1,ki,jc)
                    cc1(2,ki,l) = cc1(2,ki,l)+war*ch1(2,ki,j)
                    cc1(2,ki,lc) = cc1(2,ki,lc)+wai*ch1(2,ki,jc)
                end do
            end do
        end do

        if (1 >= ido)then
            sn = ONE/(iip * l1)
            if (na /= 1) then
                do ki=1,lid
                    cc1(1,ki,1) = sn*cc1(1,ki,1)
                    cc1(2,ki,1) = sn*cc1(2,ki,1)
                end do
                do j=2,iipph
                    jc = iipp2-j
                    do ki=1,lid
                        chold1 = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
                        chold2 = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
                        cc1(1,ki,j) = chold1
                        cc1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
                        cc1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
                        cc1(1,ki,jc) = chold2
                    end do
                end do
            else
                do ki=1,lid
                    ch1(1,ki,1) = sn*cc1(1,ki,1)
                    ch1(2,ki,1) = sn*cc1(2,ki,1)
                end do
                do j=2,iipph
                    jc = iipp2-j
                    do ki=1,lid
                        ch1(1,ki,j) = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
                        ch1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
                        ch1(1,ki,jc) = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
                        ch1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
                    end do
                end do
            end if
        else
            do ki=1,lid
                ch1(1,ki,1) = cc1(1,ki,1)
                ch1(2,ki,1) = cc1(2,ki,1)
            end do
            do j=2,iipph
                jc = iipp2-j
                do ki=1,lid
                    ch1(1,ki,j) = cc1(1,ki,j)-cc1(2,ki,jc)
                    ch1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
                    ch1(1,ki,jc) = cc1(1,ki,j)+cc1(2,ki,jc)
                    ch1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
                end do
            end do
            do i=1,ido
                cc(1,:,1,i) = ch(1,:,i,1)
                cc(2,:,1,i) = ch(2,:,i,1)
            end do
            do j=2,iip
                cc(1,:,j,1) = ch(1,:,1,j)
                cc(2,:,j,1) = ch(2,:,1,j)
            end do
            do j=2,iip
                do i=2,ido
                    cc(1,:,j,i) = wa(i,j-1,1)*ch(1,:,i,j)+wa(i,j-1,2)*ch(2,:,i,j)
                    cc(2,:,j,i) = wa(i,j-1,1)*ch(2,:,i,j)-wa(i,j-1,2)*ch(1,:,i,j)
                end do
            end do
        end if

    end subroutine c1fgkf



end submodule complex_forward_1d
