submodule (complex_transform_routines) multiple_complex_backward

contains

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
        ! Local variables
        !--------------------------------------------------
        integer (ip) :: iw1, iw2
        !--------------------------------------------------

        !
        !==> Check validity of input arguments
        !
        if (lenc < (lot-1)*jump + inc*(n-1) + 1) then
            ierror = 1
            call fft_error_handler('cfftmb ', 6)
        else if (lensav < get_complex_nd_saved_workspace_length(n)) then
            ierror = 2
            call fft_error_handler('cfftmb ', 8)
        else if (lenwrk < get_complex_nd_workspace_length(n, lot)) then
            ierror = 3
            call fft_error_handler('cfftmb ', 10)
        else if (.not. fft_consistent(inc,jump,n,lot)) then
            ierror = 4
            call fft_error_handler('cfftmb ', -1)
        else
            ierror = 0
        end if

        !
        !==> Perform transform
        !
        if (n == 1) return

        ! Set workspace index pointer
        iw1 = 2*n+1
        iw2 = iw1 + 1

        call cmfm1b(lot, jump, n, inc, c, work, wsave, wsave(iw1), wsave(iw2))

    end subroutine cfftmb



    subroutine cmfm1b(lot, jump, n, inc, c, ch, wa, fnf, fac)
        !--------------------------------------------------
        ! Dummy arguments
        !--------------------------------------------------
        real (wp) c(2,*)
        real (wp) ch(*)
        real (wp) fac(*)
        real (wp) fnf
        integer (ip) ido
        integer (ip) inc
        integer (ip) iip
        integer (ip) iw
        integer (ip) jump
        integer (ip) k1
        integer (ip) l1
        integer (ip) l2
        integer (ip) lid
        integer (ip) lot
        integer (ip) n
        integer (ip) na
        integer (ip) nbr
        integer (ip) nf
        real (wp) wa(*)

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
                    call cmf2kb(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (2)
                    call cmf2kb(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (3)
                    call cmf3kb(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (4)
                    call cmf3kb(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (5)
                    call cmf4kb(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (6)
                    call cmf4kb(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (7)
                    call cmf5kb(lot,ido,l1,na,c,jump,inc,ch,1,lot,wa(iw))
                case (8)
                    call cmf5kb(lot,ido,l1,na,ch,1,lot,c,jump,inc,wa(iw))
                case (9)
                    call cmfgkb(lot,ido,iip,l1,lid,na,c,c,jump,inc,ch,ch,1,lot,wa(iw))
                case (10)
                    call cmfgkb(lot,ido,iip,l1,lid,na,ch,ch,1,lot,c,c, &
                        jump,inc,wa(iw))
            end select

            l1 = l2
            iw = iw+(iip-1)*(2*ido)

            if (iip <= 5) na = 1-na

        end do


    end subroutine cmfm1b


    subroutine cmf2kb(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

        integer (ip) ido
        integer (ip) in1
        integer (ip) in2
        integer (ip) l1

        real (wp) cc(2,in1,l1,ido,2)
        real (wp) ch(2,in2,l1,2,ido)
        real (wp) chold1
        real (wp) chold2
        integer (ip) i
        integer (ip) im1
        integer (ip) im2
        integer (ip) k
        integer (ip) lot
        integer (ip) m1
        integer (ip) m1d
        integer (ip) m2
        integer (ip) m2s
        integer (ip) na
        real (wp) ti2
        real (wp) tr2
        real (wp) wa(ido,1,2)

        m1d = (lot-1)*im1+1
        m2s = 1-im2

        if (1 >= ido .and. na /= 1) then
            do k=1,l1
                do m1=1,m1d,im1
                    chold1 = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
                    cc(1,m1,k,1,2) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
                    cc(1,m1,k,1,1) = chold1
                    chold2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
                    cc(2,m1,k,1,2) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
                    cc(2,m1,k,1,1) = chold2
                end do
            end do
        else
            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+cc(1,m1,k,1,2)
                    ch(1,m2,k,2,1) = cc(1,m1,k,1,1)-cc(1,m1,k,1,2)
                    ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+cc(2,m1,k,1,2)
                    ch(2,m2,k,2,1) = cc(2,m1,k,1,1)-cc(2,m1,k,1,2)
                end do
            end do

            do i=2,ido
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+cc(1,m1,k,i,2)
                        tr2 = cc(1,m1,k,i,1)-cc(1,m1,k,i,2)
                        ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+cc(2,m1,k,i,2)
                        ti2 = cc(2,m1,k,i,1)-cc(2,m1,k,i,2)
                        ch(2,m2,k,2,i) = wa(i,1,1)*ti2+wa(i,1,2)*tr2
                        ch(1,m2,k,2,i) = wa(i,1,1)*tr2-wa(i,1,2)*ti2
                    end do
                end do
            end do
        end if

    end subroutine cmf2kb

    subroutine cmf3kb(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

        integer (ip) ido
        integer (ip) in1
        integer (ip) in2
        integer (ip) l1

        real (wp) cc(2,in1,l1,ido,3)
        real (wp) ch(2,in2,l1,3,ido)
        real (wp) ci2
        real (wp) ci3
        real (wp) cr2
        real (wp) cr3
        real (wp) di2
        real (wp) di3
        real (wp) dr2
        real (wp) dr3
        integer (ip) i
        integer (ip) im1
        integer (ip) im2
        integer (ip) k
        integer (ip) lot
        integer (ip) m1
        integer (ip) m1d
        integer (ip) m2
        integer (ip) m2s
        integer (ip) na
        real (wp), parameter :: SQRT3 = sqrt(THREE)
        real (wp), parameter :: TAUI =  SQRT3/2 !0.866025403784439_wp
        real (wp), parameter :: TAUR = -HALF
        real (wp) ti2
        real (wp) tr2
        real (wp) wa(ido,2,2)

        m1d = (lot-1)*im1+1
        m2s = 1-im2

        if (1 >= ido .and. na /= 1) then
            do k=1,l1
                do m1=1,m1d,im1
                    tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
                    cr2 = cc(1,m1,k,1,1)+TAUR*tr2
                    cc(1,m1,k,1,1) = cc(1,m1,k,1,1)+tr2
                    ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
                    ci2 = cc(2,m1,k,1,1)+TAUR*ti2
                    cc(2,m1,k,1,1) = cc(2,m1,k,1,1)+ti2
                    cr3 = TAUI*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
                    ci3 = TAUI*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
                    cc(1,m1,k,1,2) = cr2-ci3
                    cc(1,m1,k,1,3) = cr2+ci3
                    cc(2,m1,k,1,2) = ci2+cr3
                    cc(2,m1,k,1,3) = ci2-cr3
                end do
            end do
        else
            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,3)
                    cr2 = cc(1,m1,k,1,1)+TAUR*tr2
                    ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2
                    ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,3)
                    ci2 = cc(2,m1,k,1,1)+TAUR*ti2
                    ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2
                    cr3 = TAUI*(cc(1,m1,k,1,2)-cc(1,m1,k,1,3))
                    ci3 = TAUI*(cc(2,m1,k,1,2)-cc(2,m1,k,1,3))
                    ch(1,m2,k,2,1) = cr2-ci3
                    ch(1,m2,k,3,1) = cr2+ci3
                    ch(2,m2,k,2,1) = ci2+cr3
                    ch(2,m2,k,3,1) = ci2-cr3
                end do
            end do

            do i=2,ido
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,3)
                        cr2 = cc(1,m1,k,i,1)+TAUR*tr2
                        ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2
                        ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,3)
                        ci2 = cc(2,m1,k,i,1)+TAUR*ti2
                        ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2
                        cr3 = TAUI*(cc(1,m1,k,i,2)-cc(1,m1,k,i,3))
                        ci3 = TAUI*(cc(2,m1,k,i,2)-cc(2,m1,k,i,3))
                        dr2 = cr2-ci3
                        dr3 = cr2+ci3
                        di2 = ci2+cr3
                        di3 = ci2-cr3
                        ch(2,m2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
                        ch(1,m2,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
                        ch(2,m2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
                        ch(1,m2,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
                    end do
                end do
            end do
        end if
    end subroutine cmf3kb

    subroutine cmf4kb(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

        integer (ip) ido
        integer (ip) in1
        integer (ip) in2
        integer (ip) l1

        real (wp) cc(2,in1,l1,ido,4)
        real (wp) ch(2,in2,l1,4,ido)
        real (wp) ci2
        real (wp) ci3
        real (wp) ci4
        real (wp) cr2
        real (wp) cr3
        real (wp) cr4
        integer (ip) i
        integer (ip) im1
        integer (ip) im2
        integer (ip) k
        integer (ip) lot
        integer (ip) m1
        integer (ip) m1d
        integer (ip) m2
        integer (ip) m2s
        integer (ip) na
        real (wp) ti1
        real (wp) ti2
        real (wp) ti3
        real (wp) ti4
        real (wp) tr1
        real (wp) tr2
        real (wp) tr3
        real (wp) tr4
        real (wp) wa(ido,3,2)

        m1d = (lot-1)*im1+1
        m2s = 1-im2

        if (1 >= ido .and. na /= 1) then
            do k=1,l1
                do m1=1,m1d,im1
                    ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
                    ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
                    tr4 = cc(2,m1,k,1,4)-cc(2,m1,k,1,2)
                    ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
                    tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
                    tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
                    ti4 = cc(1,m1,k,1,2)-cc(1,m1,k,1,4)
                    tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
                    cc(1,m1,k,1,1) = tr2+tr3
                    cc(1,m1,k,1,3) = tr2-tr3
                    cc(2,m1,k,1,1) = ti2+ti3
                    cc(2,m1,k,1,3) = ti2-ti3
                    cc(1,m1,k,1,2) = tr1+tr4
                    cc(1,m1,k,1,4) = tr1-tr4
                    cc(2,m1,k,1,2) = ti1+ti4
                    cc(2,m1,k,1,4) = ti1-ti4
                end do
            end do
        else
            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ti1 = cc(2,m1,k,1,1)-cc(2,m1,k,1,3)
                    ti2 = cc(2,m1,k,1,1)+cc(2,m1,k,1,3)
                    tr4 = cc(2,m1,k,1,4)-cc(2,m1,k,1,2)
                    ti3 = cc(2,m1,k,1,2)+cc(2,m1,k,1,4)
                    tr1 = cc(1,m1,k,1,1)-cc(1,m1,k,1,3)
                    tr2 = cc(1,m1,k,1,1)+cc(1,m1,k,1,3)
                    ti4 = cc(1,m1,k,1,2)-cc(1,m1,k,1,4)
                    tr3 = cc(1,m1,k,1,2)+cc(1,m1,k,1,4)
                    ch(1,m2,k,1,1) = tr2+tr3
                    ch(1,m2,k,3,1) = tr2-tr3
                    ch(2,m2,k,1,1) = ti2+ti3
                    ch(2,m2,k,3,1) = ti2-ti3
                    ch(1,m2,k,2,1) = tr1+tr4
                    ch(1,m2,k,4,1) = tr1-tr4
                    ch(2,m2,k,2,1) = ti1+ti4
                    ch(2,m2,k,4,1) = ti1-ti4
                end do
            end do

            do i=2,ido
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ti1 = cc(2,m1,k,i,1)-cc(2,m1,k,i,3)
                        ti2 = cc(2,m1,k,i,1)+cc(2,m1,k,i,3)
                        ti3 = cc(2,m1,k,i,2)+cc(2,m1,k,i,4)
                        tr4 = cc(2,m1,k,i,4)-cc(2,m1,k,i,2)
                        tr1 = cc(1,m1,k,i,1)-cc(1,m1,k,i,3)
                        tr2 = cc(1,m1,k,i,1)+cc(1,m1,k,i,3)
                        ti4 = cc(1,m1,k,i,2)-cc(1,m1,k,i,4)
                        tr3 = cc(1,m1,k,i,2)+cc(1,m1,k,i,4)
                        ch(1,m2,k,1,i) = tr2+tr3
                        cr3 = tr2-tr3
                        ch(2,m2,k,1,i) = ti2+ti3
                        ci3 = ti2-ti3
                        cr2 = tr1+tr4
                        cr4 = tr1-tr4
                        ci2 = ti1+ti4
                        ci4 = ti1-ti4
                        ch(1,m2,k,2,i) = wa(i,1,1)*cr2-wa(i,1,2)*ci2
                        ch(2,m2,k,2,i) = wa(i,1,1)*ci2+wa(i,1,2)*cr2
                        ch(1,m2,k,3,i) = wa(i,2,1)*cr3-wa(i,2,2)*ci3
                        ch(2,m2,k,3,i) = wa(i,2,1)*ci3+wa(i,2,2)*cr3
                        ch(1,m2,k,4,i) = wa(i,3,1)*cr4-wa(i,3,2)*ci4
                        ch(2,m2,k,4,i) = wa(i,3,1)*ci4+wa(i,3,2)*cr4
                    end do
                end do
            end do
        end if

    end subroutine cmf4kb

    subroutine cmf5kb(lot, ido, l1, na, cc, im1, in1, ch, im2, in2, wa)

        integer (ip) ido
        integer (ip) in1
        integer (ip) in2
        integer (ip) l1

        real (wp) cc(2,in1,l1,ido,5)
        real (wp) ch(2,in2,l1,5,ido)
        real (wp) chold1
        real (wp) chold2
        real (wp) ci2
        real (wp) ci3
        real (wp) ci4
        real (wp) ci5
        real (wp) cr2
        real (wp) cr3
        real (wp) cr4
        real (wp) cr5
        real (wp) di2
        real (wp) di3
        real (wp) di4
        real (wp) di5
        real (wp) dr2
        real (wp) dr3
        real (wp) dr4
        real (wp) dr5
        integer (ip) i
        integer (ip) im1
        integer (ip) im2
        integer (ip) k
        integer (ip) lot
        integer (ip) m1
        integer (ip) m1d
        integer (ip) m2
        integer (ip) m2s
        integer (ip) na
        real (wp) ti2
        real (wp) ti3
        real (wp) ti4
        real (wp) ti5
        real (wp) tr2
        real (wp) tr3
        real (wp) tr4
        real (wp) tr5
        real (wp), parameter :: SQRT5 = sqrt(FIVE)
        real (wp), parameter :: SQRT5_PLUS_5 = SQRT5 + FIVE
        real (wp), parameter :: TI11 = sqrt(SQRT5_PLUS_5/2)/2             ! 0.9510565162951536_wp
        real (wp), parameter :: TI12 = sqrt(FIVE/(TWO*SQRT5_PLUS_5)) ! 0.5877852522924731_wp
        real (wp), parameter :: TR11 =  (SQRT5 - ONE)/4                 ! 0.3090169943749474_wp
        real (wp), parameter :: TR12 = -(ONE + SQRT5)/4                 !-0.8090169943749474_wp

        real (wp) wa(ido,4,2)

        m1d = (lot-1)*im1+1
        m2s = 1-im2

        if (1 >= ido .and. na /= 1) then
            do k=1,l1
                do m1=1,m1d,im1
                    ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
                    ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
                    ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
                    ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
                    tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
                    tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
                    tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
                    tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
                    chold1 = cc(1,m1,k,1,1)+tr2+tr3
                    chold2 = cc(2,m1,k,1,1)+ti2+ti3
                    cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
                    ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
                    cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
                    ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
                    cc(1,m1,k,1,1) = chold1
                    cc(2,m1,k,1,1) = chold2
                    cr5 = ti11*tr5+ti12*tr4
                    ci5 = ti11*ti5+ti12*ti4
                    cr4 = ti12*tr5-ti11*tr4
                    ci4 = ti12*ti5-ti11*ti4
                    cc(1,m1,k,1,2) = cr2-ci5
                    cc(1,m1,k,1,5) = cr2+ci5
                    cc(2,m1,k,1,2) = ci2+cr5
                    cc(2,m1,k,1,3) = ci3+cr4
                    cc(1,m1,k,1,3) = cr3-ci4
                    cc(1,m1,k,1,4) = cr3+ci4
                    cc(2,m1,k,1,4) = ci3-cr4
                    cc(2,m1,k,1,5) = ci2-cr5
                end do
            end do
        else
            do k=1,l1
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ti5 = cc(2,m1,k,1,2)-cc(2,m1,k,1,5)
                    ti2 = cc(2,m1,k,1,2)+cc(2,m1,k,1,5)
                    ti4 = cc(2,m1,k,1,3)-cc(2,m1,k,1,4)
                    ti3 = cc(2,m1,k,1,3)+cc(2,m1,k,1,4)
                    tr5 = cc(1,m1,k,1,2)-cc(1,m1,k,1,5)
                    tr2 = cc(1,m1,k,1,2)+cc(1,m1,k,1,5)
                    tr4 = cc(1,m1,k,1,3)-cc(1,m1,k,1,4)
                    tr3 = cc(1,m1,k,1,3)+cc(1,m1,k,1,4)
                    ch(1,m2,k,1,1) = cc(1,m1,k,1,1)+tr2+tr3
                    ch(2,m2,k,1,1) = cc(2,m1,k,1,1)+ti2+ti3
                    cr2 = cc(1,m1,k,1,1)+tr11*tr2+tr12*tr3
                    ci2 = cc(2,m1,k,1,1)+tr11*ti2+tr12*ti3
                    cr3 = cc(1,m1,k,1,1)+tr12*tr2+tr11*tr3
                    ci3 = cc(2,m1,k,1,1)+tr12*ti2+tr11*ti3
                    cr5 = ti11*tr5+ti12*tr4
                    ci5 = ti11*ti5+ti12*ti4
                    cr4 = ti12*tr5-ti11*tr4
                    ci4 = ti12*ti5-ti11*ti4
                    ch(1,m2,k,2,1) = cr2-ci5
                    ch(1,m2,k,5,1) = cr2+ci5
                    ch(2,m2,k,2,1) = ci2+cr5
                    ch(2,m2,k,3,1) = ci3+cr4
                    ch(1,m2,k,3,1) = cr3-ci4
                    ch(1,m2,k,4,1) = cr3+ci4
                    ch(2,m2,k,4,1) = ci3-cr4
                    ch(2,m2,k,5,1) = ci2-cr5
                end do
            end do

            do i=2,ido
                do k=1,l1
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ti5 = cc(2,m1,k,i,2)-cc(2,m1,k,i,5)
                        ti2 = cc(2,m1,k,i,2)+cc(2,m1,k,i,5)
                        ti4 = cc(2,m1,k,i,3)-cc(2,m1,k,i,4)
                        ti3 = cc(2,m1,k,i,3)+cc(2,m1,k,i,4)
                        tr5 = cc(1,m1,k,i,2)-cc(1,m1,k,i,5)
                        tr2 = cc(1,m1,k,i,2)+cc(1,m1,k,i,5)
                        tr4 = cc(1,m1,k,i,3)-cc(1,m1,k,i,4)
                        tr3 = cc(1,m1,k,i,3)+cc(1,m1,k,i,4)
                        ch(1,m2,k,1,i) = cc(1,m1,k,i,1)+tr2+tr3
                        ch(2,m2,k,1,i) = cc(2,m1,k,i,1)+ti2+ti3
                        cr2 = cc(1,m1,k,i,1)+tr11*tr2+tr12*tr3
                        ci2 = cc(2,m1,k,i,1)+tr11*ti2+tr12*ti3
                        cr3 = cc(1,m1,k,i,1)+tr12*tr2+tr11*tr3
                        ci3 = cc(2,m1,k,i,1)+tr12*ti2+tr11*ti3
                        cr5 = ti11*tr5+ti12*tr4
                        ci5 = ti11*ti5+ti12*ti4
                        cr4 = ti12*tr5-ti11*tr4
                        ci4 = ti12*ti5-ti11*ti4
                        dr3 = cr3-ci4
                        dr4 = cr3+ci4
                        di3 = ci3+cr4
                        di4 = ci3-cr4
                        dr5 = cr2+ci5
                        dr2 = cr2-ci5
                        di5 = ci2-cr5
                        di2 = ci2+cr5
                        ch(1,m2,k,2,i) = wa(i,1,1)*dr2-wa(i,1,2)*di2
                        ch(2,m2,k,2,i) = wa(i,1,1)*di2+wa(i,1,2)*dr2
                        ch(1,m2,k,3,i) = wa(i,2,1)*dr3-wa(i,2,2)*di3
                        ch(2,m2,k,3,i) = wa(i,2,1)*di3+wa(i,2,2)*dr3
                        ch(1,m2,k,4,i) = wa(i,3,1)*dr4-wa(i,3,2)*di4
                        ch(2,m2,k,4,i) = wa(i,3,1)*di4+wa(i,3,2)*dr4
                        ch(1,m2,k,5,i) = wa(i,4,1)*dr5-wa(i,4,2)*di5
                        ch(2,m2,k,5,i) = wa(i,4,1)*di5+wa(i,4,2)*dr5
                    end do
                end do
            end do
        end if

    end subroutine cmf5kb

    subroutine cmfgkb(lot, ido, iip, l1, lid, na, cc, cc1, im1, in1, &
        ch, ch1, im2, in2, wa)

        integer (ip) ido
        integer (ip) in1
        integer (ip) in2
        integer (ip) iip
        integer (ip) l1
        integer (ip) lid

        real (wp) cc(2,in1,l1,iip,ido)
        real (wp) cc1(2,in1,lid,iip)
        real (wp) ch(2,in2,l1,ido,iip)
        real (wp) ch1(2,in2,lid,iip)
        real (wp) chold1
        real (wp) chold2
        integer (ip) i
        integer (ip) idlj
        integer (ip) im1
        integer (ip) im2
        integer (ip) iipp2
        integer (ip) iipph
        integer (ip) j
        integer (ip) jc
        integer (ip) k
        integer (ip) ki
        integer (ip) l
        integer (ip) lc
        integer (ip) lot
        integer (ip) m1
        integer (ip) m1d
        integer (ip) m2
        integer (ip) m2s
        integer (ip) na
        real (wp) wa(ido,iip-1,2)
        real (wp) wai
        real (wp) war

        m1d = (lot-1)*im1+1
        m2s = 1-im2
        iipp2 = iip+2
        iipph = (iip+1)/2

        do ki=1,lid
            m2 = m2s
            do m1=1,m1d,im1
                m2 = m2+im2
                ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
                ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
            end do
        end do

        do j=2,iipph
            jc = iipp2-j
            do ki=1,lid
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch1(1,m2,ki,j) =  cc1(1,m1,ki,j)+cc1(1,m1,ki,jc)
                    ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)-cc1(1,m1,ki,jc)
                    ch1(2,m2,ki,j) =  cc1(2,m1,ki,j)+cc1(2,m1,ki,jc)
                    ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(2,m1,ki,jc)
                end do
            end do
        end do

        do j=2,iipph
            do ki=1,lid
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    cc1(1,m1,ki,1) = cc1(1,m1,ki,1)+ch1(1,m2,ki,j)
                    cc1(2,m1,ki,1) = cc1(2,m1,ki,1)+ch1(2,m2,ki,j)
                end do
            end do
        end do

        do l=2,iipph
            lc = iipp2-l
            do ki=1,lid
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    cc1(1,m1,ki,l) = ch1(1,m2,ki,1)+wa(1,l-1,1)*ch1(1,m2,ki,2)
                    cc1(1,m1,ki,lc) = wa(1,l-1,2)*ch1(1,m2,ki,iip)
                    cc1(2,m1,ki,l) = ch1(2,m2,ki,1)+wa(1,l-1,1)*ch1(2,m2,ki,2)
                    cc1(2,m1,ki,lc) = wa(1,l-1,2)*ch1(2,m2,ki,iip)
                end do
            end do
            do j=3,iipph
                jc = iipp2-j
                idlj = mod((l-1)*(j-1),iip)
                war = wa(1,idlj,1)
                wai = wa(1,idlj,2)
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        cc1(1,m1,ki,l) = cc1(1,m1,ki,l)+war*ch1(1,m2,ki,j)
                        cc1(1,m1,ki,lc) = cc1(1,m1,ki,lc)+wai*ch1(1,m2,ki,jc)
                        cc1(2,m1,ki,l) = cc1(2,m1,ki,l)+war*ch1(2,m2,ki,j)
                        cc1(2,m1,ki,lc) = cc1(2,m1,ki,lc)+wai*ch1(2,m2,ki,jc)
                    end do
                end do
            end do
        end do

        if (1 >= ido .and. na /= 1) then
            do j=2,iipph
                jc = iipp2-j
                do ki=1,lid
                    do m1=1,m1d,im1
                        chold1 = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
                        chold2 = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
                        cc1(1,m1,ki,j) = chold1
                        cc1(2,m1,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
                        cc1(2,m1,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
                        cc1(1,m1,ki,jc) = chold2
                    end do
                end do
            end do
        else
            do ki=1,lid
                m2 = m2s
                do m1=1,m1d,im1
                    m2 = m2+im2
                    ch1(1,m2,ki,1) = cc1(1,m1,ki,1)
                    ch1(2,m2,ki,1) = cc1(2,m1,ki,1)
                end do
            end do

            do j=2,iipph
                jc = iipp2-j
                do ki=1,lid
                    m2 = m2s
                    do m1=1,m1d,im1
                        m2 = m2+im2
                        ch1(1,m2,ki,j) = cc1(1,m1,ki,j)-cc1(2,m1,ki,jc)
                        ch1(1,m2,ki,jc) = cc1(1,m1,ki,j)+cc1(2,m1,ki,jc)
                        ch1(2,m2,ki,jc) = cc1(2,m1,ki,j)-cc1(1,m1,ki,jc)
                        ch1(2,m2,ki,j) = cc1(2,m1,ki,j)+cc1(1,m1,ki,jc)
                    end do
                end do
            end do

            if (ido /= 1) then
                do i=1,ido
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc(1,m1,k,1,i) = ch(1,m2,k,i,1)
                            cc(2,m1,k,1,i) = ch(2,m2,k,i,1)
                        end do
                    end do
                end do

                do j=2,iip
                    do k=1,l1
                        m2 = m2s
                        do m1=1,m1d,im1
                            m2 = m2+im2
                            cc(1,m1,k,j,1) = ch(1,m2,k,1,j)
                            cc(2,m1,k,j,1) = ch(2,m2,k,1,j)
                        end do
                    end do
                end do

                do j=2,iip
                    do i=2,ido
                        do k=1,l1
                            m2 = m2s
                            do m1=1,m1d,im1
                                m2 = m2+im2
                                cc(1,m1,k,j,i) = wa(i,j-1,1)*ch(1,m2,k,i,j) &
                                    -wa(i,j-1,2)*ch(2,m2,k,i,j)
                                cc(2,m1,k,j,i) = wa(i,j-1,1)*ch(2,m2,k,i,j) &
                                    +wa(i,j-1,2)*ch(1,m2,k,i,j)
                            end do
                        end do
                    end do
                end do
            end if
        end if

    end subroutine cmfgkb






end submodule multiple_complex_backward
