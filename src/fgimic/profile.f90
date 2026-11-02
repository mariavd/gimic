!
! Current profiles: the integrated current resolved into slices along the
! width (axis 2) and/or the height (axis 1) of the integration plane.
!
! The slices are the blocks of the composite Gauss quadrature that the
! integration grid is built from (Grid.slices sets one block per slice), so
! a profile is obtained from the column sums of a single integration over
! the whole plane; no extra work is done here beyond summing and writing.
!
! Files written (a spin suffix _alpha, _beta or _spindens is added for
! open-shell cases):
!
!   current_profile.dat   one line per slice along the width, from the start
!                         of the integration to its end: distance of the
!                         slice start from the start of the integration, net,
!                         positive (diatropic) and negative (paratropic)
!                         current through the slice in nA/T. No header, the
!                         format the current-profile scripts have always
!                         produced.
!   vertical_profile.dat  the same for the slices along the height, when the
!                         height is sliced: distance from the grid origin along
!                         the height axis, net, positive and negative current.
!   profile.dat           the same with the geometry (slice centres in space,
!                         basis vectors, field, reporting direction), the
!                         slices along the height if any, and the current
!                         per unit length at every Gauss node along the
!                         width, so that any sub-range can be re-integrated.
!   profile_points.dat    the analysis of the profile along the width: maxima
!                         and minima of the net, diatropic and paratropic
!                         current per bohr, the zero crossings of the net
!                         current, and the current integrated between them.
!   squares_profile.dat   when both axes are sliced: the net current through
!                         each rectangle, as a matrix (rows: height slices,
!                         columns: width slices) for plotting with image
!                         tools, followed by a list with coordinates.
!
module profile_module
    use globals_module
    use settings_module
    use grid_class
    use teletype_module
    implicit none

    public write_profile
    private

    integer, parameter :: PFD=71
contains

    ! blk(:,b,j): (net, positive, negative) current through column j of the
    ! grid, resolved into the blocks b of the height axis and integrated along
    ! the height, but not yet weighted along the width (see integral.f90)
    subroutine write_profile(g, blk, normal, bb, spin)
        type(grid_t), intent(in) :: g
        real(DP), dimension(:,:,:), intent(in) :: blk
        real(DP), dimension(3), intent(in) :: normal, bb
        character(*), intent(in) :: spin

        integer(I4) :: nb1, nb2, p2, j, b1, b2, order, n
        real(DP), dimension(:,:,:), allocatable :: s
        real(DP), dimension(:,:), allocatable :: node
        real(DP), dimension(3) :: tot, c
        real(DP) :: w, d1, d2, x, t, dir
        character(16) :: suffix
        logical :: from_in

        nb1=size(blk,2)
        p2=size(blk,3)
        order=g%gauss_order
        nb2=max(1, g%nslice(2))
        d1=g%l(1)/real(nb1,DP)
        d2=g%l(2)/real(nb2,DP)

        suffix=''
        if (trim(spin) /= 'total') suffix='_'//trim(spin)

        ! The grid origin is at t=0 along the width. For a bond grid that is
        ! the "out" edge (origin = centre - out*basv2), and the profile is
        ! written from the "in" edge towards it, like the slice scripts did.
        ! For other grids it runs from the origin along basv2.
        from_in = (trim(g%mode) == 'bond')

        ! slice integrals (au) and the current per unit length at the nodes
        allocate(s(3,nb1,nb2), node(3,p2))
        s=0.d0
        do j=1,p2
            w=get_weight(g,j,2)
            b2=1
            if (nb2 > 1) b2=(j-1)/order+1
            s(:,:,b2)=s(:,:,b2)+w*blk(:,:,j)
            node(:,j)=sum(blk(:,:,j),dim=2)
        end do
        tot=sum(sum(s,dim=3),dim=2)

        ! ---- current_profile.dat: the format of the slice scripts --------
        if (g%nslice(2) > 0) then
            open(PFD, file='current_profile'//trim(suffix)//'.dat')
            do n=1,nb2
                b2=slice_index(n, nb2, from_in)
                x=real(n-1,DP)*d2
                write(PFD, '(f5.2,a,f11.8,a,f11.8,a,f11.8)') x, char(9), &
                    au2si(sum(s(1,:,b2))), char(9), au2si(sum(s(2,:,b2))), &
                    char(9), au2si(sum(s(3,:,b2)))
            end do
            close(PFD)
            call msg_out('Current profile written to current_profile'//trim(suffix)//'.dat')
        end if

        if (g%nslice(2) > 2) call analyse_profile(g, s, d2, from_in, suffix)

        ! ---- vertical_profile.dat: the slices along the height -----------
        if (g%nslice(1) > 0) then
            open(PFD, file='vertical_profile'//trim(suffix)//'.dat')
            do b1=1,nb1
                x=real(b1-1,DP)*d1
                write(PFD, '(f5.2,a,f11.8,a,f11.8,a,f11.8)') x, char(9), &
                    au2si(sum(s(1,b1,:))), char(9), au2si(sum(s(2,b1,:))), &
                    char(9), au2si(sum(s(3,b1,:)))
            end do
            close(PFD)
            call msg_out('Vertical profile written to vertical_profile'//trim(suffix)//'.dat')
        end if

        ! ---- profile.dat: geometry, slices along both axes, nodes --------
        open(PFD, file='profile'//trim(suffix)//'.dat')
        write(PFD, '(a)') '# GIMIC current profile; lengths in bohr, currents in nA/T'
        write(PFD, '(a,3f14.8)') '# grid origin        :', g%origin
        write(PFD, '(a,3f14.8)') '# basv1 (height)     :', g%basv(:,1)
        write(PFD, '(a,3f14.8)') '# basv2 (width)      :', g%basv(:,2)
        write(PFD, '(a,3f14.8)') '# basv3 (normal)     :', g%basv(:,3)
        write(PFD, '(a,3f14.8)') '# lengths            :', g%l
        write(PFD, '(a,3f14.8)') '# magnetic field     :', bb
        write(PFD, '(a,3f14.8)') '# current counted along:', normal
        if (g%has_towards) then
            write(PFD, '(a,3f14.8)') '# reference centre   :', g%towards
            write(PFD, '(a)') '# sign: positive for flow in the diatropic sense around the reference centre'
        else
            write(PFD, '(a)') '# sign: flux along basv3 (no reference centre given)'
        end if
        if (from_in) then
            write(PFD, '(a)') '# x: distance from the start of the integration (the "in" edge, t = l2) towards the "out" edge'
        else
            write(PFD, '(a)') '# x: distance from the grid origin along basv2'
        end if
        ! end points of the profile on the mid-height line of the plane
        if (from_in) then
            c=slice_centre(g, 0.5d0*g%l(1), g%l(2))
            write(PFD, '(a,3f14.8)') '# profile start      :', c
            c=slice_centre(g, 0.5d0*g%l(1), 0.d0)
            write(PFD, '(a,3f14.8)') '# profile end        :', c
        else
            c=slice_centre(g, 0.5d0*g%l(1), 0.d0)
            write(PFD, '(a,3f14.8)') '# profile start      :', c
            c=slice_centre(g, 0.5d0*g%l(1), g%l(2))
            write(PFD, '(a,3f14.8)') '# profile end        :', c
        end if
        write(PFD, '(a,i6,a,f10.5,a,i6,a,f10.5,a)') '# slices: ', nb2, &
            ' along the width of', d2, ' bohr;', g%nslice(1), &
            ' along the height of', d1, ' bohr'
        write(PFD, '(a,3f14.8)') '# total current      :', au2si(tot(1)), &
            au2si(tot(2)), au2si(tot(3))
        write(PFD, '(a)') '#'

        if (g%nslice(2) > 0) then
            write(PFD, '(a)') '# Slices along the width'
            write(PFD, '(a)') '#    n     x_start    x_centre' // &
                '        cx            cy            cz' // &
                '           net           dia          para' // &
                '      dJ/dx [nA/T/bohr]'
            do n=1,nb2
                b2=slice_index(n, nb2, from_in)
                x=real(n-1,DP)*d2
                c=slice_centre(g, 0.5d0*g%l(1), (real(b2,DP)-0.5d0)*d2)
                write(PFD, '(i6,2f12.5,3f14.6,3f14.8,f16.8)') n, x, x+0.5d0*d2, c, &
                    au2si(sum(s(1,:,b2))), au2si(sum(s(2,:,b2))), &
                    au2si(sum(s(3,:,b2))), au2si(sum(s(1,:,b2)))/d2
            end do
            write(PFD, '(a)') ''
        end if

        if (g%nslice(1) > 0) then
            write(PFD, '(a)') '# Slices along the height (from the grid origin along basv1)'
            write(PFD, '(a)') '#    n     s_start    s_centre' // &
                '        cx            cy            cz' // &
                '           net           dia          para' // &
                '      dJ/ds [nA/T/bohr]'
            do b1=1,nb1
                x=real(b1-1,DP)*d1
                c=slice_centre(g, (real(b1,DP)-0.5d0)*d1, 0.5d0*g%l(2))
                write(PFD, '(i6,2f12.5,3f14.6,3f14.8,f16.8)') b1, x, x+0.5d0*d1, c, &
                    au2si(sum(s(1,b1,:))), au2si(sum(s(2,b1,:))), &
                    au2si(sum(s(3,b1,:))), au2si(sum(s(1,b1,:)))/d1
            end do
            write(PFD, '(a)') ''
        end if

        write(PFD, '(a)') '# Current per unit length at the Gauss nodes along the width [nA/T/bohr]'
        write(PFD, '(a)') '#    j           x        weight' // &
                '           net           dia          para'
        if (from_in) then
            dir=-1.d0
        else
            dir=1.d0
        end if
        do n=1,p2
            if (from_in) then
                j=p2+1-n
                t=g%l(2)-g%gdata(2)%pts(j)
            else
                j=n
                t=g%gdata(2)%pts(j)
            end if
            write(PFD, '(i6,f12.6,f14.8,3f14.8)') j, t, get_weight(g,j,2), &
                au2si(node(1,j)), au2si(node(2,j)), au2si(node(3,j))
        end do
        close(PFD)
        call msg_out('Profile details written to profile'//trim(suffix)//'.dat')

        ! ---- squares_profile.dat ---------------------------------------
        if (g%nslice(1) > 0 .and. g%nslice(2) > 0) then
            open(PFD, file='squares_profile'//trim(suffix)//'.dat')
            write(PFD, '(a)') '# GIMIC squares profile: net current [nA/T] through each rectangle'
            write(PFD, '(a,i6,a,i6,a)') '# matrix: ', nb1, ' rows (height slices, from the grid origin along basv1) x ', &
                nb2, ' columns (width slices, in the order of current_profile.dat)'
            do b1=1,nb1
                write(PFD, '(*(f14.8))') (au2si(sum(s(1,b1,slice_index(n,nb2,from_in):slice_index(n,nb2,from_in)))), n=1,nb2)
            end do
            write(PFD, '(a)') ''
            write(PFD, '(a)') '# list: x_centre  s_centre  cx  cy  cz  net  dia  para'
            do b1=1,nb1
                do n=1,nb2
                    b2=slice_index(n, nb2, from_in)
                    c=slice_centre(g, (real(b1,DP)-0.5d0)*d1, (real(b2,DP)-0.5d0)*d2)
                    write(PFD, '(2f12.5,3f14.6,3f14.8)') (real(n,DP)-0.5d0)*d2, &
                        (real(b1,DP)-0.5d0)*d1, c, au2si(s(1,b1,b2)), &
                        au2si(s(2,b1,b2)), au2si(s(3,b1,b2))
                end do
                write(PFD, '(a)') ''
            end do
            close(PFD)
            call msg_out('Squares profile written to squares_profile'//trim(suffix)//'.dat')
        end if
        call nl

        deallocate(s, node)
    end subroutine

    ! Maxima, minima and zero crossings of the profile along the width, and
    ! the current integrated between consecutive zero crossings of the net
    ! current. The extrema are located with a parabola through the slice and
    ! its two neighbours, the zero crossings by linear interpolation between
    ! slice centres; the integrated ranges consist of whole slices.
    subroutine analyse_profile(g, s, d2, from_in, suffix)
        type(grid_t), intent(in) :: g
        real(DP), dimension(:,:,:), intent(in) :: s
        real(DP), intent(in) :: d2
        logical, intent(in) :: from_in
        character(*), intent(in) :: suffix

        integer(I4) :: nb2, n, q
        real(DP), dimension(:,:), allocatable :: v, sl
        real(DP), dimension(3) :: p, tot
        real(DP) :: a, b, c, off, val, x0, xprev
        character(10), dimension(3), parameter :: qname = &
            (/ 'net       ', 'diatropic ', 'paratropic' /)
        character(7) :: kind
        real(DP), parameter :: THRS=1.d-3   ! nA/T per bohr: ignore extrema of the far-field noise

        nb2=size(s,3)
        allocate(v(3,nb2), sl(3,nb2))
        do n=1,nb2
            do q=1,3
                sl(q,n)=au2si(sum(s(q,:,slice_index(n,nb2,from_in))))
            end do
            v(:,n)=sl(:,n)/d2
        end do

        open(PFD, file='profile_points'//trim(suffix)//'.dat')
        write(PFD, '(a)') '# Analysis of the current profile along the width; x in bohr from the start of the profile'
        write(PFD, '(a)') '# extremum  quantity        x            cx            cy            cz     dJ/dx [nA/T/bohr]'
        call nl
        call msg_out('Current profile: extrema of dJ/dx [nA/T/bohr]')
        do q=1,3
            do n=2,nb2-1
                a=v(q,n-1); b=v(q,n); c=v(q,n+1)
                if (.not.((b > a .and. b >= c) .or. (b < a .and. b <= c))) cycle
                if (abs(b) < THRS) cycle
                off=0.d0
                val=b
                if (abs(a-2.d0*b+c) > 1.d-14) then
                    off=0.5d0*d2*(a-c)/(a-2.d0*b+c)
                    val=b-(a-c)**2/(8.d0*(a-2.d0*b+c))
                end if
                x0=(real(n,DP)-0.5d0)*d2+off
                p=point_at(g, x0, from_in)
                if (b > a) then
                    kind='maximum'
                else
                    kind='minimum'
                end if
                write(PFD, '(2x,a7,2x,a10,f12.4,3f14.6,f16.6)') kind, qname(q), x0, p, val
                write(str_g, '(2x,a10,1x,a7,f9.3,a,3f8.3,a,f11.4)') qname(q), kind, x0, &
                    ' at (', p, ' )', val
                call msg_out(str_g)
            end do
        end do

        write(PFD, '(a)') '#'
        write(PFD, '(a)') '# zero crossings of the net current, and the current [nA/T] integrated between them (whole slices)'
        write(PFD, '(a)') '#    x_from      x_to           net           dia          para     [zero crossing at x_to: cx cy cz]'
        call nl
        call msg_out('Current profile: zero crossings of the net current and the current between them [nA/T]')
        tot=0.d0
        xprev=0.d0
        do n=1,nb2
            tot=tot+sl(:,n)
            if (n < nb2) then
                if (v(1,n)*v(1,n+1) < 0.d0) then
                    x0=(real(n,DP)-0.5d0)*d2+d2*v(1,n)/(v(1,n)-v(1,n+1))
                    p=point_at(g, x0, from_in)
                    write(PFD, '(2f11.4,3f14.6,4x,3f14.6)') xprev, x0, tot, p
                    write(str_g, '(2x,f8.3,a,f8.3,a,3f12.5)') xprev, ' to', x0, ':', tot
                    call msg_out(str_g)
                    xprev=x0
                    tot=0.d0
                end if
            end if
        end do
        x0=real(nb2,DP)*d2
        write(PFD, '(2f11.4,3f14.6)') xprev, x0, tot
        write(str_g, '(2x,f8.3,a,f8.3,a,3f12.5)') xprev, ' to', x0, ':', tot
        call msg_out(str_g)
        close(PFD)
        call msg_out('Profile analysis written to profile_points'//trim(suffix)//'.dat')
        deallocate(v, sl)
    end subroutine

    ! point on the mid-height line of the plane at the profile coordinate x
    function point_at(g, x, from_in) result(c)
        type(grid_t), intent(in) :: g
        real(DP), intent(in) :: x
        logical, intent(in) :: from_in
        real(DP), dimension(3) :: c

        if (from_in) then
            c=slice_centre(g, 0.5d0*g%l(1), g%l(2)-x)
        else
            c=slice_centre(g, 0.5d0*g%l(1), x)
        end if
    end function

    ! block index along the width of the n:th slice in writing order
    function slice_index(n, nb2, from_in) result(b2)
        integer(I4), intent(in) :: n, nb2
        logical, intent(in) :: from_in
        integer(I4) :: b2

        if (from_in) then
            b2=nb2+1-n
        else
            b2=n
        end if
    end function

    ! point of the plane at the distances s along basv1 and t along basv2
    ! from the grid origin (and in the middle of the third axis)
    function slice_centre(g, s, t) result(c)
        type(grid_t), intent(in) :: g
        real(DP), intent(in) :: s, t
        real(DP), dimension(3) :: c

        c=g%origin+s*g%basv(:,1)+t*g%basv(:,2)+0.5d0*g%l(3)*g%basv(:,3)
    end function

end module

! vim:et:sw=4:ts=4
