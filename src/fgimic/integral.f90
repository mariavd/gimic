! Integrate the current density using quadrature.
!
! Written by Jonas Juselius, University of Helsinki, 2003 (I think).
! Improvements by jj, University of Tromsø, 2011 (I'm sure).
!

module integral_class
    use globals_module
    use settings_module
    use grid_class
    use jfield_class
    use dens_class
    use jtensor_class
    use gaussint_module
    use lip_module
    use teletype_module
    use parallel_module
    use magnet_module
    use profile_module
    ! ACID stuff
    use acid_module
    use tensor_module
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none

    type integral_t
        type(grid_t), pointer :: grid
    end type

    public new_integral, del_integral, integral_t
    public integrate_tensor_field, integrate_current
    public integrate_modulus
    public integrate_acid

    private

    character(8) :: spin='total'
contains
    subroutine new_integral(this, grid)
        type(integral_t) :: this
        type(grid_t), target :: grid

        this%grid=>grid
    end subroutine

    subroutine del_integral(this)
        type(integral_t) :: this
        nullify(this%grid)
    end subroutine

    ! The direction along which the current through the plane is reported.
    ! Without a reference centre it is the grid normal basv3 (for bond grids
    ! the direction from atom 1 to atom 2). With Grid.towards the normal is
    ! chosen such that a current flowing in the diatropic sense around the
    ! reference centre C is positive: the diatropic flow at a point P is
    ! along B x (P - C), so the sign of (B x (P - C)).basv3 decides. A plane
    ! that reaches past C sees that flow from both sides, so P is taken at
    ! the end of the profile: the "out" edge of a bond grid (the grid
    ! origin side), the edge far from the origin for other grids. Only the
    ! sign is taken from this, nothing else changes.
    function reporting_normal(g, bb) result(n)
        type(grid_t) :: g
        real(DP), dimension(3), intent(in) :: bb
        real(DP), dimension(3) :: n

        real(DP), dimension(3) :: c, r, t
        real(DP) :: x, t_end
        logical, save :: announced=.false.

        n=get_grid_normal(g)
        if (.not.g%has_towards) return

        if (trim(g%mode) == 'bond') then
            t_end=0.d0
        else
            t_end=g%l(2)
        end if
        c=g%origin+0.5d0*g%l(1)*g%basv(:,1)+t_end*g%basv(:,2)+0.5d0*g%l(3)*g%basv(:,3)
        r=c-g%towards
        t=cross_product(bb, r)
        x=dot_product(t, n)
        if (abs(x) < 1.d-8*sqrt(sum(bb**2))*sqrt(sum(r**2))) then
            call msg_warn('towards: the reference centre lies along the field &
                &from the end of the plane; sign convention undefined, &
                &keeping the grid normal')
            return
        end if
        if (x < 0.d0) n=-n
        if (.not.announced) then
            if (x < 0.d0) then
                call msg_info('Plane normal reversed: currents are positive &
                    &when they flow in the diatropic sense around the &
                    &reference centre')
            else
                call msg_info('Currents are positive when they flow in the &
                    &diatropic sense around the reference centre')
            end if
            announced=.true.
        end if
    end function

    subroutine integrate_current(this, mol, xdens, spinn)
        type(integral_t), intent(inout) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spinn

        integer(I4) :: i, j, k, p1, p2, p3, lo, hi, nb1, b, order, ndone, nstep
        real(DP), dimension(3) :: normal, rr, center, bb, tot
        real(DP) :: w, jp, r, bound
        real(DP), dimension(3) :: jvec
        real(DP), dimension(9) :: tt
        real(DP), dimension(:,:), allocatable :: coords, ttb, rsum
        ! blk(:,b,j): current through column j of the grid, resolved into the
        ! blocks b of the height axis (net, positive, negative contribution),
        ! integrated along the height and the third axis but not yet
        ! weighted along the width. Each column is written by exactly one
        ! thread, so no reduction is needed, the totals are summed in a
        ! fixed order afterwards (independent of the number of threads), and
        ! the current profile is a sum over blocks of columns.
        real(DP), dimension(:,:,:), allocatable :: blk
        type(jtensor_t) :: jt
        integer(8) :: t0, t1, rate

        if (present(spinn)) spin = spinn

        if (settings%is_uhf) then
            select case(spin)
                case('total')
                    call msg_note("Integrating total density")
                case('alpha')
                    call msg_note("Integrating alpha density")
                case('beta')
                    call msg_note("Integrating beta density")
                case('spindens')
                    call msg_note("Integrating spin density")
                case default
                    call msg_error("Invalid spin: " // spin)
                    stop
            end select
        end if

        call get_grid_size(this%grid, p1, p2, p3)
        call get_magnet(this%grid, bb)

        normal=reporting_normal(this%grid, bb)

        bound=this%grid%radius
        if (bound < 1.d+10) then
            write(str_g, *) 'Integration bound set to radius ', bound
            call msg_out(str_g)
        end if

        call grid_center(this%grid,center)
        call schedule(p2, lo, hi)
        if (mpi_rank == 0) then
            lo = 1
        end if

        ! blocks along the height: one per slice, or a single one
        nb1=max(1, this%grid%nslice(1))
        order=this%grid%gauss_order
        allocate(blk(3,nb1,p2))
        blk=0.d0
        ndone=0
        nstep=max(1, (p2*p3)/10)
        call system_clock(t0, rate)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(i,j,k,b,r,rr,w,jp,tt,jvec,jt,coords,ttb,rsum,t1) &
!$OMP SHARED(p1,p2,p3,this,center,spin,bb,normal,mol,xdens,lo,hi,bound) &
!$OMP SHARED(nb1,order,blk,ndone,nstep,t0,rate,mpi_rank)
        call new_jtensor(jt, mol, xdens)
        allocate(coords(3,p1), ttb(9,p1), rsum(3,nb1))
        do k=1,p3
            ! columns far from the molecule are cheap (screening): balance
            !$OMP DO SCHEDULE(DYNAMIC,4)
            do j=lo,hi
                rsum=0.d0
                do i=1,p1
                    coords(:,i)=gridpoint(this%grid, i, j, k)
                end do
                call ctensor_batch(jt, coords, ttb, spin)
                do i=1,p1
                    rr=coords(:,i)
                    r=sqrt(sum((rr-center)**2))
                    tt=ttb(:,i)
                    jvec=matmul(reshape(tt,(/3,3/)),bb)
                    if ( r > bound ) then
                        jp=0.d0
                    else
                        w=get_weight(this%grid, i, 1)
                        jp=dot_product(normal,jvec)*w
                    end if
                    b=1
                    if (nb1 > 1) b=(i-1)/order+1
                    rsum(1,b)=rsum(1,b)+jp
                    if (jp > 0.d0) then
                        ! positive contribution
                        rsum(2,b)=rsum(2,b)+jp
                    else
                        ! negative contribution
                        rsum(3,b)=rsum(3,b)+jp
                    end if
                end do
                w=get_weight(this%grid,k,3)
                blk(:,:,j)=blk(:,:,j)+w*rsum

                !$OMP ATOMIC CAPTURE
                ndone=ndone+1
                i=ndone
                !$OMP END ATOMIC
                ! progress report on stderr (stdout holds the results)
                if (p2*p3 >= 100 .and. mod(i,nstep) == 0 .and. mpi_rank == 0) then
                    !$OMP CRITICAL(progress)
                    call system_clock(t1)
                    write(error_unit, '(a,i8,a,i8,a,f10.1,a)') ' columns done:', i, &
                        ' of', p2*p3, '   wall time', &
                        real(t1-t0,DP)/real(rate,DP), ' s'
                    flush(error_unit)
                    !$OMP END CRITICAL(progress)
                end if
            end do
            !$OMP END DO
        end do
        deallocate(coords, ttb, rsum)
        call del_jtensor(jt)
!$OMP END PARALLEL

        ! MPI: every rank filled its own columns, sum them onto rank 0
        call collect_array(blk, size(blk))

        ! totals, summed over the columns in a fixed order
        tot=0.d0
        do j=1,p2
            w=get_weight(this%grid,j,2)
            tot=tot+w*sum(blk(:,:,j),dim=2)
        end do

        call nl
        call msg_out(repeat('*', 60))
        write(str_g, '(a,f13.6)') '   Induced current (au)    :', tot(1)
        call msg_out(str_g)
        write(str_g, '(a,f13.6,a,f11.6,a)') &
            '      Positive contribution:', tot(2), '  (',au2si(tot(2)),' )'
        call msg_out(str_g)
        write(str_g, '(a,f13.6,a,f11.6,a)') &
            '      Negative contribution:', tot(3), '  (',au2si(tot(3)),' )'
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f13.6)') '   Induced current (nA/T)  :', au2si(tot(1))
        call msg_out(str_g)
        write(str_g, '(a,f13.6)') '      (conversion factor)  :', au2si(1.d0)
        call msg_out(str_g)
        call msg_out(repeat('*', 60))
        call nl

        if (mpi_rank == 0 .and. any(this%grid%nslice > 0)) then
            call write_profile(this%grid, blk, normal, bb, spin)
        end if
        deallocate(blk)

        spin = 'total'
    end subroutine

    ! integrate the modulus of the current, retaining the sign
    ! test version
    subroutine integrate_modulus(this, mol, xdens, spinn)
        type(integral_t), intent(inout) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spinn

        integer(I4) :: i, j, k, p1, p2, p3, lo, hi
        real(DP), dimension(3) :: normal, rr, center, bb
        real(DP) :: psum, nsum, w, jp, r, bound, sgn
        real(DP) :: psum2, nsum2
        real(DP) :: psum3, nsum3
        real(DP) :: xsum, xsum2, xsum3
        real(DP), dimension(3) :: jvec
        real(DP), dimension(9) :: tt
        real(DP), dimension(:,:), allocatable :: coords, ttb
        type(jtensor_t) :: jt

        if (present(spinn)) spin = spinn

        if (settings%is_uhf) then
            select case(spin)
                case('total')
                    call msg_note("Integrating total density")
                case('alpha')
                    call msg_note("Integrating alpha density")
                case('beta')
                    call msg_note("Integrating beta density")
                case('spindens')
                    call msg_note("Integrating spin density")
                case default
                    call msg_error("Invalid spin: " // spin)
                    stop
            end select
        end if

        call get_grid_size(this%grid, p1, p2, p3)
        call get_magnet(this%grid, bb)

        normal=reporting_normal(this%grid, bb)

        bound=this%grid%radius
        if (bound < 1.d+10) then
            write(str_g, *) 'Integration bound set to radius ', bound
            call msg_out(str_g)
        end if
        call grid_center(this%grid,center)
        call schedule(p2, lo, hi)
        if (mpi_rank == 0) then
            lo = 1
        end if

        xsum3=0.d0
        psum3=0.d0
        nsum3=0.d0

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(i,j,k,r,rr,sgn,xsum,psum,nsum,xsum2,psum2,nsum2) &
!$OMP PRIVATE(jt,w,jp,tt,jvec,coords,ttb) &
!$OMP SHARED(p1,p2,p3,this,center,spin,bb,normal,lo,hi,bound,xdens,mol) &
!$OMP REDUCTION(+:xsum3,psum3,nsum3)
        call new_jtensor(jt, mol, xdens)
        allocate(coords(3,p1), ttb(9,p1))
        do k=1,p3
            xsum2=0.d0
            psum2=0.d0
            nsum2=0.d0

            !$OMP DO
            do j=lo,hi
                xsum=0.d0
                psum=0.d0
                nsum=0.d0
                do i=1,p1
                    coords(:,i)=gridpoint(this%grid, i, j, k)
                end do
                call ctensor_batch(jt, coords, ttb, spin)
                do i=1,p1
                    rr=coords(:,i)
                    r=sqrt(sum((rr-center)**2))
                    ! points outside the integration radius do not contribute;
                    ! skip them before touching their tensor
                    if (r > bound) cycle
                    tt=ttb(:,i)
                    jvec=matmul(reshape(tt,(/3,3/)),bb)
                    w=get_weight(this%grid, i, 1)
                    jp=dot_product(normal,jvec)
                    if (abs(jp) < 1.d-12) then ! prob. parallel component
                        sgn=0.d0
                    else if (jp > 0) then
                        sgn=1.d0
                    else
                        sgn=-1.d0
                    end if
                    jp=sgn*sqrt(sum(jvec**2))
                    xsum=xsum+jp*w
                    if (jp > 0.d0) then
                        psum=psum+jp*w
                    else
                        nsum=nsum+jp*w
                    end if
                end do
                w=get_weight(this%grid,j,2)
                xsum2=xsum2+xsum*w
                psum2=psum2+psum*w
                nsum2=nsum2+nsum*w
            end do
            !$OMP END DO

            w = get_weight(this%grid,k,3)
            xsum3 = xsum3 + w*xsum2
            psum3 = psum3 + w*psum2
            nsum3 = nsum3 + w*nsum2

            ! TODO old code used collect_sum
            ! this will be needed for MPI

        end do
        deallocate(coords, ttb)
        call del_jtensor(jt)
!$OMP END PARALLEL

        call nl
        call msg_out(repeat('*', 60))
        write(str_g, '(a,f13.6)') 'Induced mod current (au)   :', xsum3
        call msg_out(str_g)
        write(str_g, '(a,f13.6,a,f11.6,a)') &
            '      Positive contribution:', psum3, '  (',au2si(psum3),' )'
        call msg_out(str_g)
        write(str_g, '(a,f13.6,a,f11.6,a)') &
            '      Negative contribution:', nsum3, '  (',au2si(nsum3),' )'
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f13.6)') 'Induced mod current (nA/T) :', au2si(xsum3)
        call msg_out(str_g)
        write(str_g, '(a,f13.6)') '      (conversion factor)  :', au2si(1.d0)
        call msg_out(str_g)
        call msg_out(repeat('*', 60))
        call nl
        spin = 'total'
    end subroutine

    subroutine integrate_tensor_field(this, mol, xdens)
        type(integral_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens

        integer(I4) :: i, j, k, p1, p2, p3
        real(DP), dimension(9) :: xsum
        real(DP), dimension(:,:), allocatable  :: jt1, jt2, jt3, coords
        type(jtensor_t) :: jt

        !call jfield_eta(this%jf)
        call get_grid_size(this%grid, p1, p2, p3)

        call new_jtensor(jt, mol, xdens)
        allocate(jt1(9,p1))
        allocate(jt2(9,p2))
        allocate(jt3(9,p3))
        allocate(coords(3,p1))

        do k=1,p3
            do j=1,p2
                do i=1,p1
                    coords(:,i)=gridpoint(this%grid, i, j, k)
                end do
                call ctensor_batch(jt, coords, jt1, 'total')
                jt2(:,j)=int_t_1d(jt1,this%grid,1)
            end do
            jt3(:,k)=int_t_1d(jt2,this%grid,2)
        end do
        xsum=int_t_1d(jt3,this%grid,3)

        call print_tensor_int(xsum)
        call del_jtensor(jt)
        deallocate(jt1, jt2, jt3, coords)
    end subroutine

    function int_t_1d(jt, grid, axis) result(xsum)
        real(DP), dimension(:,:), intent(in) :: jt
        type(grid_t), intent(in) :: grid
        integer(I4), intent(in) :: axis
        real(DP), dimension(9) :: xsum

        integer(I4) :: i, pts
        real(DP) :: w

        pts=size(jt,2)

        xsum=0.d0
        do i=1,pts
            w=get_weight(grid, i, axis)
            xsum=xsum+jt(:,i)*w
        end do
    end function

    subroutine print_tensor_int(xsum)
        real(DP), dimension(3,3) :: xsum

        integer(I4) :: i
        print *, xsum(1,:)
        print *, xsum(2,:)
        print *, xsum(3,:)

        call nl
        call msg_out(repeat('*', 70))
        write(str_g, '(a,3f13.6)') '   Induced current (au)    :', xsum(1,:)
        call msg_out(str_g)
        write(str_g, '(a,3f13.6)') '                            ', xsum(2,:)
        call msg_out(str_g)
        write(str_g, '(a,3f13.6)') '                            ', xsum(3,:)
        call msg_out(str_g)
        call nl
        write(str_g, '(a,3f13.6)') '   Induced current (nA/T)  :', &
            (au2si(xsum(1,i)), i=1,3)
        call msg_out(str_g)
        write(str_g, '(a,3f13.6)') '                            ', &
            (au2si(xsum(2,i)), i=1,3)
        call msg_out(str_g)
        write(str_g, '(a,3f13.6)') '                            ', &
            (au2si(xsum(3,i)), i=1,3)
        call msg_out(str_g)
        call nl
        write(str_g, '(a,f13.6)') '        Conversion Factor   :', au2si(1.d0)
        call msg_out(str_g)
        call msg_out(repeat('*', 70))
        call nl
    end subroutine

    ! do here ACID integration here !
    subroutine integrate_acid(this, mol, xdens)
    use acid_module
    ! this subroutine is based on integrate_current !
        type(integral_t), intent(inout) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens

        integer(I4) :: i, j, k, p1, p2, p3, lo, hi
        real(DP), dimension(3) :: rr, center
        real(DP) :: w, r, bound
        real(DP) :: xsum, xsum2, xsum3
        real(DP) :: val, acid_val
        real(DP), dimension(9) :: tt
        real(DP), dimension(:,:), allocatable :: coords, ttb
        type(jtensor_t) :: jt

        if (settings%is_uhf) then
            select case(spin)
                case('total')
                    call msg_note("Integrating total density")
                case('alpha')
                    call msg_note("Integrating alpha density")
                case('beta')
                    call msg_note("Integrating beta density")
                case('spindens')
                    call msg_note("Integrating spin density")
                case default
                    call msg_error("Invalid spin: " // spin)
                    stop
            end select
        end if

        call get_grid_size(this%grid, p1, p2, p3)

        bound=1.d+10
        bound=this%grid%radius
        if (bound < 1.d+10) then
            write(str_g, *) 'Integration bound set to radius ', bound
            call msg_out(str_g)
        end if

        call grid_center(this%grid,center)
        call schedule(p2, lo, hi)
        if (mpi_rank == 0) then
            lo = 1
        end if

        xsum3=0.d0

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(i,j,k,r,rr,xsum,xsum2) &
!$OMP PRIVATE(jt,w,tt,val,coords,ttb) &
!$OMP SHARED(p1,p2,p3,this,center,spin,mol,xdens,lo,hi,bound) &
!$OMP REDUCTION(+:xsum3)
        call new_jtensor(jt, mol, xdens)
        allocate(coords(3,p1), ttb(9,p1))
        do k=1,p3
            xsum2=0.d0

            !$OMP DO
            do j=lo,hi
                xsum=0.d0
                do i=1,p1
                    coords(:,i)=gridpoint(this%grid, i, j, k)
                end do
                call ctensor_batch(jt, coords, ttb, spin)
                do i=1,p1
                    rr=coords(:,i)
                    r=sqrt(sum((rr-center)**2))
                    tt=ttb(:,i)
                    ! attention: output of get_acid is in au !
                    val = get_acid(tt)
                    if ( r > bound ) then
                        w=0.d0
                    else
                        w=get_weight(this%grid, i, 1)
                    end if
                    xsum=xsum+val*w
                end do
                w=get_weight(this%grid,j,2)
                xsum2=xsum2+xsum*w
            end do
            !$OMP END DO

            w = get_weight(this%grid,k,3)
            xsum3 = xsum3 + w*xsum2

            ! TODO old code used collect_sum
            ! this will be needed for MPI
        end do
        deallocate(coords, ttb)
        call del_jtensor(jt)
!$OMP END PARALLEL
        acid_val = dsqrt(xsum3)

        call nl
        call msg_out(repeat('*', 60))
        write(str_g, '(a,f13.6)') '   ACID (au) sqrt(delta J^2):', acid_val
        call msg_out(str_g)
        write(str_g, '(a,f13.6)') '   ACID (nA/T)              :', au2si(acid_val)
        call msg_out(str_g)
        call nl
        call msg_out(repeat('*', 60))
        call nl
    end subroutine

end module

! vim:et:sw=4:ts=4
