!
! Tropicity of the current density at a point, by following the streamline
! through that point.
!
! The tropicity is a non-local property: it is the sense in which the
! streamline through a point circulates with respect to the magnetic field.
! Here the streamline is integrated with an adaptive Runge-Kutta method,
! the current density being evaluated analytically at every stage (no
! interpolation on a grid). A streamline is classified as soon as it has
! made one turn, which is detected in two ways:
!
!   * it comes back to its starting point. The sense of circulation is the
!     sign of the area it encloses, projected on the magnetic field,
!     A = 1/2 sum (r - r0) x dr . B
!   * the direction of the current, projected on the plane perpendicular to
!     the field, has turned by more than a full turn. This classifies
!     streamlines that spiral and never close.
!
! A streamline that does neither within the maximum path length, that runs
! into a region of vanishing current, or that encloses no area with respect
! to the field (a loop in a plane containing the field) is left undefined.
!
! Convention: a diatropic current circulates in the right-handed sense
! around the magnetic field (A > 0), i.e. clockwise when looking along the
! field, and a paratropic current in the left-handed sense (A < 0). This is
! the sign of (B x r).J for a current circulating around the origin, as
! used for the signed modulus of the current density (jmod.vti).
!
module tropicity_module
    use globals_module
    use settings_module
    use jtensor_class
    use teletype_module
    implicit none

    integer(I4), parameter :: TROP_UNDEF = 0
    integer(I4), parameter :: TROP_DIA = 1
    integer(I4), parameter :: TROP_PARA = -1

    real(DP), parameter :: PI = 3.141592653589793d0
    ! currents weaker than this (au) have no direction worth following
    real(DP), parameter :: JMIN = 1.d-12
    ! a closed streamline must enclose at least this fraction of the area
    ! of a circle of the same perimeter, projected on the field
    real(DP), parameter :: FLAT = 1.d-2
    ! the projected direction of the current is followed only where it is
    ! not nearly parallel to the field
    real(DP), parameter :: PMIN = 1.d-1
    ! turn of the projected direction that classifies an open streamline
    real(DP), parameter :: FULL_TURN = 2.d0*PI*1.02d0

    ! smallest step, as a fraction of the largest one
    real(DP), parameter :: HMIN_FAC = 1.d-3
    ! steps (accepted or not) after which a streamline is given up
    integer(I4), parameter :: MAXITER = 20000

    public tropicity_batch, tropicity_summary
    public TROP_UNDEF, TROP_DIA, TROP_PARA
    private

contains
    ! Classify the points r0(3,n). jt must be private to the calling thread.
    !
    ! All streamlines of the batch are advanced together, so that the
    ! current density is evaluated for many points at a time. Each has its
    ! own step, chosen by the embedded Dormand-Prince 5(4) pair to keep the
    ! estimated error of a step below settings%trop_accuracy.
    subroutine tropicity_batch(jt, r0, bb, spin, label)
        type(jtensor_t) :: jt
        real(DP), dimension(:,:), intent(in) :: r0
        real(DP), dimension(3), intent(in) :: bb
        character(*), intent(in) :: spin
        integer(I4), dimension(:), intent(out) :: label

        ! Dormand-Prince coefficients
        real(DP), parameter :: A21=1.d0/5.d0
        real(DP), parameter :: A31=3.d0/40.d0, A32=9.d0/40.d0
        real(DP), parameter :: A41=44.d0/45.d0, A42=-56.d0/15.d0, &
            A43=32.d0/9.d0
        real(DP), parameter :: A51=19372.d0/6561.d0, A52=-25360.d0/2187.d0, &
            A53=64448.d0/6561.d0, A54=-212.d0/729.d0
        real(DP), parameter :: A61=9017.d0/3168.d0, A62=-355.d0/33.d0, &
            A63=46732.d0/5247.d0, A64=49.d0/176.d0, A65=-5103.d0/18656.d0
        real(DP), parameter :: B1=35.d0/384.d0, B3=500.d0/1113.d0, &
            B4=125.d0/192.d0, B5=-2187.d0/6784.d0, B6=11.d0/84.d0
        real(DP), parameter :: E1=71.d0/57600.d0, E3=-71.d0/16695.d0, &
            E4=71.d0/1920.d0, E5=-17253.d0/339200.d0, E6=22.d0/525.d0, &
            E7=-1.d0/40.d0

        integer(I4) :: n, na, ia, ip, iter, st
        real(DP) :: hmax, hmin, tol, etol, maxlen, d, pn, err, hh
        real(DP), dimension(3) :: b, dr, rel, p, q
        logical :: done, stuck
        integer(I4), dimension(:), allocatable :: idx
        real(DP), dimension(:), allocatable :: area, turn, h, plen
        logical, dimension(:), allocatable :: away, hasp
        real(DP), dimension(:,:), allocatable :: r, rs, pold, tt
        real(DP), dimension(:,:,:), allocatable :: k

        n=size(r0,2)
        label=TROP_UNDEF
        if (n < 1) return

        hmax=settings%trop_step
        hmin=HMIN_FAC*hmax
        tol=settings%trop_tol
        etol=settings%trop_accuracy
        maxlen=settings%trop_length
        b=bb/sqrt(sum(bb**2))

        allocate(idx(n), area(n), turn(n), h(n), plen(n), away(n), hasp(n))
        allocate(r(3,n), rs(3,n), pold(3,n), tt(9,n), k(3,7,n))

        na=n
        do ip=1,n
            idx(ip)=ip
        end do
        r=r0
        area=0.d0
        turn=0.d0
        plen=0.d0
        h=min(hmax, 0.1d0)
        away=.false.
        hasp=.false.

        call direction(r(:,1:na), 1)

        do iter=1,MAXITER
            if (na < 1) exit
            do ip=1,na
                rs(:,ip)=r(:,ip)+h(ip)*A21*k(:,1,ip)
            end do
            call direction(rs(:,1:na), 2)
            do ip=1,na
                rs(:,ip)=r(:,ip)+h(ip)*(A31*k(:,1,ip)+A32*k(:,2,ip))
            end do
            call direction(rs(:,1:na), 3)
            do ip=1,na
                rs(:,ip)=r(:,ip)+h(ip)*(A41*k(:,1,ip)+A42*k(:,2,ip)+ &
                    A43*k(:,3,ip))
            end do
            call direction(rs(:,1:na), 4)
            do ip=1,na
                rs(:,ip)=r(:,ip)+h(ip)*(A51*k(:,1,ip)+A52*k(:,2,ip)+ &
                    A53*k(:,3,ip)+A54*k(:,4,ip))
            end do
            call direction(rs(:,1:na), 5)
            do ip=1,na
                rs(:,ip)=r(:,ip)+h(ip)*(A61*k(:,1,ip)+A62*k(:,2,ip)+ &
                    A63*k(:,3,ip)+A64*k(:,4,ip)+A65*k(:,5,ip))
            end do
            call direction(rs(:,1:na), 6)
            do ip=1,na
                rs(:,ip)=r(:,ip)+h(ip)*(B1*k(:,1,ip)+B3*k(:,3,ip)+ &
                    B4*k(:,4,ip)+B5*k(:,5,ip)+B6*k(:,6,ip))
            end do
            ! the last stage is the first one of the next step
            call direction(rs(:,1:na), 7)

            ia=0
            do ip=1,na
                done=.false.
                ! a vanishing current at any stage: nothing to follow
                stuck=.false.
                hh=h(ip)
                err=0.d0
                do st=1,7
                    if (sum(k(:,st,ip)**2) < 0.5d0) stuck=.true.
                end do
                if (stuck) then
                    done=.true.
                else
                    err=hh*sqrt(sum((E1*k(:,1,ip)+E3*k(:,3,ip)+E4*k(:,4,ip)+ &
                        E5*k(:,5,ip)+E6*k(:,6,ip)+E7*k(:,7,ip))**2))
                    h(ip)=hh*min(5.d0, max(0.2d0, &
                        0.9d0*(etol/max(err, 1.d-30))**0.2d0))
                    h(ip)=min(hmax, max(hmin, h(ip)))
                end if

                if (.not.stuck .and. (err <= etol .or. hh <= hmin)) then
                    ! step accepted
                    dr=rs(:,ip)-r(:,ip)
                    rel=r(:,ip)-r0(:,idx(ip))
                    area(ip)=area(ip)+0.5d0*dot_product(cross(rel,dr), b)
                    plen(ip)=plen(ip)+hh

                    ! turn of the direction of the current around the field
                    p=k(:,7,ip)-dot_product(k(:,7,ip),b)*b
                    pn=sqrt(sum(p**2))
                    if (pn > PMIN) then
                        p=p/pn
                        if (hasp(ip)) then
                            q=pold(:,ip)
                            turn(ip)=turn(ip)+atan2( &
                                dot_product(cross(q,p), b), dot_product(q,p))
                        end if
                        pold(:,ip)=p
                        hasp(ip)=.true.
                    end if

                    ! back at the starting point?
                    d=segment_distance(rel, dr)
                    if (.not.away(ip)) then
                        if (sqrt(sum((rel+dr)**2)) > 2.d0*tol) away(ip)=.true.
                    else if (d < tol) then
                        if (abs(area(ip)) > FLAT*plen(ip)**2/(4.d0*PI)) then
                            label(idx(ip))=sense(area(ip))
                        end if
                        done=.true.
                    end if
                    if (.not.done .and. abs(turn(ip)) > FULL_TURN) then
                        label(idx(ip))=sense(turn(ip))
                        done=.true.
                    end if
                    if (plen(ip) > maxlen) done=.true.
                    r(:,ip)=rs(:,ip)
                    k(:,1,ip)=k(:,7,ip)
                end if

                ! keep the unfinished streamlines at the front
                if (.not.done) then
                    ia=ia+1
                    if (ia /= ip) then
                        idx(ia)=idx(ip)
                        r(:,ia)=r(:,ip)
                        k(:,1,ia)=k(:,1,ip)
                        h(ia)=h(ip)
                        plen(ia)=plen(ip)
                        area(ia)=area(ip)
                        turn(ia)=turn(ip)
                        away(ia)=away(ip)
                        hasp(ia)=hasp(ip)
                        pold(:,ia)=pold(:,ip)
                    end if
                end if
            end do
            na=ia
        end do

        deallocate(idx, area, turn, h, plen, away, hasp)
        deallocate(r, rs, pold, tt, k)

    contains
        ! stage st: unit vector along the current at the points x; zero
        ! where the current vanishes
        subroutine direction(x, st)
            real(DP), dimension(:,:), intent(in) :: x
            integer(I4), intent(in) :: st

            integer(I4) :: i, m
            real(DP) :: jn
            real(DP), dimension(3) :: jv

            m=size(x,2)
            call ctensor_batch(jt, x, tt(:,1:m), spin)
            do i=1,m
                jv=matmul(reshape(tt(:,i),(/3,3/)), bb)
                jn=sqrt(sum(jv**2))
                if (jn > JMIN) then
                    k(:,st,i)=jv/jn
                else
                    k(:,st,i)=0.d0
                end if
            end do
        end subroutine
    end subroutine

    function sense(x) result(l)
        real(DP), intent(in) :: x
        integer(I4) :: l

        if (x > 0.d0) then
            l=TROP_DIA
        else
            l=TROP_PARA
        end if
    end function

    ! distance from the origin to the segment a -> a+d
    function segment_distance(a, d) result(dist)
        real(DP), dimension(3), intent(in) :: a, d
        real(DP) :: dist

        real(DP) :: s, d2

        d2=sum(d**2)
        if (d2 > 0.d0) then
            s=-dot_product(a,d)/d2
            s=max(0.d0, min(1.d0, s))
        else
            s=0.d0
        end if
        dist=sqrt(sum((a+s*d)**2))
    end function

    function cross(a, b) result(c)
        real(DP), dimension(3), intent(in) :: a, b
        real(DP), dimension(3) :: c

        c(1)=a(2)*b(3)-a(3)*b(2)
        c(2)=a(3)*b(1)-a(1)*b(3)
        c(3)=a(1)*b(2)-a(2)*b(1)
    end function

    subroutine tropicity_summary(label)
        integer(I4), dimension(:), intent(in) :: label

        integer(I4) :: n

        n=size(label)
        if (n < 1) return
        write(str_g, '(a,i0,a)') 'Tropicity assigned by streamline tracing (', &
            n, ' points)'
        call msg_out(str_g)
        write(str_g, '(a,i10,a,f6.2,a)') '   diatropic :', &
            count(label == TROP_DIA), '  (', &
            100.d0*real(count(label == TROP_DIA), DP)/real(n, DP), ' %)'
        call msg_out(str_g)
        write(str_g, '(a,i10,a,f6.2,a)') '   paratropic:', &
            count(label == TROP_PARA), '  (', &
            100.d0*real(count(label == TROP_PARA), DP)/real(n, DP), ' %)'
        call msg_out(str_g)
        write(str_g, '(a,i10,a,f6.2,a)') '   undefined :', &
            count(label == TROP_UNDEF), '  (', &
            100.d0*real(count(label == TROP_UNDEF), DP)/real(n, DP), ' %)'
        call msg_out(str_g)
        call nl
    end subroutine
end module

! vim:et:sw=4:ts=4
