!
! Basis function evaluator
!

module caos_module
    use globals_module
    use settings_module
    use gtodefs_module
    use basis_class
    use cao2sao_class
    implicit none

    public  cgto, dcgto, cgto_dr
    private
contains

    ! Evaluate one contracted CAO and its gradient in one go: the
    ! exponentials are computed once (cgto + 3 x dcgto compute them four
    ! times) and the Cartesian powers by multiplication rather than pow().
    ! val(i) and dval(i,1:3) equal what cgto and dcgto return.
    subroutine cgto_dr(r, ctr, val, dval)
        real(DP), dimension(3), intent(in) :: r
        type(contraction_t), intent(in), target :: ctr
        real(DP), dimension(:), intent(out) :: val
        real(DP), dimension(:,:), intent(out) :: dval

        real(DP), dimension(:,:), pointer :: f
        real(DP), dimension(0:MAX_L,3) :: rp   ! rp(n,ax) = r(ax)**n
        real(DP) :: bfval, dbfval, p, q
        integer(I4), dimension(3) :: n, m
        integer(I4) :: i, ax, l

        l=ctr%l
        call get_gto_nlm(l, f)
        call cao2(ctr, sum(r**2), bfval, dbfval)

        rp(0,:)=D1
        do i=1,l
            rp(i,:)=rp(i-1,:)*r
        end do

        do i=1,ctr%nccomp
            n=nint(f(:,i))
            p=rp(n(1),1)*rp(n(2),2)*rp(n(3),3)
            val(i)=p*bfval
            ! d/dr_ax = n_ax r_ax**(n_ax-1) (...) bfval - 2 r_ax r**n dbfval
            do ax=1,3
                if (n(ax) > 0) then
                    m=n
                    m(ax)=m(ax)-1
                    q=rp(m(1),1)*rp(m(2),2)*rp(m(3),3)
                    dval(i,ax)=n(ax)*q*bfval-2.d0*r(ax)*p*dbfval
                else
                    dval(i,ax)=-2.d0*r(ax)*p*dbfval
                end if
            end do
        end do
    end subroutine

    subroutine cgto(r, ctr, val)
        real(DP), dimension(:), intent(in) :: r
        type(contraction_t), intent(in), target :: ctr
        real(DP), dimension(:), intent(out) :: val

        real(DP) :: p
        real(DP) :: q
        integer(I4) :: i
        real(DP), dimension(:,:), pointer :: f
        real(DP) :: rr2

        rr2=sum(r**2)

        call get_gto_nlm(ctr%l, f)
        q=cao(ctr, rr2)
        do i=1,ctr%nccomp
            p=product(r**f(:,i))*q
!            if (abs(p) > 1.d-20) val(i)=p
            val(i)=p
        end do
    end subroutine

    subroutine dcgto(r, ctr, ax, val)
        real(DP), dimension(:), intent(in) :: r
        type(contraction_t), intent(in), target :: ctr
        integer(I4), intent(in) :: ax
        real(DP), dimension(:), intent(out) :: val

        real(DP) :: bfval, dbfval
        real(DP), dimension(:,:), pointer :: f
        real(DP), dimension(3) :: df
        real(DP) :: up, down
        integer(I4) :: i, j
        real(DP) :: rr2

        rr2=sum(r**2)

        call get_gto_nlm(ctr%l,f)
        call cao2(ctr, rr2, bfval, dbfval)
        do i=1,ctr%nccomp
            df=f(:,i)
            df(ax)=df(ax)-1.d0
            if (df(ax) < D0) df(ax)=D0
            down=f(ax,i)*product(r**df)*bfval
            up=2.d0*r(ax)*product(r**f(:,i))*dbfval
            val(i)=down-up
        end do
    end subroutine

    ! Evaluate one contracted CAO
    function cao(cc, rr2) result(ff)
        type(contraction_t), intent(in)  :: cc
        real(DP), intent(in) :: rr2
        real(DP) :: ff

        integer(I4) :: i

        ff=D0
        do i=1,cc%npf
            ff=ff+cc%ncc(i)*exp(-cc%xp(i)*rr2)
        end do
    end function

    ! Evaluate one differentiated CAO
    function dcao(cc, rr2) result(ff)
        type(contraction_t), intent(in) :: cc
        real(DP), intent(in) :: rr2
        real(DP) :: ff

        integer(I4) :: i

        ff=D0
        do i=1,cc%npf
            ff=ff-cc%xp(i)*cc%ncc(i)*exp(-cc%xp(i)*rr2)
        end do
    end function

    subroutine cao2(cc, rr2, vcao, vdcao)
        type(contraction_t), intent(in)  :: cc
        real(DP), intent(in) :: rr2
        real(DP), intent(out) :: vcao, vdcao

        integer(I4) :: i
        real(DP) :: q

        vcao=D0
        vdcao=D0
        do i=1,cc%npf
            q=cc%ncc(i)*exp(-cc%xp(i)*rr2)
!            if (abs(q) < 1.d-20) cycle
            vcao=vcao+q
            vdcao=vdcao+cc%xp(i)*q
        end do
    end subroutine

end module

! vim:et:sw=4:ts=4
