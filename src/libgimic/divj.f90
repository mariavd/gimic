!
!  Calculate the divergence of a vector field.
!  Useful for checking the convergce of the actual gauge invariance,
!  in contrast to the basis set convergence and gauge independence.
!

module divj_class
    use globals_m
    use grid_class
    use jfield_class
    use jtensor_class
    use teletype_m
    use magnet_m
    use parallel_m
    implicit none

    type divj_t
        real(DP), dimension(:,:), pointer :: buf
        type(grid_t), pointer :: grid
        type(jtensor_t), pointer :: jt
        real(DP), dimension(3) :: bb
    end type

    public divj_direct, divj_direct_plt, new_divj, del_divj
    public divj_plot, set_divj, divj_t, divj
    private
    
    real(DP), dimension(5), parameter :: wgt=(/2.d0,-16.d0,0.d0,16.d0,-2.d0/)
    real(DP), parameter :: step=1.d-3
    real(DP), parameter :: hx=1.d0/(24.d0*step)
    real(DP), parameter :: hy=1.d0/(24.d0*step)
    real(DP), parameter :: hz=1.d0/(24.d0*step)

contains
    subroutine new_divj(this, grid, jt)
        type(divj_t), intent(inout) :: this
        type(grid_t), target :: grid
        type(jtensor_t), target :: jt

        integer(I4) :: djrl, p1, p2
        logical :: foo_p

        this%bb=D0
        p1=0
        
        call push_section(input, 'divj')
        call get_magnet(grid, this%bb)
        call pop_section(input)

        call get_grid_size(grid, p1, p2)
        djrl=p1*p2*DP
        if (master_p) then
            open(DIVJFD, file='DIVJ', access='direct', recl=djrl)
        end if

        this%grid=>grid
        this%jt=>jt
        allocate(this%buf(p1,p2))
    end subroutine

    subroutine del_divj(this)
        type(divj_t), intent(inout) :: this

        deallocate(this%buf)
        if (master_p) then
            close(DIVJFD)
        end if
    end subroutine

    subroutine set_divj(this, k)
        type(divj_t), intent(in) :: this
        integer(I4), intent(in) :: k

        write(DIVJFD, rec=k) this%buf
    end subroutine

    subroutine divj_plot(this, gopen_file)
        type(divj_t), intent(inout) :: this
        character(*), intent(in) :: gopen_file
        
        integer(I4) :: i,j,p1,p2,p3
        real(DP) :: amax
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf

        call get_grid_size(this%grid, p1, p2, p3)

        buf=>this%buf

        amax=D0
        read(DIVJFD, rec=1) this%buf
        str_g='DIVJPLT'
        open(DJPFD, file=trim(str_g))
        do j=1,p2
            do i=1,p1
                rr=gridpoint(this%grid, i, j, 1)
                write(DJPFD, '(4f19.8)') rr, buf(i,j)
                if (abs(buf(i,j)) > amax) amax=abs(buf(i,j))
            end do
            write(DJPFD, *) 
        end do
        close(DJPFD)
        write(str_g, '(a,e19.12)') 'Max divergence:', amax
        call msg_info(str_g)
        call divj_gopenmol(this, gopen_file)
    end subroutine

    subroutine divj_direct_plt(this)
        type(divj_t), intent(inout) :: this

        integer(I4) :: i, j, k, p1, p2, p3
        real(DP) :: div, amax
        real(DP), dimension(3) :: rr 
        
        if (mpirun_p) then
            call msg_error('divj_direct_plt(): does not work in parallel')
            stop
        end if
        call eta(this%jt, this%grid, 12.d0)
        call get_grid_size(this%grid, p1, p2, p3)
    
        open(DJPFD, file='DIVJPLT')
        amax=D0
        do k=1,p3
            do j=1,p2
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    call divergence(this, rr, div)
!                    call divergence2(this, rr, div)
                    write(DJPFD, '(4f19.8)') rr, div
                    if (abs(div) > amax) amax=abs(div)
                end do
                write(DJPFD, *) 
            end do
        end do
        write(str_g, '(a,e19.12)')'Max divergence:', amax
        call msg_note(str_g)
        call nl
        close(DJPFD)
    end subroutine

    subroutine divj_direct(this, k)
        type(divj_t) :: this
        integer(I4), intent(in) :: k

        integer(I4) :: i, j, p1, p2
        integer(I4) :: lo, hi
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf

        call get_grid_size(this%grid, p1, p2)
        call schedule(p2, lo, hi)

        buf=>this%buf

        do j=1,p2
            do i=1,p1
                rr=gridpoint(this%grid, i, j, k)
                call divergence(this,rr,buf(i,j))
!                call divergence2(this,rr,buf(i,j))
            end do
        end do
        call gather_data(buf,buf(:,lo:hi))
    end subroutine

    subroutine divj(this)
        type(divj_t) :: this

        integer(I4) :: i, j, k, p, p1, p2, p3
        integer(I4) :: lo, hi
        real(DP), dimension(3) :: rr
        real(DP), dimension(:,:), pointer :: buf

        call get_grid_size(this%grid, p1, p2, p3)

        call schedule(p2, lo, hi)

        buf=>this%buf
    
        do k=1,p3
            do j=lo, hi
                do i=1,p1
                    rr=gridpoint(this%grid, i, j, k)
                    call divergence(this,rr,buf(i,j))
!                    call divergence2(this,rr,buf(i,j))
                end do
            end do
            call gather_data(buf,buf(:,lo:hi))
            if (master_p) write(DIVJFD, rec=k) this%buf
        end do
    end subroutine

    subroutine divergence2(this, rr, div)
        type(divj_t), intent(in) :: this
        real(DP), dimension(3), intent(in) :: rr
        real(DP), intent(out) :: div

        integer(I4) :: q
        type(tensor_t) :: jtxp, jtyp, jtzp
        type(tensor_t) :: jtxd, jtyd, jtzd
        real(DP) :: djx,djy,djz
        real(DP), dimension(5) :: jx,jy,jz
        real(DP), dimension(3) :: tvec

        jx=D0; jy=D0; jz=D0
        do q=-2,2
            call ctensor2(this%jt, rr+(/step*real(q),D0,D0/), jtxp,jtxd, 'total')
            call ctensor2(this%jt, rr+(/D0,step*real(q),D0/), jtyp,jtyd, 'total')
            call ctensor2(this%jt, rr+(/D0,D0,step*real(q)/), jtzp,jtzd, 'total')
            call jvector(jtxp%t, jtxd%t, this%bb, tvec)
            jx(q+3)=tvec(1)
            call jvector(jtyp%t, jtyd%t, this%bb, tvec)
            jy(q+3)=tvec(2)
            call jvector(jtzp%t, jtzd%t, this%bb, tvec)
            jz(q+3)=tvec(3)
        end do
        djx=hx*sum(wgt*jx)
        djy=hy*sum(wgt*jy)
        djz=hz*sum(wgt*jz)
        div=djx+djy+djz
    end subroutine

    subroutine divergence(this, rr, div)
        type(divj_t), intent(in) :: this
        real(DP), dimension(3), intent(in) :: rr
        real(DP), intent(out) :: div

        integer(I4) :: q
        type(tensor_t) :: jtx, jty, jtz
        real(DP) :: djx,djy,djz
        real(DP), dimension(5) :: jx,jy,jz
        real(DP), dimension(3) :: tvec

        jx=D0; jy=D0; jz=D0
        do q=-2,2
            if ( q == 0 ) cycle
            call ctensor(this%jt, rr+(/step*real(q),D0,D0/), jtx, 'total')
            call ctensor(this%jt, rr+(/D0,step*real(q),D0/), jty, 'total')
            call ctensor(this%jt, rr+(/D0,D0,step*real(q)/), jtz, 'total')
            tvec=matmul(jtx%t, this%bb)
            jx(q+3)=tvec(1)
            tvec=matmul(jty%t, this%bb)
            jy(q+3)=tvec(2)
            tvec=matmul(jtz%t, this%bb)
            jz(q+3)=tvec(3)
        end do
        djx=hx*sum(wgt*jx)
        djy=hy*sum(wgt*jy)
        djz=hz*sum(wgt*jz)
        div=djx+djy+djz
    end subroutine

    subroutine divj_gopenmol(this, gopen_file)
        type(divj_t) :: this
        character(*), intent(in) :: gopen_file

        integer(I4) :: surface, rank, p1, p2, p3
        integer(I4) :: i, j, k, l
        real(SP), dimension(3) :: qmin, qmax
        real(DP), dimension(:,:), pointer :: buf

        buf=>this%buf
        if (trim(gopen_file) == '') return
        open(GOPFD,file=trim(gopen_file),access='direct',recl=4)

        surface=200
        rank=3

        call get_grid_size(this%grid, p1, p2, p3)
        qmin=real(gridpoint(this%grid,1,1,1)*AU2A)
        qmax=real(gridpoint(this%grid,p1,p2,p3)*AU2A)

        write(GOPFD,rec=1) rank
        write(GOPFD,rec=2) surface
        write(GOPFD,rec=3) p3
        write(GOPFD,rec=4) p2
        write(GOPFD,rec=5) p1
        write(GOPFD,rec=6) qmin(3)
        write(GOPFD,rec=7) qmax(3)
        write(GOPFD,rec=8) qmin(2)
        write(GOPFD,rec=9) qmax(2)
        write(GOPFD,rec=10) qmin(1)
        write(GOPFD,rec=11) qmax(1)

!        write(100,*) rank
!        write(100,*) surface
!        write(100,*) p3
!        write(100,*) p2
!        write(100,*) p1
!        write(100,*) qmin(3)
!        write(100,*) qmax(3)
!        write(100,*) qmin(2)
!        write(100,*) qmax(2)
!        write(100,*) qmin(1)
!        write(100,*) qmax(1)

        l=12
        do k=1,p3
            read(DIVJFD, rec=k) buf
            do j=1,p2
                do i=1,p1
                    write(GOPFD,rec=l) real(buf(i,j))
!                    write(100,*) real(buf(i,j))
                    l=l+1
                end do
            end do
        end do

        close(GOPFD)
    end subroutine
end module

! vim:et:sw=4:ts=4