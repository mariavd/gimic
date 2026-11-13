!
! Divergence of the current density, div J, by central finite differences.
!
! For every grid point r the current density is evaluated at the six points
! r +- h*e_x, r +- h*e_y, r +- h*e_z, and
!
!   div J(r) = sum_a [J_a(r + h*e_a) - J_a(r - h*e_a)] / (2h)
!
! The stencil is built around each point separately, so the result does not
! depend on the kind of grid: evenly spaced, quadrature or read from a file.
! The step h is Advanced.divj_step (in bohr).
!
! Parallelised with OpenMP only.
!
module divj_class
    use globals_module
    use settings_module
    use jtensor_class
    use jfield_class
    use dens_class
    use grid_class
    use basis_class
    use teletype_module
    implicit none

    type divj_t
        real(DP), dimension(3) :: b
        real(DP) :: step
        ! in the order of the grid points
        real(DP), dimension(:), allocatable :: div
        type(grid_t), pointer :: grid
    end type

    public divj_t, new_divj, del_divj, calc_divj, divj_plots
    private

contains

    subroutine new_divj(this, g, magnet)
        type(divj_t) :: this
        type(grid_t), target :: g
        real(DP), dimension(3), intent(in) :: magnet

        integer(I4) :: i, j, k

        this%b=magnet
        this%grid=>g
        this%step=settings%divj_step

        if (this%step <= 0.d0) then
            call msg_error('new_divj(): Advanced.divj_step must be positive!')
            stop 1
        end if

        call get_grid_size(g, i, j, k)
        allocate(this%div(i*j*k))
        this%div=0.d0

        write(str_g, '(a,es10.3,a)') 'Central differences, step = ', &
            this%step, ' bohr'
        call msg_info(str_g)
        call nl
    end subroutine

    subroutine del_divj(this)
        type(divj_t), intent(inout) :: this

        deallocate(this%div)
        nullify(this%grid)
    end subroutine

    subroutine calc_divj(this, mol, xdens, spin)
        type(divj_t) :: this
        type(molecule_t) :: mol
        type(dens_t) :: xdens
        character(*), optional :: spin

        integer(I4) :: i, j, k, p1, p2, p3
        integer(I4) :: n, nn, ip, a, q
        integer(I4) :: npts
        real(DP) :: h, jp, jm
        real(DP), dimension(3) :: r
        real(DP), dimension(:,:), allocatable :: coords, tens
        character(8) :: spincase
        type(jtensor_t) :: jt
        type(jfield_t) :: jf

        if (present(spin)) then
            spincase=spin
        else
            spincase='total'
        end if

        h=this%step
        call get_grid_size(this%grid, p1, p2, p3)
        npts=p1*p2*p3


        ! six evaluations of J per grid point
        call new_jfield(jf, this%grid, this%b)
        call jfield_eta(jf, mol, xdens, 6.d0)
        call del_jfield(jf)

        !$omp parallel default(none) &
        !$omp private(jt,coords,tens,n,nn,ip,i,j,k,a,q,r,jp,jm) &
        !$omp shared(this,mol,xdens,spincase,npts,h)
        call new_jtensor(jt, mol, xdens)
        allocate(coords(3,6*JT_BLOCK), tens(9,6*JT_BLOCK))

        !$omp do schedule(static)
        do n=1,npts,JT_BLOCK
            nn=min(JT_BLOCK, npts-n+1)
            ! point ip, direction a: +h at column 6*(ip-1)+2*a-1, -h at the
            ! next one
            do ip=1,nn
                call get_grid_index(this%grid, n+ip-1, i, j, k)
                r=gridpoint(this%grid, i, j, k)
                do a=1,3
                    q=6*(ip-1)+2*a-1
                    coords(:,q)=r
                    coords(a,q)=r(a)+h
                    coords(:,q+1)=r
                    coords(a,q+1)=r(a)-h
                end do
            end do
            call ctensor_batch(jt, coords(:,1:6*nn), tens(:,1:6*nn), spincase)
            do ip=1,nn
                this%div(n+ip-1)=0.d0
                do a=1,3
                    q=6*(ip-1)+2*a-1
                    ! J_a = sum_b T(a,b) B_b, T(a,b) = tens(a+3*(b-1))
                    jp=dot_product(tens(a:9:3,q), this%b)
                    jm=dot_product(tens(a:9:3,q+1), this%b)
                    this%div(n+ip-1)=this%div(n+ip-1)+(jp-jm)/(2.d0*h)
                end do
            end do
        end do
        !$omp end do

        deallocate(coords, tens)
        call del_jtensor(jt)
        !$omp end parallel
    end subroutine

    ! Write div J to divj.txt (all grids), divj.vti (evenly spaced grids) and
    ! divj.vtu (file and quadrature grids, if the mesh grid.1.ele is there).
    ! Coordinates and div J are in atomic units.
    subroutine divj_plots(this, tag)
        use vtkplot_module
        type(divj_t), intent(inout) :: this
        character(*), optional :: tag

        integer(I4) :: i, j, k, p1, p2, p3, n, np, fd
        integer(int32) :: ncells, dummy1, dummy2
        integer(int32), allocatable :: cells(:,:)
        real(DP), allocatable :: val(:,:,:)
        real(DP), allocatable :: xyz(:,:)
        real(DP) :: dmax, drms
        real(DP), dimension(3) :: rmax
        character(BUFLEN) :: base
        logical :: is_regular, elements_exists

        if (present(tag)) then
            base='divj' // trim(tag)
        else
            base='divj'
        end if

        call get_grid_size(this%grid, p1, p2, p3)
        np=p1*p2*p3

        call getfd(fd)
        open(fd, file=trim(base) // '.txt', status='unknown')
        dmax=-1.d0
        drms=0.d0
        rmax=0.d0
        do n=1,np
            call get_grid_index(this%grid, n, i, j, k)
            write(fd, '(3f16.10,es20.10)') gridpoint(this%grid, i, j, k), &
                this%div(n)
            drms=drms+this%div(n)**2
            if (abs(this%div(n)) > dmax) then
                dmax=abs(this%div(n))
                rmax=gridpoint(this%grid, i, j, k)
            end if
        end do
        call closefd(fd)
        if (np > 0) drms=sqrt(drms/real(np, DP))

        if (present(tag)) then
            call msg_out('Divergence of the current density (' // &
                trim(tag) // '), atomic units:')
        else
            call msg_out('Divergence of the current density, atomic units:')
        end if
        write(str_g, '(a,es14.6)') '   max |div J| = ', dmax
        call msg_out(str_g)
        write(str_g, '(a,3f12.6)') '      at <x,y,z> = ', rmax
        call msg_out(str_g)
        write(str_g, '(a,es14.6)') '   rms  div J  = ', drms
        call msg_out(str_g)
        call msg_out('   written to ' // trim(base) // '.txt')
        call nl

        is_regular = trim(this%grid%mode) /= 'file' .and. &
            trim(this%grid%gtype) == 'even'

        if (is_regular) then
            allocate(val(p1,p2,p3))
            val=reshape(this%div(1:np), (/p1,p2,p3/))
            call write_vtk_imagedata(trim(base) // '.vti', this%grid, val)
            deallocate(val)
        else
            ! tetrahedral mesh: the number of cells is first on the first line
            inquire(file='grid.1.ele', exist=elements_exists)
            if (elements_exists) then
                open(GRIDELE, file='grid.1.ele')
                read(GRIDELE, *) ncells, dummy1, dummy2
                allocate(cells(4,ncells))
                do n=1,ncells
                    read(GRIDELE, *) dummy1, cells(1:4,n)
                end do
                close(GRIDELE)
                allocate(xyz(3,np))
                do n=1,np
                    call get_grid_index(this%grid, n, i, j, k)
                    xyz(:,n)=gridpoint(this%grid, i, j, k)
                end do
                call write_vtk_scalar_unstructuredgrid(trim(base) // '.vtu', &
                    xyz, this%div(1:np), cells)
                deallocate(xyz, cells)
            end if
        end if
    end subroutine

end module

! vim:et:sw=4:ts=4
