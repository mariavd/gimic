!
! Current tensor evaluation with OpenMP target offload (GPU).
!
! The basis set and the density matrices are flattened into plain arrays
! and kept on the device for the whole run. Grid points are processed in
! super-blocks of up to npb points:
!
!   1. one kernel evaluates basis values, gradients and (with GIAOs) the
!      magnetic derivatives for every point of the block,
!   2. each density matrix is applied with one matrix-matrix product,
!   3. one kernel per magnetic-field component forms the dot products and
!      assembles the tensor.
!
! The arithmetic mirrors jtensor_class (cgto_dr, d2fdrdb, contract_batch)
! term by term, so results agree with the CPU path up to summation order
! in the matrix products.
!
! Without a device the target regions run on the host, which is what the
! tests use (GIMIC_OFFLOAD=host). Selection at run time:
!   GIMIC_OFFLOAD=auto (default)  use the GPU if one is present
!   GIMIC_OFFLOAD=host            run this code path on the CPU
!   GIMIC_OFFLOAD=off             use the CPU path in jtensor_class
!
module jtensor_offload
    use globals_module
    use settings_module
    use basis_class
    use dens_class
    use gtodefs_module
    use teletype_module
#ifdef HAVE_OPENMP
    use omp_lib
#endif
    implicit none

    public offload_init, offload_finalize, offload_active, ctensor_offload
    private

    integer(I4), parameter :: MAXCOMP=21      ! Cartesian components for l=MAX_L
    integer(I4), parameter :: OFFLOAD_NPB=4096 ! max points per super-block

    logical, save :: active=.false.
    logical, save :: initialized=.false.
    logical, save :: giao=.true.
    integer(I4), save :: nv=0, natoms=0, nctr=0, npf=0, nspin=1, npb=0

    ! flattened basis
    real(DP), allocatable, save :: at_coord(:,:)                 ! (3,natoms)
    integer(I4), allocatable, save :: ctr_atom(:), ctr_l(:), ctr_ncomp(:)
    integer(I4), allocatable, save :: ctr_pos(:), ctr_pf0(:), ctr_npf(:)
    real(DP), allocatable, save :: ctr_thrs(:)                   ! screening radius
    real(DP), allocatable, save :: pf_xp(:), pf_ncc(:)           ! primitives
    integer(I4), save :: nlm(3,MAXCOMP,0:MAX_L)                  ! Cartesian powers

    ! densities: dmat0 = D, pt(:,:,b,:) = transpose(P_b), per spin
    real(DP), allocatable, save :: dmat0(:,:,:), pt(:,:,:,:)

    ! per-super-block work arrays
    real(DP), allocatable, save :: rr(:,:)                       ! (3,npb) points
    real(DP), allocatable, save :: bf(:,:)                       ! (nv,npb)
    real(DP), allocatable, save :: dr(:,:,:)                     ! (nv,3,npb)
    real(DP), allocatable, save :: db(:,:,:)                     ! (nv,npb,3)
    real(DP), allocatable, save :: d2m(:,:,:)                     ! (nv,9,npb)
    real(DP), allocatable, save :: denbf(:,:), pdbf(:,:), dendb(:,:)
    real(DP), allocatable, save :: diapam(:)
    real(DP), allocatable, save :: ct(:,:)                       ! (9,npb)

contains

    function offload_active() result(r)
        logical :: r
        r=active
    end function

    subroutine offload_init(mol, xdens)
        type(molecule_t), target :: mol
        type(dens_t), target :: xdens

        character(16) :: mode
        integer(I4) :: ndev, i, j, k, c, l, s, b, idx2, status
        type(atom_t), pointer :: atom
        type(basis_t), pointer :: basis
        type(contraction_t), pointer :: ctr
        real(DP), dimension(:,:), pointer :: dmat, gt
        real(DP) :: mem

        active=.false.
        mode='auto'
        call get_environment_variable('GIMIC_OFFLOAD', mode, status=status)
        if (status /= 0) mode='auto'

        ndev=0
#ifdef HAVE_OPENMP
        ndev=omp_get_num_devices()
#endif
        select case (trim(mode))
            case ('off', 'OFF', '0', 'no', 'NO', 'false', 'FALSE')
                return
            case ('host', 'HOST')
                active=.true.
                call msg_info('Offload: running the offload code path on the host')
            case default
                active = (ndev > 0)
                if (active) then
                    write(str_g, '(a,i2,a)') 'Offload: ', ndev, ' device(s) found'
                    call msg_info(str_g)
                end if
        end select
        if (.not. active) return

        if (settings%use_spherical) then
            call msg_info('Offload: spherical basis not supported, using the CPU path')
            active=.false.
            return
        end if

        giao=settings%use_giao
        nv=get_ncgto(mol)
        natoms=get_natoms(mol)

        ! count contractions and primitives
        nctr=0
        npf=0
        do i=1,natoms
            call get_atom(mol, i, atom)
            call get_basis(atom, basis)
            do j=1,get_nctr(basis)
                call get_contraction(atom, j, ctr)
                nctr=nctr+1
                npf=npf+ctr%npf
            end do
        end do

        allocate(at_coord(3,natoms))
        allocate(ctr_atom(nctr), ctr_l(nctr), ctr_ncomp(nctr), ctr_pos(nctr))
        allocate(ctr_pf0(nctr), ctr_npf(nctr), ctr_thrs(nctr))
        allocate(pf_xp(npf), pf_ncc(npf))

        ! same ordering and offsets as bfeval_dr
        c=0
        k=0
        idx2=0
        do i=1,natoms
            call get_atom(mol, i, atom)
            call get_coord(atom, at_coord(:,i))
            call get_basis(atom, basis)
            do j=1,get_nctr(basis)
                call get_contraction(atom, j, ctr)
                c=c+1
                ctr_atom(c)=i
                ctr_l(c)=ctr%l
                ctr_ncomp(c)=ctr%nccomp
                ctr_pos(c)=idx2+get_ctridx(basis, j)
                ctr_thrs(c)=basis%thrs(j)
                ctr_pf0(c)=k+1
                ctr_npf(c)=ctr%npf
                pf_xp(k+1:k+ctr%npf)=ctr%xp(1:ctr%npf)
                pf_ncc(k+1:k+ctr%npf)=ctr%ncc(1:ctr%npf)
                k=k+ctr%npf
            end do
            idx2=idx2+get_ncgto(basis)
        end do

        nlm=0
        do l=0,MAX_L
            call get_gto_nlm(l, gt)
            do i=1,size(gt,2)
                nlm(:,i,l)=nint(gt(:,i))
            end do
        end do

        ! density matrices (spin index 1 = alpha, 2 = beta)
        nspin=1
        if (settings%is_uhf) nspin=2
        allocate(dmat0(nv,nv,nspin), pt(nv,nv,3,nspin))
        do s=1,nspin
            call get_dens(xdens, dmat, s)
            dmat0(:,:,s)=dmat(1:nv,1:nv)
            do b=1,3
                call get_pdens(xdens, b, dmat, s)
                pt(:,:,b,s)=transpose(dmat(1:nv,1:nv))
            end do
        end do

        ! block size: about 1.5 GB of work arrays at most
        npb=int(1.5d9/(real(nv,DP)*19.0d0*8.0d0))
        npb=max(64, min(OFFLOAD_NPB, npb))
        allocate(rr(3,npb), bf(nv,npb), dr(nv,3,npb), db(nv,npb,3), d2m(nv,9,npb))
        allocate(denbf(nv,npb), pdbf(nv,npb), dendb(nv,npb), diapam(npb), ct(9,npb))

        mem=real(nv,DP)*real(nv,DP)*8.0d0*4.0d0*nspin + real(nv,DP)*19.0d0*8.0d0*npb
        write(str_g, '(a,i6,a,f8.1,a)') 'Offload: block size ', npb, &
            ' points, device memory ', mem/1.0d6, ' MB'
        call msg_info(str_g)
        call nl

        !$omp target enter data map(to: at_coord, ctr_atom, ctr_l, ctr_ncomp, ctr_pos, &
        !$omp&  ctr_pf0, ctr_npf, ctr_thrs, pf_xp, pf_ncc, nlm, dmat0, pt) &
        !$omp&  map(alloc: rr, bf, dr, db, d2m, denbf, pdbf, dendb, diapam, ct)
        initialized=.true.
    end subroutine

    subroutine offload_finalize()
        if (.not. initialized) return
        !$omp target exit data map(delete: at_coord, ctr_atom, ctr_l, ctr_ncomp, ctr_pos, &
        !$omp&  ctr_pf0, ctr_npf, ctr_thrs, pf_xp, pf_ncc, nlm, dmat0, pt, &
        !$omp&  rr, bf, dr, db, d2m, denbf, pdbf, dendb, diapam, ct)
        deallocate(at_coord, ctr_atom, ctr_l, ctr_ncomp, ctr_pos, ctr_pf0, ctr_npf, ctr_thrs)
        deallocate(pf_xp, pf_ncc, dmat0, pt)
        deallocate(rr, bf, dr, db, d2m, denbf, pdbf, dendb, diapam, ct)
        initialized=.false.
        active=.false.
    end subroutine

!
! Same as ctensor_batch: r(3,npts) -> j(9,npts) for the spin case op.
!
    subroutine ctensor_offload(r, j, op)
        real(DP), dimension(:,:), intent(in) :: r
        real(DP), dimension(:,:), intent(out) :: j
        character(*), intent(in) :: op

        integer(I4) :: i0, i1, n, npts

        if (.not. initialized) then
            call msg_error('ctensor_offload(): not initialized')
            stop
        end if
        if ((op == 'beta' .or. op == 'spindens') .and. .not.settings%is_uhf) then
            call msg_error('ctensor_offload(): ' // trim(op) // &
                ' requested, but not open-shell system!')
            stop
        end if

        npts=size(r,2)
        do i0=1,npts,npb
            n=min(npb, npts-i0+1)
            i1=i0+n-1
            rr(:,1:n)=r(:,i0:i1)
            !$omp target update to(rr(:,1:n))
            call eval_block(n)
            select case (op)
                case ('alpha')
                    call contract_block(n, 1)
                    !$omp target update from(ct(:,1:n))
                    j(:,i0:i1)=ct(:,1:n)
                case ('beta')
                    call contract_block(n, 2)
                    !$omp target update from(ct(:,1:n))
                    j(:,i0:i1)=ct(:,1:n)
                case ('total')
                    call contract_block(n, 1)
                    !$omp target update from(ct(:,1:n))
                    j(:,i0:i1)=ct(:,1:n)
                    if (settings%is_uhf) then
                        call contract_block(n, 2)
                        !$omp target update from(ct(:,1:n))
                        j(:,i0:i1)=j(:,i0:i1)+ct(:,1:n)
                    end if
                case ('spindens')
                    call contract_block(n, 1)
                    !$omp target update from(ct(:,1:n))
                    j(:,i0:i1)=ct(:,1:n)
                    call contract_block(n, 2)
                    !$omp target update from(ct(:,1:n))
                    j(:,i0:i1)=j(:,i0:i1)-ct(:,1:n)
                case default
                    call msg_error('ctensor_offload(): unknown spin case ' // trim(op))
                    stop
            end select
        end do
    end subroutine

!
! Basis values, gradients and magnetic derivatives for np points
! (cgto_dr + mkdbop + dfdb + d2fdrdb, per point).
!
    subroutine eval_block(np)
        integer(I4), intent(in) :: np

        integer(I4) :: ip, c, i, jp, l, nc, pos, ia, idx
        integer(I4) :: n1, n2, n3
        real(DP) :: x, y, z, dx, dy, dz, r2, q, bfval, dbfval, p, v
        real(DP) :: g1, g2, g3, rx, ry, rz, b1, b2, b3
        real(DP) :: rp(0:MAX_L,3)
        logical :: use_giao

        use_giao=giao

        !$omp target teams distribute parallel do &
        !$omp&  private(c,i,jp,l,nc,pos,ia,idx,n1,n2,n3,x,y,z,dx,dy,dz,r2,q, &
        !$omp&          bfval,dbfval,p,v,g1,g2,g3,rx,ry,rz,b1,b2,b3,rp)
        do ip=1,np
            x=rr(1,ip)
            y=rr(2,ip)
            z=rr(3,ip)
            do c=1,nctr
                ia=ctr_atom(c)
                rx=at_coord(1,ia)
                ry=at_coord(2,ia)
                rz=at_coord(3,ia)
                dx=x-rx
                dy=y-ry
                dz=z-rz
                r2=dx*dx+dy*dy+dz*dz
                pos=ctr_pos(c)
                nc=ctr_ncomp(c)
                l=ctr_l(c)

                ! screening, as filter_screened
                if (sqrt(r2) > ctr_thrs(c)) then
                    do i=0,nc-1
                        bf(pos+i,ip)=D0
                        dr(pos+i,1,ip)=D0
                        dr(pos+i,2,ip)=D0
                        dr(pos+i,3,ip)=D0
                        if (use_giao) then
                            db(pos+i,ip,1)=D0
                            db(pos+i,ip,2)=D0
                            db(pos+i,ip,3)=D0
                            d2m(pos+i,1:9,ip)=D0
                        end if
                    end do
                    cycle
                end if

                ! radial part (cao2)
                bfval=D0
                dbfval=D0
                do jp=ctr_pf0(c),ctr_pf0(c)+ctr_npf(c)-1
                    q=pf_ncc(jp)*exp(-pf_xp(jp)*r2)
                    bfval=bfval+q
                    dbfval=dbfval+pf_xp(jp)*q
                end do

                rp(0,1)=D1
                rp(0,2)=D1
                rp(0,3)=D1
                do i=1,l
                    rp(i,1)=rp(i-1,1)*dx
                    rp(i,2)=rp(i-1,2)*dy
                    rp(i,3)=rp(i-1,3)*dz
                end do

                if (use_giao) then
                    ! mkdbop: r x R for this atom
                    b1=y*rz-z*ry
                    b2=z*rx-x*rz
                    b3=x*ry-y*rx
                end if

                do i=1,nc
                    idx=pos+i-1
                    n1=nlm(1,i,l)
                    n2=nlm(2,i,l)
                    n3=nlm(3,i,l)
                    p=rp(n1,1)*rp(n2,2)*rp(n3,3)
                    v=p*bfval
                    bf(idx,ip)=v
                    ! gradient, as cgto_dr
                    if (n1 > 0) then
                        q=rp(n1-1,1)*rp(n2,2)*rp(n3,3)
                        g1=n1*q*bfval-2.d0*dx*p*dbfval
                    else
                        g1=-2.d0*dx*p*dbfval
                    end if
                    if (n2 > 0) then
                        q=rp(n1,1)*rp(n2-1,2)*rp(n3,3)
                        g2=n2*q*bfval-2.d0*dy*p*dbfval
                    else
                        g2=-2.d0*dy*p*dbfval
                    end if
                    if (n3 > 0) then
                        q=rp(n1,1)*rp(n2,2)*rp(n3-1,3)
                        g3=n3*q*bfval-2.d0*dz*p*dbfval
                    else
                        g3=-2.d0*dz*p*dbfval
                    end if
                    dr(idx,1,ip)=g1
                    dr(idx,2,ip)=g2
                    dr(idx,3,ip)=g3
                    if (use_giao) then
                        ! dfdb
                        db(idx,ip,1)=b1*v
                        db(idx,ip,2)=b2*v
                        db(idx,ip,3)=b3*v
                        ! d2fdrdb
                        d2m(idx,1,ip)=g1*b1
                        d2m(idx,2,ip)=g2*b1+rz*v
                        d2m(idx,3,ip)=g3*b1-ry*v
                        d2m(idx,4,ip)=g1*b2-rz*v
                        d2m(idx,5,ip)=g2*b2
                        d2m(idx,6,ip)=g3*b2+rx*v
                        d2m(idx,7,ip)=g1*b3+ry*v
                        d2m(idx,8,ip)=g2*b3-rx*v
                        d2m(idx,9,ip)=g3*b3
                    end if
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine

!
! y(1:nv,1:np) = a(1:nv,1:nv) x(1:nv,1:np) on the device. Plain kernel:
! consecutive threads take consecutive rows i, so a(i,k) is read coalesced
! and x(k,ip) is shared within a team. A vendor BLAS (cublasDgemm through
! use_device_addr) can replace this later.
!
    subroutine dev_gemm(a, x, y, np)
        real(DP), dimension(:,:), intent(in) :: a, x
        real(DP), dimension(:,:), intent(inout) :: y
        integer(I4), intent(in) :: np

        integer(I4) :: i, k, ip, n
        real(DP) :: s

        n=nv
        !$omp target teams distribute parallel do collapse(2) private(k,s)
        do ip=1,np
            do i=1,n
                s=D0
                do k=1,n
                    s=s+a(i,k)*x(k,ip)
                end do
                y(i,ip)=s
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine

!
! contract_batch for the block: density products and per-point dots.
!
    subroutine contract_block(np, is)
        integer(I4), intent(in) :: np, is

        integer(I4) :: ip, i, b, m, k, n
        real(DP) :: s, ppd, p1, p2, dp1, dp2, dp3
        logical :: use_giao, use_para, use_dia

        n=nv
        use_giao=giao
        use_para=settings%use_paramag
        use_dia=settings%use_diamag
        if (.not.use_para) your_results_are_questionable=.true.
        if (.not.use_dia) your_results_are_questionable=.true.

        call dev_gemm(dmat0(:,:,is), bf, denbf, np)

        !$omp target teams distribute parallel do private(i,s)
        do ip=1,np
            s=D0
            do i=1,n
                s=s+denbf(i,ip)*bf(i,ip)
            end do
            diapam(ip)=s
        end do
        !$omp end target teams distribute parallel do

        do b=1,3
            call dev_gemm(pt(:,:,b,is), bf, pdbf, np)
            if (use_giao) call dev_gemm(dmat0(:,:,is), db(:,:,b), dendb, np)

            !$omp target teams distribute parallel do private(i,m,k,ppd,p1,p2)
            do ip=1,np
                do m=1,3
                    k=(b-1)*3+m
                    ppd=D0
                    p1=D0
                    p2=D0
                    do i=1,n
                        ppd=ppd+pdbf(i,ip)*dr(i,m,ip)
                    end do
                    if (use_giao) then
                        do i=1,n
                            p1=p1+dendb(i,ip)*dr(i,m,ip)
                            p2=p2+denbf(i,ip)*d2m(i,k,ip)
                        end do
                        ct(k,ip)=ZETA*ppd+ZETA*(-p1+p2)
                    else
                        ct(k,ip)=ZETA*ppd
                    end if
                end do
            end do
            !$omp end target teams distribute parallel do
        end do

        ! diamagnetic contribution, off-diagonal only (see contract_batch)
        !$omp target teams distribute parallel do private(dp1,dp2,dp3)
        do ip=1,np
            dp1=diapam(ip)*DP50*rr(1,ip)
            dp2=diapam(ip)*DP50*rr(2,ip)
            dp3=diapam(ip)*DP50*rr(3,ip)
            if (.not.use_para) ct(1:9,ip)=D0
            if (.not.use_dia) then
                dp1=D0
                dp2=D0
                dp3=D0
            end if
            ct(4,ip)=ct(4,ip)+dp3
            ct(7,ip)=ct(7,ip)-dp2
            ct(2,ip)=ct(2,ip)-dp3
            ct(8,ip)=ct(8,ip)+dp1
            ct(3,ip)=ct(3,ip)+dp2
            ct(6,ip)=ct(6,ip)-dp1
        end do
        !$omp end target teams distribute parallel do
    end subroutine

end module
