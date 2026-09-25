!
! This is the actual work horse, calculates the current tensor
! for a particular spin case. The cost of calling jtensor twice is
! very small, since no exponentials have to be calculated.
!
module jtensor_class
    use globals_module
    use settings_module
    use basis_class
    use bfeval_class
    use dens_class
    use teletype_module
    implicit none

!    intrinsic dtime

    type jtensor_t
        private
        type(molecule_t), pointer :: mol
        type(dens_t), pointer :: xdens
        type(bfeval_t) :: basv
        real(DP), dimension(:), pointer :: pdbf, denbf, dendb
        real(DP), dimension(3) :: rho
        ! scratch mem
        real(DP), dimension(:), pointer :: bfvec
        real(DP), dimension(:,:), pointer :: dbvec, drvec, d2fvec, dbop
        real(DP), dimension(:,:), pointer :: aodens, pdens
        ! block scratch for ctensor_batch: basis functions and their
        ! derivatives for up to JT_BLOCK grid points side by side, so that
        ! the density matrices are applied with one matrix-matrix product
        ! per block instead of one matrix-vector product per point
        integer(I4) :: nv
        real(DP), dimension(:,:), allocatable :: bfm, denbfm, pdbfm, dendbm
        real(DP), dimension(:,:,:), allocatable :: drm, dbm, d2m
        ! Screening: only the basis functions whose value or gradient is
        ! non-zero somewhere in the block ("active") take part in the
        ! products. act(1:na) lists them; the slabs are compacted to those
        ! rows and dsub holds the corresponding sub-matrix of a density.
        integer(I4) :: na
        integer(I4), dimension(:), allocatable :: act
        real(DP), dimension(:,:), allocatable :: dsub, psub
    end type

    public new_jtensor, del_jtensor, jtensor, get_jvector
    public ctensor, ctensor_batch, jvector
    public jtensor_t, JT_BLOCK

    private
    integer(I4), parameter :: NOTIFICATION=1000
    ! number of grid points contracted per matrix-matrix product
    integer(I4), parameter :: JT_BLOCK=64

#ifdef HAVE_BLAS
    ! Explicit interfaces (assumed-size dummies, as in reference BLAS).
    interface
        function ddot(n, x, incx, y, incy)
            import :: DP
            integer, intent(in) :: n, incx, incy
            real(DP), intent(in) :: x(*), y(*)
            real(DP) :: ddot
        end function
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: DP
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(DP), intent(in) :: alpha, beta
            real(DP), intent(in) :: a(lda, *), x(*)
            real(DP), intent(inout) :: y(*)
        end subroutine
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: DP
            character, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(DP), intent(in) :: alpha, beta
            real(DP), intent(in) :: a(lda, *), b(ldb, *)
            real(DP), intent(inout) :: c(ldc, *)
        end subroutine
    end interface
#endif

contains
    ! set up memory (once) for the different components
    subroutine new_jtensor(this, mol, xdens)
        type(jtensor_t) :: this
        type(molecule_t), target :: mol
        type(dens_t), target :: xdens
        integer(I4) ::  ncgto, nccgto

        ncgto=get_ncgto(mol)
        nccgto=get_nccgto(mol)

        this%mol=>mol
        this%xdens=>xdens
        call new_bfeval(this%basv, this%mol)

        ! intermediates
        allocate(this%denbf(ncgto))
        allocate(this%dendb(ncgto))
        allocate(this%pdbf(ncgto))

        ! block scratch (sizes mirror those in new_bfeval)
        this%nv=0
        this%na=0
        allocate(this%act(ncgto))
        allocate(this%dsub(ncgto,ncgto))
        allocate(this%psub(ncgto,ncgto))
        allocate(this%bfm(ncgto,JT_BLOCK))
        allocate(this%denbfm(ncgto,JT_BLOCK))
        allocate(this%pdbfm(ncgto,JT_BLOCK))
        allocate(this%drm(nccgto,3,JT_BLOCK))
        if (settings%use_giao) then
            allocate(this%dendbm(ncgto,JT_BLOCK))
            allocate(this%dbm(ncgto,JT_BLOCK,3))
            allocate(this%d2m(ncgto,9,JT_BLOCK))
        end if
    end subroutine

    subroutine del_jtensor(this)
        type(jtensor_t) :: this

        call del_bfeval(this%basv)
        deallocate(this%denbf)
        deallocate(this%pdbf)
        deallocate(this%dendb)
        deallocate(this%bfm, this%denbfm, this%pdbfm, this%drm)
        deallocate(this%act, this%dsub, this%psub)
        if (allocated(this%dendbm)) deallocate(this%dendbm, this%dbm, this%d2m)
    end subroutine

!
! Same as ctensor, for a batch of points r(3,npts) -> j(9,npts).
! Points are processed JT_BLOCK at a time; see contract_batch.
!
    subroutine ctensor_batch(this, r, j, op)
        type(jtensor_t) :: this
        real(DP), dimension(:,:), intent(in) :: r
        real(DP), dimension(:,:), intent(out) :: j
        character(*) :: op

        integer(I4) :: i0, i1, n, npts
        real(DP), dimension(9,JT_BLOCK) :: jt2

        npts=size(r,2)
        do i0=1,npts,JT_BLOCK
            n=min(JT_BLOCK, npts-i0+1)
            i1=i0+n-1
            call eval_basis_block(this, n, r(:,i0:i1))
            select case (op)
                case ('alpha')
                    call contract_batch(this, n, r(:,i0:i1), j(:,i0:i1), spin_a)
                case ('beta')
                    if (settings%is_uhf) then
                        call contract_batch(this, n, r(:,i0:i1), j(:,i0:i1), spin_b)
                    else
                        call msg_error('ctensor_batch(): &
                        &beta current requested, but not open-shell system!')
                        stop
                    end if
                case ('total')
                    call contract_batch(this, n, r(:,i0:i1), j(:,i0:i1), spin_a)
                    if (settings%is_uhf) then
                        call contract_batch(this, n, r(:,i0:i1), jt2(:,1:n), spin_b)
                        j(:,i0:i1)=j(:,i0:i1)+jt2(:,1:n)
                    end if
                case ('spindens')
                    if (.not.settings%is_uhf) then
                        call msg_error('ctensor_batch(): &
                        &spindens requested, but not open-shell system!')
                        stop
                    end if
                    call contract_batch(this, n, r(:,i0:i1), j(:,i0:i1), spin_a)
                    call contract_batch(this, n, r(:,i0:i1), jt2(:,1:n), spin_b)
                    j(:,i0:i1)=j(:,i0:i1)-jt2(:,1:n)
            end select
        end do
    end subroutine

    ! evaluate basis functions and derivatives for n points into the
    ! block scratch arrays
    subroutine eval_basis_block(this, n, r)
        type(jtensor_t) :: this
        integer(I4), intent(in) :: n
        real(DP), dimension(:,:), intent(in) :: r

        integer(I4) :: ip, nv

        do ip=1,n
            if (settings%use_giao) then
                call calc_basis(this%basv, r(:,ip), this%bfvec, this%drvec, &
                this%dbvec, this%d2fvec)
            else
                call calc_basis(this%basv, r(:,ip), this%bfvec, this%drvec)
            end if
            nv=size(this%bfvec)
            this%nv=nv
            this%bfm(1:nv,ip)=this%bfvec
            this%drm(1:nv,:,ip)=this%drvec(1:nv,:)
            if (settings%use_giao) then
                this%dbm(1:nv,ip,:)=this%dbvec(1:nv,:)
                this%d2m(1:nv,:,ip)=this%d2fvec(1:nv,:)
            end if
        end do
        call compact_block(this, n)
    end subroutine

    ! find the functions with a non-zero value or gradient anywhere in the
    ! block and move their rows to the front of the slabs (act is increasing,
    ! so this is safe in place). The GIAO slabs are built from value and
    ! gradient, so they vanish on the same rows.
    subroutine compact_block(this, n)
        type(jtensor_t) :: this
        integer(I4), intent(in) :: n

        integer(I4) :: i, k, na, nv

        nv=this%nv
        na=0
        do i=1,nv
            if (any(this%bfm(i,1:n) /= D0) .or. any(this%drm(i,:,1:n) /= D0)) then
                na=na+1
                this%act(na)=i
            end if
        end do
        this%na=na
        if (na == nv) return

        do k=1,na
            i=this%act(k)
            if (i == k) cycle
            this%bfm(k,1:n)=this%bfm(i,1:n)
            this%drm(k,:,1:n)=this%drm(i,:,1:n)
            if (settings%use_giao) then
                this%dbm(k,1:n,:)=this%dbm(i,1:n,:)
                this%d2m(k,:,1:n)=this%d2m(i,:,1:n)
            end if
        end do
    end subroutine

    ! sub(1:na,1:na) = d(act,act)
    subroutine gather_sub(this, d, sub)
        type(jtensor_t) :: this
        real(DP), dimension(:,:), intent(in) :: d
        real(DP), dimension(:,:), intent(out) :: sub

        integer(I4) :: k, l, jl

        do l=1,this%na
            jl=this%act(l)
            do k=1,this%na
                sub(k,l)=d(this%act(k),jl)
            end do
        end do
    end subroutine

    ! y(1:na,1:n) = op(d) x(1:na,1:n) over the active functions, where
    ! d is the density matrix restricted to them (dsub, or the full matrix
    ! when nothing is screened)
    subroutine dens_times(this, trans, d, x, y, n)
        type(jtensor_t) :: this
        character, intent(in) :: trans
        real(DP), dimension(:,:), intent(in) :: d, x
        real(DP), dimension(:,:), intent(inout) :: y
        integer(I4), intent(in) :: n

        integer(I4) :: na

        na=this%na
#ifdef HAVE_BLAS
        call dgemm(trans, 'n', na, n, na, 1.0d0, d, size(d,1), &
            x, size(x,1), 0.0d0, y, size(y,1))
#else
        if (trans == 't') then
            y(1:na,1:n)=matmul(transpose(d(1:na,1:na)), x(1:na,1:n))
        else
            y(1:na,1:n)=matmul(d(1:na,1:na), x(1:na,1:n))
        end if
#endif
    end subroutine

!
! Batched version of contract(): the same contractions for n points at
! once. The three density-matrix products, which dominate the cost, become
! matrix-matrix products so that each density matrix is read once per
! block rather than once per point. Results equal contract() up to
! summation order.
!
    subroutine contract_batch(this, n, r, ct, spin)
        type(jtensor_t) :: this
        integer(I4), intent(in) :: n
        real(DP), dimension(:,:), intent(in) :: r
        real(DP), dimension(:,:), intent(out) :: ct
        integer(I4), intent(in) :: spin

        integer(I4) :: b, m, k, ip, nv
        logical :: screened
        real(DP) :: prsp1, prsp2, ppd
        real(DP), dimension(JT_BLOCK) :: diapam
        real(DP), dimension(3,JT_BLOCK) :: dpd

        nv=this%na
        screened = (this%na < this%nv)
        call get_dens(this%xdens, this%aodens, spin)
        if (screened) then
            call gather_sub(this, this%aodens, this%dsub)
            call dens_times(this, 'n', this%dsub, this%bfm, this%denbfm, n)
        else
            call dens_times(this, 'n', this%aodens, this%bfm, this%denbfm, n)
        end if
        do ip=1,n
            diapam(ip)=dot_product(this%denbfm(1:nv,ip), this%bfm(1:nv,ip))
        end do

        do b=1,3 ! dB <x,y,z>
            if (settings%use_giao) then
                if (screened) then
                    call dens_times(this, 'n', this%dsub, this%dbm(:,:,b), this%dendbm, n)
                else
                    call dens_times(this, 'n', this%aodens, this%dbm(:,:,b), this%dendbm, n)
                end if
            end if
            call get_pdens(this%xdens, b, this%pdens, spin)
            if (screened) then
                call gather_sub(this, this%pdens, this%psub)
                call dens_times(this, 't', this%psub, this%bfm, this%pdbfm, n)
            else
                call dens_times(this, 't', this%pdens, this%bfm, this%pdbfm, n)
            end if
            do ip=1,n
                dpd(b,ip)=diapam(ip)*DP50*r(b,ip) ! diamag. contr. to J
                do m=1,3 ! dm <x,y,z>
                    k=(b-1)*3+m
                    ppd=dot_product(this%pdbfm(1:nv,ip), this%drm(1:nv,m,ip))
                    ct(k,ip)=ZETA*ppd
                    if (settings%use_giao) then
                        prsp1=-dot_product(this%dendbm(1:nv,ip), this%drm(1:nv,m,ip)) ! (-i)**2=-1
                        prsp2=dot_product(this%denbfm(1:nv,ip), this%d2m(1:nv,k,ip))
                        ct(k,ip)=ct(k,ip)+ZETA*(prsp1+prsp2)
                    end if
                end do
            end do
        end do

        ! annihilate paramagnetic contribution
        if (.not.settings%use_paramag) then
            ct(:,1:n)=D0
            your_results_are_questionable = .true.
        end if
        ! annihilate diamagnetic  contribution
        if (.not.settings%use_diamag) then
            dpd=D0
            your_results_are_questionable = .true.
        end if

        ! the diamagnetic probability density only contributes off-diagonal;
        ! ct(m,b) is stored at ct(m+3*(b-1)), as in contract()
        do ip=1,n
            ct(4,ip)=ct(4,ip)+dpd(3,ip) ! (1,2)
            ct(7,ip)=ct(7,ip)-dpd(2,ip) ! (1,3)
            ct(2,ip)=ct(2,ip)-dpd(3,ip) ! (2,1)
            ct(8,ip)=ct(8,ip)+dpd(1,ip) ! (2,3)
            ct(3,ip)=ct(3,ip)+dpd(2,ip) ! (3,1)
            ct(6,ip)=ct(6,ip)-dpd(1,ip) ! (3,2)
        end do
    end subroutine

    subroutine ctensor(this, r, j, op)
        type(jtensor_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(9), intent(inout) :: j
        character(*) :: op

        real(DP), dimension(9) :: jt1, jt2

        select case (op)
            case ('alpha')
                call jtensor(this, r, j, spin_a)
            case ('beta')
                if (settings%is_uhf) then
                    call jtensor(this, r, j, spin_b)
                else
                    call msg_error('ctensor(): &
                    &beta current requested, but not open-shell system!')
                    stop
                end if
            case ('total')
                if (settings%is_uhf) then
                    call jtensor(this, r, jt1, spin_a)
                    call jtensor(this, r, jt2, spin_b)
                    j=jt1+jt2
                else
                    call jtensor(this, r, j, spin_a)
                end if
            case ('spindens')
                if (.not.settings%is_uhf) then
                    call msg_error('ctensor(): &
                    &spindens requested, but not open-shell system!')
                    stop
                end if
                call jtensor(this, r, jt1, spin_a)
                call jtensor(this, r, jt2, spin_b)
                j=jt1-jt2
        end select
    end subroutine

    subroutine jtensor(this, r, j, spin)
        type(jtensor_t) :: this
        real(DP), dimension(:), intent(in) :: r
        real(DP), dimension(9), intent(inout) :: j
        integer(I4) :: spin

        integer(I4) :: i, b
        integer(I4), save :: notify=1

        this%rho=DP50*r ! needed for diamag. contr.
        if (settings%use_giao) then
            call calc_basis(this%basv, r, this%bfvec, this%drvec, &
            this%dbvec, this%d2fvec)
        else
            call calc_basis(this%basv, r, this%bfvec, this%drvec)
        end if

        call contract(this, j, spin)
    end subroutine

    subroutine jvector(this, r, bb, jv, op)
        type(jtensor_t) :: this
        real(DP), dimension(3), intent(in) :: r
        real(DP), dimension(:), intent(in) :: bb
        real(DP), dimension(:), intent(out) :: jv
        character(*) :: op

        real(DP), dimension(9) :: j
        call ctensor(this, r, j, op)
        jv=matmul(reshape(j,(/3,3/)), bb)
    end subroutine

    subroutine get_jvector(pj, dj, bb, jv)
        real(DP), dimension(:), intent(in) :: pj, dj
        real(DP), dimension(:), intent(in) :: bb
        real(DP), dimension(:), intent(out) :: jv

        jv=matmul(reshape(pj+dj,(/3,3/)), bb)
    end subroutine

!
! Contract all contributions to J. This is where all the actual work is done.
!
    subroutine contract(this, ct, spin)
        type(jtensor_t) :: this
        real(DP), dimension(3,3), intent(out) :: ct
        integer(I4), intent(in) :: spin

        integer(I4) :: b, m, n, k, ii,jj
        integer(I4) :: vec_size
        real(DP) :: prsp1, prsp2       ! paramagnetic wavefunction response
        real(DP) :: ppd                ! paramagnetic probability density
        real(DP), dimension(3) :: dpd  ! diamagnetic probability density
        real(DP) :: diapam

        call get_dens(this%xdens, this%aodens, spin)
#ifdef HAVE_BLAS
        vec_size = size(this%bfvec)
        call dgemv('n', vec_size, vec_size, 1.0d0, this%aodens, vec_size, this%bfvec, 1, 0.0d0, this%denbf, 1)
        diapam = ddot(vec_size, this%denbf, 1, this%bfvec, 1)
#else
        this%denbf=matmul(this%bfvec, this%aodens)
        diapam=dot_product(this%denbf, this%bfvec)
#endif

        k=1
        do b=1,3! dB <x,y,z>
            ! get perturbed densities: x,y,z
            call get_pdens(this%xdens, b, this%pdens, spin)
#ifdef HAVE_BLAS
            call dgemv('t', vec_size, vec_size, 1.0d0, this%pdens, vec_size, this%bfvec, 1, 0.0d0, this%pdbf, 1)
#else
            this%pdbf=matmul(this%bfvec, this%pdens)
#endif
            if (settings%use_giao) then
#ifdef HAVE_BLAS
              call dgemv('n', vec_size, vec_size, 1.0d0, this%aodens, vec_size, this%dbvec(:,b), 1, 0.0d0, this%dendb, 1)
#else
              this%dendb=matmul(this%dbvec(:,b), this%aodens)
#endif
            end if
            dpd(b)=diapam*this%rho(b) ! diamag. contr. to J
            do m=1,3 !dm <x,y,z>

              ! we zero these out to avoid compiler warning us that these may be used uninitialized
              prsp1 = 0.0d0
              prsp2 = 0.0d0

              if (settings%use_giao) then
#ifdef HAVE_BLAS
                ! (-i)**2 = -1
                prsp1 = -ddot(vec_size, this%dendb, 1, this%drvec(:, m), 1)
                prsp2 = ddot(vec_size, this%denbf, 1, this%d2fvec(:, k), 1)
#else
                prsp1=-dot_product(this%dendb, this%drvec(:,m)) ! (-i)**2=-1
                prsp2=dot_product(this%denbf, this%d2fvec(:,k))
#endif
              end if
#ifdef HAVE_BLAS
              ppd = ddot(vec_size, this%pdbf, 1, this%drvec(:, m), 1)
#else
              ppd=dot_product(this%pdbf, this%drvec(:,m))
#endif
              ct(m,b)=ZETA*ppd
              if (settings%use_giao) ct(m,b)=ct(m,b)+ZETA*(prsp1+prsp2)
              k=k+1
            end do
        end do


! Calculate the total J tensor. The diamagnetic pobability density only
! contributes to the diagonal.

        ! annihilate paramagnetic contribution
        if (.not.settings%use_paramag) then
            ct=D0
            your_results_are_questionable = .true.
        end if
        ! annihilate diamagnetic  contribution
        if (.not.settings%use_diamag) then
            dpd=D0
            your_results_are_questionable = .true.
        end if

        ct(1,2)=ct(1,2)+dpd(3)
        ct(1,3)=ct(1,3)-dpd(2)
        ct(2,1)=ct(2,1)-dpd(3)
        ct(2,3)=ct(2,3)+dpd(1)
        ct(3,1)=ct(3,1)+dpd(2)
        ct(3,2)=ct(3,2)-dpd(1)

    end subroutine

end module
