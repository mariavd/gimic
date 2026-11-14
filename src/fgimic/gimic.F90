!
! GIMIC - a pretty advanced 'Hello World!' program.
!

#include "config.h"

program gimic
    use globals_module
    use getkw_class
    use settings_module
    use teletype_module
    use basis_class
    use timer_module
    use magnet_module
    use parallel_module
    use cao2sao_class
    use grid_class
    use basis_class
    use dens_class
    use jtensor_class
    use jtensor_offload
    use jfield_class
    use divj_class
    use caos_module
    use gaussint_module
    ! acid and jav stuff
    use acid_module
    !$ use omp_lib
    implicit none

    type(molecule_t) :: mol
    type(getkw_t) :: input
    type(grid_t) :: grid
    type(jfield_t) :: jf
    type(dens_t) :: xdens
    real(DP), dimension(3) :: magnet

    integer(I4) :: nthreads

    call start_klocka()
    call new_getkw(input)
    call set_debug_level(3)

    mpi_rank = init_mpi(settings%is_mpirun)

    call program_header
    ! the team size actually obtained, not the requested maximum
    nthreads=1
    !$omp parallel
    !$omp master
    !$ nthreads=omp_get_num_threads()
    !$omp end master
    !$omp end parallel
    !$ write(str_g, '(a,i0)') 'OpenMP threads: ', nthreads
    !$ call msg_out(str_g)
    call initialize()
    call driver()
    call finalize()

    call stockas_klocka()
    call msg_out('Hello World! (tm)')
    call nl
    call msg_out('done.')
    call nl

    if (your_results_are_questionable) then
        call msg_error('Your results are questionable.')
        call nl
    end if

contains
    subroutine initialize()
        integer(I4) :: hostnm, ierr, i
        character(BUFLEN) :: fdate, sys
        character(BUFLEN) :: vti_format


        call getkw(input, 'debug', debug_level)
        call set_debug_level(debug_level)

        call getkw(input, 'title', settings%title)
        call getkw(input, 'basis', settings%basis)
        call getkw(input, 'xdens', settings%xdens)
        call getkw(input, 'density', settings%density)
        call getkw(input, 'magnet_axis', settings%magnet_axis)
        call getkw(input, 'magnet', settings%magnet)
        call getkw(input, 'openshell', settings%is_uhf)
        call getkw(input, 'dryrun', settings%dryrun)
        call getkw(input, 'calc', settings%calc)
        call getkw(input, 'xdens', settings%xdens)
        call getkw(input, 'density', settings%density)
        call getkw(input, 'mofile', settings%mofile)
        call getkw(input, 'mos', settings%morange)

        settings%use_spherical=.false.
        settings%screen_thrs = SCREEN_THRS
        call getkw(input, 'Advanced.spherical', settings%use_spherical)
        call getkw(input, 'Advanced.GIAO', settings%use_giao)
        call getkw(input, 'Advanced.diamag', settings%use_diamag)
        call getkw(input, 'Advanced.paramag',settings%use_paramag)
        call getkw(input, 'Advanced.screening', settings%use_screening)
        call getkw(input, 'Advanced.screening_thrs', settings%screen_thrs)
        call getkw(input, 'Advanced.divj_step', settings%divj_step)
        vti_format = 'binary'
        call getkw(input, 'Advanced.vti_format', vti_format)
        ! not case sensitive
        do i=1,len_trim(vti_format)
            if (vti_format(i:i) >= 'A' .and. vti_format(i:i) <= 'Z') then
                vti_format(i:i) = achar(iachar(vti_format(i:i)) + 32)
            end if
        end do
        select case (trim(adjustl(vti_format)))
        case ('binary')
            settings%vti_ascii = .false.
        case ('ascii')
            settings%vti_ascii = .true.
        case default
            call msg_error('vti_format must be binary or ascii, not ' // &
                trim(vti_format))
            stop 1
        end select

        call getkw(input, 'Essential.acid', settings%acid)
        call getkw(input, 'Essential.jmod', settings%jmod)
        call getkw(input, 'Essential.prop', settings%prop)
        call getkw(input, 'Essential.tropicity', settings%tropicity)
        call getkw(input, 'Tropicity.max_step', settings%trop_step)
        call getkw(input, 'Tropicity.accuracy', settings%trop_accuracy)
        call getkw(input, 'Tropicity.tolerance', settings%trop_tol)
        call getkw(input, 'Tropicity.max_length', settings%trop_length)
        call getkw(input, 'Tropicity.capture', settings%trop_capture)
        if (settings%tropicity .and. settings%is_uhf) then
            call msg_warn('The tropicity is not available for open-shell &
                &systems, skipped')
            settings%tropicity=.false.
        end if

        ierr=hostnm(sys)
        if (mpi_rank == 0) then
            write(str_g, '(a,i3,a,a)') 'MPI master', mpi_rank,' on ', trim(sys)
            call msg_debug(str_g,2)
        else
            call set_teletype_unit(DEVNULL)
            write(str_g, '(a,i3,a,a)') 'MPI slave', mpi_rank,' on ', trim(sys)
            call msg_debug(str_g,2)
        end if

        call msg_out(fdate())

        call msg_out(' TITLE: '// trim(settings%title))
        call nl

        if (.not.settings%use_giao) then
            call msg_info('GIAOs not used!')
            call nl
        end if

        if (.not.settings%use_diamag) then
            call msg_info( 'Diamagnetic contributions not calculated!')
            call nl
        end if
        if (.not.settings%use_paramag) then
            call msg_info( 'Paramagnetic contributions not calculated!')
            call nl
        end if
        if ((.not.settings%use_diamag).and.(.not.settings%use_paramag)) then
            call msg_out( '    ...this does not make sense...')
            call nl
            call msg_critical( '    PLEASE SEEK PROFESSIONAL HELP, ASAP!  ')
            call nl
            stop
        end if
    end subroutine

    subroutine finalize()
        call del_basis(mol)
        call stop_mpi()
        call del_getkw(input)
    end subroutine

    subroutine driver
        type(cao2sao_t) :: c2s

        if (settings%use_screening) then
            call new_basis(mol, settings%basis, settings%screen_thrs)
        else
            call new_basis(mol, settings%basis, -1.d0)
        end if

        if (.not. settings%dryrun) then
            if (settings%use_spherical) then
                call new_c2sop(c2s,mol)
                call set_c2sop(mol, c2s)
            end if

            call new_dens(xdens, mol)
            call read_dens(xdens, settings%xdens)
            call offload_init(mol, xdens)
        end if

        call new_grid(grid, input, mol)
        if (mpi_rank == 0) then
            call plot_grid_xyz(grid, 'grid.xyz', mol)
        end if

        call get_magnet(grid, magnet)

        if (settings%is_uhf) then
            call msg_info('Open-shell calculation')
        else
            call msg_info('Closed-shell calculation')
        end if
        call nl

        if (settings%dryrun) then
            call msg_note('Dry run, not calculating ...')
            call nl
        end if

        if (settings%calc(1:5) == 'cdens') then
            call run_cdens(jf,mol,xdens)
        else if (settings%calc(1:8) == 'integral') then
            call run_integral()
        else if (settings%calc(1:4) == 'divj') then
            call run_divj(mol,xdens)
        else
            call msg_error('gimic(): Unknown operation!')
        end if

        if (.not. settings%dryrun) then
            if (settings%use_spherical) then
                call del_c2sop(c2s)
            end if
            call offload_finalize()
            call del_dens(xdens)
        end if
        call del_grid(grid)
    end subroutine

    subroutine run_cdens(jf,mol,xdens)
        type(jfield_t) :: jf
        type(dens_t) :: xdens
        type(molecule_t) :: mol
        call msg_out('Calculating current density')
        call msg_out('*****************************************')
        call new_jfield(jf, grid, magnet)

        if (settings%dryrun) return

        call calc_jvectors(jf, mol, xdens)
        call jvector_plots(jf)
        if (settings%tropicity) call tropicity_plots(jf, mol, xdens)

        if (settings%is_uhf) then
            call calc_jvectors(jf, mol, xdens, 'alpha')
            call jvector_plots(jf, 'alpha')

            call calc_jvectors(jf, mol, xdens, 'beta')
            call jvector_plots(jf, 'beta')

            call calc_jvectors(jf, mol, xdens, 'spindens')
            call jvector_plots(jf, 'spindens')
        endif
        call del_jfield(jf)
    end subroutine

    subroutine run_divj(mol,xdens)
        type(dens_t) :: xdens
        type(molecule_t) :: mol
        type(divj_t) :: dj

        call msg_out('Calculating divergence of the current density')
        call msg_out('*****************************************')
        call new_divj(dj, grid, magnet)

        if (.not. settings%dryrun) then
            call calc_divj(dj, mol, xdens)
            call divj_plots(dj)

            if (settings%is_uhf) then
                call calc_divj(dj, mol, xdens, 'alpha')
                call divj_plots(dj, 'alpha')

                call calc_divj(dj, mol, xdens, 'beta')
                call divj_plots(dj, 'beta')

                call calc_divj(dj, mol, xdens, 'spindens')
                call divj_plots(dj, 'spindens')
            end if
        end if
        call del_divj(dj)
    end subroutine

    subroutine run_integral()
        use integral_class
        type(integral_t) :: it
        type(cao2sao_t) :: c2s
        logical :: have_dens

        call msg_out('Integrating current density')
        call msg_out('*****************************************')
        call new_integral(it, grid)

        if (settings%dryrun) then
            ! estimate the run time when the density is there to be read
            inquire(file=settings%xdens, exist=have_dens)
            if (have_dens) then
                if (settings%use_spherical) then
                    call new_c2sop(c2s,mol)
                    call set_c2sop(mol, c2s)
                end if
                call new_dens(xdens, mol)
                call read_dens(xdens, settings%xdens)
                call integral_eta(grid, mol, xdens)
                call del_dens(xdens)
                if (settings%use_spherical) call del_c2sop(c2s)
            end if
            return
        end if

        call integral_eta(grid, mol, xdens)

        if (settings%jmod) then
          call msg_note('Integrating |J|')
          call integrate_modulus(it, mol, xdens)
          if (settings%is_uhf) then
              call integrate_modulus(it, mol, xdens, 'alpha')
              call integrate_modulus(it, mol, xdens, 'beta')
              call integrate_modulus(it, mol, xdens, 'spindens')
          end if
          call nl
        else
          write(*,*) "Jmod integration skipped."
        end if

        call msg_note('Integrating current')
        call integrate_current(it, mol, xdens)
        if (settings%is_uhf) then
            call integrate_current(it, mol, xdens, 'alpha')
            call integrate_current(it, mol, xdens, 'beta')
            call integrate_current(it, mol, xdens, 'spindens')
        end if
        call nl
        if (settings%acid) then
        ! add here ACID integration !
           call msg_note('Integrating ACID density')
           call integrate_acid(it, mol, xdens)
           call nl
        end if

        call del_integral(it)
    end subroutine

    subroutine program_header
call nl
call msg_out('****************************************************************')
call msg_out('***                                                          ***')
call msg_out('***           GIMIC '// PROJECT_VERSION // &
                       ' (' // GIT_HASH // ')                          ***')
call msg_out('***           Written by Jonas Juselius and co-authors       ***')
call msg_out('***                                                          ***')
call msg_out('***  This software is copyright (c) 2003-2023 by             ***')
call msg_out('***  GIMIC authors.                                          ***')
call msg_out('***                                                          ***')
call msg_out('***  You are free to distribute this software under the      ***')
call msg_out('***  terms of the GNU General Public License                 ***')
call msg_out('****************************************************************')
call nl
    end subroutine

end program

! vim:et:sw=4:ts=4
