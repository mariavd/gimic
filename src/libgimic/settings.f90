!
!  Tunable parameters
!
module settings_module
    use kinds_module
    implicit none

    integer(I4) :: debug_level=0

    type settings_t
        character(256) :: title
        character(128) :: debug_file
        logical :: dryrun
        logical :: show_axis
        character(16) :: calc
        character(128) :: basis
        character(128) :: xdens
        character(128) :: mofile
        character(128) :: density
        real(DP), dimension(3) :: magnet
        character(2) :: magnet_axis
        logical :: is_uhf=.false.
        logical :: is_mpirun=.false.
        logical :: use_spherical=.false.
        logical :: use_giao=.true.
        logical :: use_diamag=.true.
        logical :: use_paramag=.true.
        logical :: acid=.false.
        logical :: jmod=.false.
        logical :: prop=.false.
        ! streamline tracing for the tropicity, see tropicity_module
        logical :: tropicity=.false.
        real(DP) :: trop_step=0.5d0      ! largest step along the streamline (bohr)
        real(DP) :: trop_accuracy=1.d-5  ! error allowed in one step (bohr)
        real(DP) :: trop_tol=0.1d0       ! closure tolerance (bohr)
        real(DP) :: trop_length=100.d0   ! maximum path length (bohr)
        ! grid points within this distance of a streamline, in units of the
        ! grid spacing, take its label (0: follow every grid point)
        real(DP) :: trop_capture=0.5d0
        real(DP) :: screen_thrs
        logical :: use_screening
        integer(I4), dimension(2) :: morange
        integer(I4) :: lip_order
        real(DP) :: divj_step=1.d-3
        ! .vti files as text instead of binary (Advanced.vti_format=ascii)
        logical :: vti_ascii=.false.
    end type

    type(settings_t), save :: settings
end module

! vim:et:sw=4:ts=4
