!============================================================================
!
module cpl_parameters

use ice_kinds_mod

implicit none

    integer(kind=int_kind) :: nt_cells                   ! nx_global x ny_global 

    integer, parameter :: MAX_COUPLING_FIELDS = 32
    character(len=16), dimension(MAX_COUPLING_FIELDS) :: fields_from_atm
    character(len=16), dimension(MAX_COUPLING_FIELDS) :: fields_from_ocn
    character(len=16), dimension(MAX_COUPLING_FIELDS) :: fields_to_ocn

    integer :: num_fields_from_atm
    integer :: num_fields_from_ocn
    integer :: num_fields_to_ocn

    integer(kind=int_kind) :: il_out            ! format io unit(s) for coupling cpu(s)

    integer(kind=int_kind) :: num_cpl_ai    ! num of (a2i) cpl periods for this run
    integer(kind=int_kind) :: num_cpl_io    ! num of (i2o) cpl periods each atm_ice_timestep
    integer(kind=int_kind) :: num_ice_io    ! ice time loop iteration number per ice_ocean_timestep

    real(kind=dbl_kind) :: meltlimit = 50.  !12/03/2008: set max melt
    real(kind=dbl_kind) :: ocn_albedo = 0.06

    logical :: &                         !pop_icediag is as that for ocn model, if true
       pop_icediag    = .false. , &      !    ice formation from ocn is via POP approach 
       use_ocnslope   = .false. , &      !if .t. use the sea srf tilt passed from ocn
       use_umask      = .false. , &      !if .t. use the pre-processed umask (be careful!)
       ice_pressure_on = .true. , &
       ice_fwflux     = .false. , &
       rotate_winds   = .false. , &      !.t. if oasis sends U,V as scalars. 20090319 
       limit_icemelt  = .false. , &      !.f. no limit to ice melt .         20090320
       use_core_nyf_runoff = .false. , &      !.t. use core Normal Year Forcing runoff data (remapped) 20090718
       use_core_iaf_runoff = .false. , &      !.t. use core Inter-Annual Forcing runoff data (remapped) 20120302
       cst_ocn_albedo = .true.  , &      !.t. use constant ocean albedo (e.g., 0.06, to 0.1)
       chk_frzmlt_sst = .false. , &      !      otherwise use alfa = 0.069 - 0.011 cos(2phi)
       chk_a2i_fields = .false. , &      !      as in Large & Yeager (2009).
       chk_i2a_fields = .false. , &
       chk_i2o_fields = .false. , &
       chk_o2i_fields = .false. , &
       gfdl_surface_flux = .true., &     !.t. use gfdl ocean surface flux calculation (dec2009)
       chk_gfdl_roughness = .false., &      !.t. output u_star & roughness once a cpl interval (jan2010)
       debug_output = .false.

    real(kind=dbl_kind) :: precip_factor = 1.0   !test the precip (temporary use)

    character(len=1024) :: accessom2_config_dir = '../'

    namelist/coupling_nml/       &
             pop_icediag,    &
             use_ocnslope,   &
             use_umask,      &
             rotate_winds,   &
             ice_pressure_on, &
             ice_fwflux,     &
             limit_icemelt,  &  
             meltlimit,      &
             precip_factor,  &
             cst_ocn_albedo, &
             ocn_albedo,     &
             gfdl_surface_flux, &
             chk_gfdl_roughness, &
             chk_frzmlt_sst, &
             use_core_nyf_runoff, &
             use_core_iaf_runoff, &
             chk_a2i_fields, &
             chk_i2a_fields, &
             chk_i2o_fields, &
             chk_o2i_fields, &
             accessom2_config_dir, &
             debug_output, &
             fields_from_atm, &
             fields_to_ocn, &
             fields_from_ocn

    integer(kind=int_kind) :: iniday, inimon, iniyear   !from inidate
    real(kind=dbl_kind) :: coef_ic    !dt_ice/dt_cpl_io, for i2o fields tavg 
    real(kind=dbl_kind) :: frazil_factor = 0.5

    ! frazil_factor is associated with the difference between ocean 
    ! model and ice model time-stepping: for mom4, two-level frog-leap
    ! is used and frazil heat flux is calculated and accumulated with
    ! frazil_factor = 1, which is supposed to be used for a ice model
    ! with the same two-level time-stepping scheme such as SIS. but 
    ! cice uses forward time-stepping, which means we need 'correct'
    ! the received frazil energy by multiplying 0.5...
!---------------------------------------------------------------------------------------

contains

subroutine read_namelist_parameters()

    use ice_exit
    use ice_fileunits

    integer (int_kind) :: nml_error, i

    do i=1, MAX_COUPLING_FIELDS
        fields_from_atm(i) = char(0)
        fields_from_ocn(i) = char(0)
        fields_to_ocn(i) = char(0)
    enddo

    ! all processors read the namelist
    call get_fileunit(nu_nml)
    open(unit=nu_nml,file="input_ice.nml",form="formatted",status="old",iostat=nml_error)

    if (nml_error /= 0) then
       nml_error = -1
    else
       nml_error =  1
    endif
    do while (nml_error > 0)
       read(nu_nml, nml=coupling_nml,iostat=nml_error)
       if (nml_error > 0) read(nu_nml,*)  ! for Nagware compiler
    end do
    if (nml_error == 0) close(nu_nml)

    call release_fileunit(nu_nml)

    if (nml_error /= 0) then
        call abort_ice('ice: error reading coupling_nml')
    endif

    num_fields_from_atm = 0
    do i=1, MAX_COUPLING_FIELDS
        if (fields_from_atm(i) /= CHAR(0)) then
            num_fields_from_atm = num_fields_from_atm + 1
        else
            exit
        endif
    enddo

    num_fields_from_ocn = 0
    do i=1, MAX_COUPLING_FIELDS
        if (fields_from_ocn(i) /= CHAR(0)) then
            num_fields_from_ocn = num_fields_from_ocn + 1
        else
            exit
        endif
    enddo

    num_fields_to_ocn = 0
    do i=1, MAX_COUPLING_FIELDS
        if (fields_to_ocn(i) /= CHAR(0)) then
            num_fields_to_ocn = num_fields_to_ocn + 1
        else
            exit
        endif
    enddo


endsubroutine read_namelist_parameters

subroutine get_cpl_timecontrol(atm_ice_timestep, ice_ocean_timestep)

    use ice_calendar, only : dt, npt

    integer, intent(in) :: atm_ice_timestep, ice_ocean_timestep

    ! make sure runtime is mutliple of atm_ice_timestep, atm_ice_timestep
    ! is mutliple of ice_ocean_timestep, and ice_ocean_timestep
    ! is mutliple of dt_cice!
    num_cpl_ai = (npt*dt)/atm_ice_timestep
    num_cpl_io = atm_ice_timestep/ice_ocean_timestep
    num_ice_io = ice_ocean_timestep/dt

    coef_ic = real(dt)/real(ice_ocean_timestep)

end subroutine get_cpl_timecontrol

end module cpl_parameters
