!============================================================================
!
module cpl_parameters

use ice_kinds_mod

implicit none

        integer(kind=int_kind) :: nt_cells                   ! nx_global x ny_global 
                                                     ! assigned in prism_init	
        integer(kind=int_kind), parameter :: jpfldout = 16   ! total number of fields sent
        integer(kind=int_kind), parameter :: jpfldin  = 17   ! total number of fields rcvd 

        integer(kind=int_kind), parameter :: n_a2i = 10      ! number of a2i fields
        integer(kind=int_kind), parameter :: n_o2i = 9       ! number of o2i fields
        integer(kind=int_kind), parameter :: n_i2a = 1       ! number of i2a fields
        integer(kind=int_kind), parameter :: n_i2o = 17      ! number of i2o fields

!
character(len=8), dimension(jpfldout) :: cl_writ ! Symb names fields sent
character(len=8), dimension(jpfldin)  :: cl_read ! Symb names fields rcvd
        integer(kind=int_kind) :: il_out                 ! format io unit(s) for coupling cpu(s)
!
 
        integer(kind=int_kind) :: num_cpl_ai    ! num of (a2i) cpl periods for this run
        integer(kind=int_kind) :: num_cpl_io    ! num of (i2o) cpl periods each dt_cpl_ai
        integer(kind=int_kind) :: num_ice_io    ! ice time loop iteration number per dt_cpl_io

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
   chk_gfdl_roughness = .false.      !.t. output u_star & roughness once a cpl interval (jan2010)

        integer(kind=int_kind) :: dt_cpl_ai = 21600    !atm<==>ice coupling interval (seconds) 
        integer(kind=int_kind) :: dt_cpl_io = 3600    !ice<==>ocn coupling interval (seconds)

    real(kind=dbl_kind) :: precip_factor = 1.0   !test the precip (temporary use)

namelist/coupling_nml/       &
         dt_cpl_ai,      &
         dt_cpl_io,      &
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
         chk_o2i_fields

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

subroutine get_cpl_timecontrol

use ice_exit
use ice_fileunits
use ice_calendar, only : dt, npt

implicit none

    integer (int_kind) :: nml_error       ! namelist read error flag

! all processors read the namelist--

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

! * make sure runtime is mutliple of dt_cpl_ai, dt_cpl_ai is mutliple of dt_cpl_io, 
! * and dt_cpl_io is mutliple of dt_cice!
num_cpl_ai = (npt*dt)/dt_cpl_ai
num_cpl_io = dt_cpl_ai/dt_cpl_io
num_ice_io = dt_cpl_io/dt

coef_ic = real(dt)/real(dt_cpl_io)

end subroutine get_cpl_timecontrol

end module cpl_parameters
