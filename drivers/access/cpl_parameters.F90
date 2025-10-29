!============================================================================
!
module cpl_parameters
!
!----------------------------------------------------------------------------

use ice_kinds_mod

implicit none

#ifdef __INTEL_COMPILER
! for intel runtime errors
! see https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2025-2/list-of-runtime-error-messages.html
include "for_iosdef.for"
#endif

integer(kind=int_kind) :: il_im, il_jm, il_imjm    ! il_im=nx_global, il_jm=ny_global 
                                                   ! assigned in prism_init
integer (kind=int_kind) :: xdim, ydim
!integer(kind=int_kind), parameter :: nsend = 50   ! maxium no of flds sent allowed
!integer(kind=int_kind), parameter :: nrecv = 50   ! maxium no of flds rcvd allowed
integer(kind=int_kind) :: nsend_i2a, nsend_i2o
integer(kind=int_kind) :: nrecv_a2i, nrecv_o2i 
integer(kind=int_kind), parameter :: jpfldout = 39 ! actual number of fields sent
integer(kind=int_kind), parameter :: jpfldin  = 35 ! actual number of fields rcvd 

character(len=8), dimension(jpfldout) :: cl_writ ! Symb names fields sent
character(len=8), dimension(jpfldin)  :: cl_read ! Symb names fields rcvd

integer(kind=int_kind), dimension(jpfldout) :: il_var_id_out ! ID for fields sent 
integer(kind=int_kind), dimension(jpfldin)  :: il_var_id_in  ! ID for fields rcvd

character(len=6), parameter :: cp_modnam='cicexx' ! Component model name

integer(kind=int_kind) :: il_out        ! format io unit(s) for coupling cpu(s)
integer(kind=int_kind) :: il_commlocal  ! Component internal communicator 
!integer(kind=int_kind) :: il_rank      ! local procesor id 
!integer(kind=int_kind) :: il_master=0	! master_task id
 
integer(kind=int_kind) :: num_cpl_ai    ! num of (a2i) cpl periods
integer(kind=int_kind) :: num_ice_ai    ! ice time loop iter num per a2i cpl interval

real(kind=dbl_kind) :: meltlimit = -200.    !12/03/2008: set max melt
real(kind=dbl_kind) :: ocn_albedo = 0.06 ! for compability with AusCOM
character(len=256) :: inputdir   = 'INPUT'
character(len=256) :: restartdir = 'RESTART'
logical :: &                         !pop_icediag is as that for ocn model, if true
   pop_icediag    = .true.  , &      !    ice formation from ocn is via POP approach 
   use_ocnslope   = .false. , &      !if .t. use the sea srf tilt passed from ocn
   use_umask      = .false. , &      !if .t. use the pre-processed umask (be careful!)
   ice_pressure_on = .true. , &
   air_pressure_on = .false., &
   ice_fwflux     = .false. , &
   rotate_winds   = .false. , &      !.t. if oasis sends U,V as scalars. 20090319
   limit_icemelt  = .false. , &
   limit_stflx    = .false. , &       !.t. set limit for the salt flux to ocn (switch 20111108)
   use_core_runoff = .true. , &      !.t. use core runoff data (remapped) 20090718
   cst_ocn_albedo = .true.  , &      !.t. use constant ocean albedo (e.g., 0.06, to 0.1)
   gbm2pericearea = .true.  , &      !.t. do GBM to per ice area conversion in set_sfcflux
   do_scale_fluxes = .true. , &      !.t. call scale_fluxes in routine coupling_prep.
   extreme_test   = .false. , &      !.t. output extreme forcing data (in set_sfcflux)
   imsk_evap      = .true.  , &
   chk_a2i_fields = .false. , &
   chk_i2a_fields = .false. , &
   chk_i2o_fields = .false. , &
   chk_o2i_fields = .false.
integer(kind=int_kind) :: jobnum = 1           !1 for initial, >1 restart
integer(kind=int_kind) :: init_date = 00010101 !beginning date of this EXP (yyyymmdd)
integer(kind=int_kind) :: iniday = 1, &        ! beginning date of this run. Read from restart
                          inimon = 1, &
                          iniyear = 1
integer(kind=int_kind) :: dt_cice = 3600       !time step of this model      (seconds) 
integer(kind=int_kind) :: dt_cpl_ai = 21600    !atm<==>ice coupling interval (seconds) 
integer(kind=int_kind) :: dt_cpl_io = -99      !ice<==>ocn coupling interval (seconds).
                                               !Hardwired to equal dt_cice and should not
                                               !be set in namelist.
real(kind=dbl_kind)    :: runtime0 = 0.0       !accumulated runtime from init_date to
                                               !run start date
integer(kind=int_kind) :: runtime = 86400      !the time length for this run segment (s)

!20100305: Harry Henden suggests turning off ocean current into UM might reduce the 
!          tropical cooling bias:
real(kind=dbl_kind) :: ocn_ssuv_factor = 1.0  ! 0.0 -- turn off the ocn_current into UM.
real(kind=dbl_kind) :: iostress_factor = 1.0  ! 0.0 -- turn off stresses into MOM4.
!
!20171227: Adding options for land ice discharge as iceberg melt (0,1,2,3,4)
integer(kind=int_kind) :: iceberg = 2 
!Allow scaling: factor for the iceberg waterflux (won't change water mass budget)
real(kind=dbl_kind) :: &
         iceberg_rate_s = 0.5, &      !rate of "iceberg" taken from the runoff off Antarctica
         iceberg_rate_n = 0.5, &      !........................................... Greenland        
         iceberg_lh = 1.0             !iceberg latent heat (=0 if CABLE already calculated melting)
             
logical :: runoff_lh = .true.         !allow runoff to carry LH when discharged into ocean
                                      !which would lead to ocean surface cooling,
                                      !when .false., only carry LH to areas where 
                                      !runoff spread by the lice (iceberg) mask
integer(kind=int_kind) :: &
         iceberg_je_s = 70, &   !(iceberg_js_s=1, always)
         runoff_je_s  = 45, &   !(runoff_js_s =1, always)
         iceberg_js_n = 201, &  !(iceberg_je_n=300, always)
         runoff_is_n  = 222, &  !------
         runoff_ie_n  = 270, &  !These 4 indices define the
         runoff_js_n  = 230, &  !Greenland runoff domain
         runoff_je_n  = 300     !-----
!202412: add option for "fixing" ocean water mass imbalance: ESM1.5 sees ~ 0.18543417E+08 kg/s
!        (annual mean) water loss from ocean in a 100-year test run (liceA0G0), meaning a net
!        loss rate of 0.18543417E+08/0.36133599E+15(ocean-surface-area) = 0.513190E-07 kg/m2/s.
!        which will be compensated for by adding this much of waterflux to global lprec field--
real(kind=dbl_kind) :: &
         add_lprec = 0.513190E-07      !kg/m2/s. ==> set to 0.0 if no fixin!
             
namelist/coupling/       &
         jobnum,         &
         runtime,        &
         dt_cice,        &
         dt_cpl_ai,      &
         inputdir,       &
         restartdir,     &
         pop_icediag,    &
         use_ocnslope,   &
         use_umask,      &
         ice_pressure_on, &
         air_pressure_on, &
         ice_fwflux,     &
         rotate_winds,   &
         limit_icemelt,  &
         limit_stflx,    &
         meltlimit,      &
         use_core_runoff, &
         gbm2pericearea,  &
         do_scale_fluxes, &
         extreme_test,   &
         imsk_evap,      &
         iceberg,        &
         iceberg_rate_s, &
         iceberg_rate_n, &
         iceberg_lh,   &
         iceberg_je_s, &
         runoff_je_s,  &
         iceberg_js_n, &
         runoff_is_n,  &
         runoff_ie_n,  &
         runoff_js_n,  &
         runoff_je_n,  &
         runoff_lh,    &
         add_lprec,    &
         ocn_ssuv_factor,&
         iostress_factor,&
         chk_a2i_fields, &
         chk_i2a_fields, &
         chk_i2o_fields, &
         chk_o2i_fields

real(kind=dbl_kind) :: coef_ai    !dt_ice/dt_cpl_ai, for i2a fields tavg

real(kind=dbl_kind) :: frazil_factor = 0.5
         !frazil_factor is associated with the difference between ocean 
         ! model and ice model time-stepping: for mom4, two-level frog-leap
         ! is used and frazil heat flux is calculated and accumulated with
         ! frazil_factor = 1, which is supposed to be used for a ice model
         ! with the same two-level time-stepping scheme such as SIS. but 
         ! cice uses forward time-stepping, which means we need 'correct'
         ! the received frazil energy by multiplying 0.5...
!---------------------------------------------------------------------------------------

logical :: newstep_ai = .false.         !20171024: for land ice availiblity control

contains

!=======================================================================================
subroutine get_cpl_timecontrol_simple

implicit none

! all processors read the namelist--
open(unit=99,file="input_ice.nml",form="formatted",status="old")
read (99, coupling)
close(unit=99)

! *** make sure dt_cpl_ai is multiple of dt_cpl_io, and dt_cpl_io if multiple of dt_ice ***

!hardrwire dt_cpl_io == dt_cice
dt_cpl_io = dt_cice

num_cpl_ai = runtime/dt_cpl_ai
num_ice_ai = dt_cpl_ai/dt_cice

coef_ai = float(dt_cice)/float(dt_cpl_ai)

return
end subroutine get_cpl_timecontrol_simple

!===============================================================================

subroutine get_cpl_timecontrol
use ice_exit, only: abort_ice
use ice_fileunits, only: nu_nml, ice_stderr, ice_stdout, get_fileunit, release_fileunit
use ice_communicate, only: my_task, master_task

implicit none

integer (int_kind) :: nml_error       ! namelist read error flag
character (len=256) :: errstr, tmpstr         ! For holding namelist read errors

! all processors read the namelist--

call get_fileunit(nu_nml)
open(unit=nu_nml,file="input_ice.nml",form="formatted",status="old",iostat=nml_error, iomsg=errstr)
!
if (my_task == master_task) then
   write(ice_stdout,*)'CICE: input_ice.nml opened at unit = ', nu_nml
endif
!
if (nml_error /= 0) then
   write(tmpstr, '(a,i3,a)') 'CICE: ERROR failed to open input_ice.nml. Error code: ', nml_error, &
                            '  - ' // trim(errstr)
   call abort_ice(trim(tmpstr))
else
   nml_error =  1
endif

do while (nml_error > 0)
   read(nu_nml, nml=coupling,iostat=nml_error,iomsg=errstr)
   ! check if error
   if (nml_error /= 0) then
      if (my_task == master_task) then
         ! backspace and re-read erroneous line
         backspace(nu_nml)
         read(nu_nml,fmt=*) tmpstr
#ifdef __INTEL_COMPILER
         if (nml_error == FOR$IOS_INVREFVAR) then
           write(ice_stderr,*)'CICE: Invalid reference to variable '//trim(tmpstr)
           write(ice_stderr,*)'CICE: is '//trim(tmpstr)//' deprecated ?'
         endif
#endif
         call abort_ice('CICE ERROR in input_ice.nml when' // &
            ' reading ' // trim(tmpstr) // ' - ' //errstr)
       endif
   endif
end do
if (nml_error == 0) close(nu_nml)

if (my_task == master_task) then
   write(6,coupling)
endif

call release_fileunit(nu_nml)

if (nml_error /= 0) then
   if (my_task == master_task) then
      call abort_ice('ice: error reading coupling namelist in "input_ice.nml"')
   endif
endif

!hardwire dt_cpl_io == dt_cice
dt_cpl_io = dt_cice

! * make sure runtime is mutliple of dt_cpl_ai, dt_cpl_ai is mutliple of dt_cpl_io, 
! * and dt_cpl_io is mutliple of dt_cice!
num_cpl_ai = runtime/dt_cpl_ai
num_ice_ai = dt_cpl_ai/dt_cice

coef_ai = float(dt_cice)/float(dt_cpl_ai)

return
end subroutine get_cpl_timecontrol

!=======================================================================================
end module cpl_parameters
