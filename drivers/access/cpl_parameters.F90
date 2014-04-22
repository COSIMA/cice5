!============================================================================
!
module cpl_parameters
!
!----------------------------------------------------------------------------

use ice_kinds_mod

implicit none

integer(kind=int_kind) :: il_im, il_jm, il_imjm    ! il_im=nx_global, il_jm=ny_global 
                                                   ! assigned in prism_init
integer (kind=int_kind) :: xdim, ydim
!integer(kind=int_kind), parameter :: nsend = 50   ! maxium no of flds sent allowed
!integer(kind=int_kind), parameter :: nrecv = 50   ! maxium no of flds rcvd allowed
integer(kind=int_kind) :: nsend_i2a, nsend_i2o
integer(kind=int_kind) :: nrecv_a2i, nrecv_o2i 
integer(kind=int_kind), parameter :: jpfldout = 37 ! actual number of fields sent
integer(kind=int_kind), parameter :: jpfldin  = 35 ! actual umber of fields rcvd 

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
integer(kind=int_kind) :: num_cpl_io    ! num of (i2o) cpl periods
integer(kind=int_kind) :: num_ice_io    ! ice time loop iter num per i2o cpl interval

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
integer(kind=int_kind) :: inidate = 01010101   !beginning date of this run (yyyymmdd)
integer(kind=int_kind) :: init_date = 00010101 !beginning date of this EXP (yyyymmdd)
integer(kind=int_kind) :: dt_cice = 3600       !time step of this model      (seconds) 
integer(kind=int_kind) :: dt_cpl_ai = 21600    !atm<==>ice coupling interval (seconds) 
integer(kind=int_kind) :: dt_cpl_io = 21600    !ice<==>ocn coupling interval (seconds)
integer(kind=int_kind) :: caltype = 0          !calendar type: 0 (365daye/yr, 'Juilian' ) 
                                               ! 1 (365/366 days/yr, 'Gregorian')
                                               ! n (n days/month)
!integer(kind=int_kind) :: runtime0    !accumulated run time by the end of last run (s)   
real(kind=dbl_kind) :: runtime0 = 0.0  !  can be too large as int to read in correctly!
integer(kind=int_kind) :: runtime = 86400      !the time length for this run segment (s)

!20100305: Harry Henden suggests turning off ocean current into UM might reduce the 
!          tropical cooling bias:
real(kind=dbl_kind) :: ocn_ssuv_factor = 1.0  ! 0.0 -- turn off the ocn_current into UM.
real(kind=dbl_kind) :: iostress_factor = 1.0  ! 0.0 -- turn off stresses into MOM4.
             
namelist/coupling/       &
         caltype,        &
         jobnum,         &
         inidate,        &
         init_date,      &
         runtime0,       &   
         runtime,        &
         dt_cice,        &
         dt_cpl_ai,      &
         dt_cpl_io,      &
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
         ocn_ssuv_factor,&
         iostress_factor,&
         chk_a2i_fields, &
         chk_i2a_fields, &
         chk_i2o_fields, &
         chk_o2i_fields

integer(kind=int_kind) :: iniday, inimon, iniyear   !from inidate
real(kind=dbl_kind) :: coef_io    !dt_ice/dt_cpl_io, for i2o fields tavg 
real(kind=dbl_kind) :: coef_ia    !dt_ice/dt_cpl_ai, for i2a fields tavg
real(kind=dbl_kind) :: coef_cpl   !dt_cpl_io/dt_cpl_ai, for ocn fields tavg 

real(kind=dbl_kind) :: frazil_factor = 0.5
         !frazil_factor is associated with the difference between ocean 
         ! model and ice model time-stepping: for mom4, two-level frog-leap
         ! is used and frazil heat flux is calculated and accumulated with
         ! frazil_factor = 1, which is supposed to be used for a ice model
         ! with the same two-level time-stepping scheme such as SIS. but 
         ! cice uses forward time-stepping, which means we need 'correct'
         ! the received frazil energy by multiplying 0.5...
!---------------------------------------------------------------------------------------

contains

!=======================================================================================
subroutine get_cpl_timecontrol_simple

implicit none

! all processors read the namelist--
open(unit=99,file="input_ice.nml",form="formatted",status="old")
read (99, coupling)
close(unit=99)
! *** make sure dt_cpl_ai is multiple of dt_cpl_io, and dt_cpl_io if multiple of dt_ice ***
num_cpl_ai = runtime/dt_cpl_ai
num_cpl_io = dt_cpl_ai/dt_cpl_io
num_ice_io = dt_cpl_io/dt_cice

coef_io = float(dt_cice)/float(dt_cpl_io)
coef_ia = float(dt_cice)/float(dt_cpl_ai)
coef_cpl = float(dt_cpl_io)/float(dt_cpl_ai)

iniday  = mod(inidate, 100)
inimon  = mod( (inidate - iniday)/100, 100)
iniyear = inidate / 10000

return
end subroutine get_cpl_timecontrol_simple

!===============================================================================
subroutine get_cpl_timecontrol

use ice_exit
use ice_fileunits

implicit none

integer (int_kind) :: nml_error       ! namelist read error flag

! all processors read the namelist--

call get_fileunit(nu_nml)
open(unit=nu_nml,file="input_ice.nml",form="formatted",status="old",iostat=nml_error)
!
write(6,*)'CICE: input_ice.nml opened at unit = ', nu_nml
!
if (nml_error /= 0) then
   nml_error = -1
else
   nml_error =  1
endif
do while (nml_error > 0)
   read(nu_nml, nml=coupling,iostat=nml_error)
   if (nml_error > 0) read(nu_nml,*)  ! for Nagware compiler
end do
if (nml_error == 0) close(nu_nml)

write(6,coupling)

call release_fileunit(nu_nml)

if (nml_error /= 0) then
   !!!call abort_ice('ice: error reading coupling')
   write(6, *)
   write(6, *)'XXX Warning: after reading coupling, nml_error = ',nml_error
   write(6, *)
endif

! * make sure runtime is mutliple of dt_cpl_ai, dt_cpl_ai is mutliple of dt_cpl_io, 
! * and dt_cpl_io is mutliple of dt_cice!
num_cpl_ai = runtime/dt_cpl_ai
num_cpl_io = dt_cpl_ai/dt_cpl_io
num_ice_io = dt_cpl_io/dt_cice

coef_io = float(dt_cice)/float(dt_cpl_io)
coef_ia = float(dt_cice)/float(dt_cpl_ai)
coef_cpl = float(dt_cpl_io)/float(dt_cpl_ai)

iniday  = mod(inidate, 100)
inimon  = mod( (inidate - iniday)/100, 100)
iniyear = inidate / 10000

return
end subroutine get_cpl_timecontrol

!=======================================================================================
end module cpl_parameters
