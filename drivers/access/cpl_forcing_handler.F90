MODULE cpl_forcing_handler
!
! It contains subroutines handling coupling fields. 
!
use ice_blocks
use ice_forcing
use ice_read_write
use ice_domain_size
use ice_domain,    only : distrb_info, nblocks 
use ice_flux            !forcing data definition (Tair, Qa, uocn, etc.)
                        !Tn_top, keffn_top ...(for multilayer configuration)   
!use ice_state,     only : aice, aicen, trcr, trcrn, nt_hpnd   !ice concentration and tracers
use ice_state,     only : aice, aicen, trcr !!!, trcrn, nt_hpnd, nt_Tsfc   !ice concentration and tracers
use ice_state,     only: uvel, vvel, vsnon, vicen
use ice_gather_scatter
use ice_broadcast
use ice_constants
use ice_grid,      only : tmask, to_ugrid
use ice_communicate, only : my_task, master_task
!use ice_ocean,     only : cprho
use ice_exit,      only : abort_ice
use ice_shortwave, only : apeffn
use ice_grid, only: tarea
use ice_calendar, only: month

use cpl_parameters
use cpl_netcdf_setup
use cpl_arrays_setup

use ice_calendar, only: month

implicit none

real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
  aiiu       ! ice fraction on u-grid 

contains

!=================================================
subroutine get_core_runoff(fname, vname, nrec)
! read in the remapped core runoff data (S.Marsland) which will be used to replace
! the ncep2 runoff sent from matm via coupler 

implicit none

character*(*), intent(in) :: fname, vname
integer(kind=int_kind), intent(in) :: nrec
logical :: dbug
integer(kind=int_kind) :: ncid

dbug = .true.

if ( file_exist(fname) ) then
  call ice_open_nc(fname, ncid)
  call ice_read_global_nc(ncid, nrec, vname, gwork, dbug)
  call scatter_global(core_runoff, gwork, master_task, distrb_info, &
                      field_loc_center, field_type_scalar)
  if (my_task == master_task) call ice_close_nc(ncid)
else
  if (my_task == 0) write(il_out,*) '(get_core_runoff) file doesnt exist: ', fname
  stop 'CICE stopped: core runoff (remapped) file not found.'
endif

return
end subroutine get_core_runoff

!=================================================
subroutine get_time0_sstsss(fname, nmonth)

! This routine is to be used only once at the beginning at an exp.

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nmonth
logical :: dbug
integer(kind=int_kind) :: ncid

dbug = .true.
!dbug = .false.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(get_time0_sstsss) opening ncfile: ',fname
  endif
  call ice_open_nc(fname, ncid)
  if (my_task==0) then
    write(il_out,*) '(get_time0_sstsss) reading in initial SST...'
  endif
  call ice_read_nc(ncid, nmonth, 'TEMP', sst, dbug)
  call gather_global(gwork, sst, master_task, distrb_info)
  if (my_task==0) then
    write(il_out,*) '(get_time0_sstsss) reading in initial SSS...' 
  endif
  call ice_read_nc(ncid, nmonth, 'SALT', sss, dbug)
  call gather_global(gwork, sss, master_task, distrb_info)
  if (my_task == master_task) call ice_close_nc(ncid)
else
  if (my_task==0) then
    write(il_out,*) '(get_time0_sstsss) file doesnt exist: ', fname
  endif
  call abort_ice('CICE stopped--initial SST and SSS ncfile not found.')  
endif

return
end subroutine get_time0_sstsss

!=================================================
! temporary use ...
subroutine read_access_a2i_data(fname,nrec,istep) 

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nrec,istep
logical :: dbug
integer(kind=int_kind) :: ncid

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(read_access_a2i_data) opening ncfile: ',fname
  endif
  call ice_open_nc(fname, ncid)
  if (my_task==0) then
    write(il_out,*) '(read_access_a2i_data) reading a2i forcing data...'
  endif
  call ice_read_nc(ncid, nrec, 'thflx_i', um_thflx, dbug)
  call ice_read_nc(ncid, nrec, 'pswflx_i', um_pswflx, dbug)
  call ice_read_nc(ncid, nrec, 'runoff_i', um_runoff, dbug)
  call ice_read_nc(ncid, nrec, 'wme_i', um_wme, dbug)
  call ice_read_nc(ncid, nrec, 'rain_i', um_rain, dbug)
  call ice_read_nc(ncid, nrec, 'snow_i', um_snow, dbug)
  call ice_read_nc(ncid, nrec, 'evap_i', um_evap, dbug)
  call ice_read_nc(ncid, nrec, 'lhflx_i', um_lhflx, dbug)
  call ice_read_nc(ncid, nrec, 'tmlt01_i', um_tmlt(:,:,1,:), dbug)
  call ice_read_nc(ncid, nrec, 'tmlt02_i', um_tmlt(:,:,2,:), dbug)
  call ice_read_nc(ncid, nrec, 'tmlt03_i', um_tmlt(:,:,3,:), dbug)
  call ice_read_nc(ncid, nrec, 'tmlt04_i', um_tmlt(:,:,4,:), dbug)
  call ice_read_nc(ncid, nrec, 'tmlt05_i', um_tmlt(:,:,5,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt01_i', um_bmlt(:,:,1,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt02_i', um_bmlt(:,:,2,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt03_i', um_bmlt(:,:,3,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt04_i', um_bmlt(:,:,4,:), dbug)
  call ice_read_nc(ncid, nrec, 'bmlt05_i', um_bmlt(:,:,5,:), dbug)
  call ice_read_nc(ncid, nrec, 'taux_i', um_taux, dbug)
  call ice_read_nc(ncid, nrec, 'tauy_i', um_tauy, dbug)
  call ice_read_nc(ncid, nrec, 'swflx_i', um_swflx,  dbug)
  call ice_read_nc(ncid, nrec, 'lwflx_i', um_lwflx,  dbug)
  call ice_read_nc(ncid, nrec, 'shflx_i', um_shflx,  dbug)
  call ice_read_nc(ncid, nrec, 'press_i', um_press,  dbug)
  call ice_read_nc(ncid, nrec, 'co2_ai', um_co2,  dbug)
  call ice_read_nc(ncid, nrec, 'wnd_ai', um_wnd,  dbug)
  if (my_task == master_task) call ice_close_nc(ncid)
else
  if (my_task==0) then
    write(il_out,*) '(ed_access_a2i_data file doesnt exist: ', fname
  endif
  call abort_ice('CICE stopped--ACCESS fields_a2i ncfile not found.')
endif

call check_a2i_fields(istep)

end subroutine read_access_a2i_data

!=================================================
subroutine read_restart_i2a(fname, sec) !'i2a.nc', 0)

! read ice to atm coupling fields from restart file, and send to atm module

implicit none
character*(*), intent(in) :: fname
integer :: sec

integer(kind=int_kind) :: ncid
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(read_restart_i2a) reading in i2a fields......'
  endif
  call ice_open_nc(fname, ncid)
  call ice_read_nc(ncid, 1, 'icecon01',   ia_aicen(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'icecon02',   ia_aicen(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'icecon03',   ia_aicen(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'icecon04',   ia_aicen(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'icecon05',   ia_aicen(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk01',   ia_snown(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk02',   ia_snown(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk03',   ia_snown(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk04',   ia_snown(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'snwthk05',   ia_snown(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk01',   ia_thikn(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk02',   ia_thikn(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk03',   ia_thikn(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk04',   ia_thikn(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'icethk05',   ia_thikn(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'isst_ia',    ia_sst,   dbug)
  call ice_read_nc(ncid, 1, 'uvel_ia',    ia_uvel,  dbug)
  call ice_read_nc(ncid, 1, 'vvel_ia',    ia_vvel,  dbug)
  call ice_read_nc(ncid, 1, 'co2_i2',     ia_co2,   dbug)
  call ice_read_nc(ncid, 1, 'co2fx_i2',   ia_co2fx, dbug)

  if (my_task == master_task) then
    call ice_close_nc(ncid)
    write(il_out,*) '(read_restart_i2a) has read in 18 i2a fields.'
  endif

else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (read_restart_i2a) not found file *** ',fname
  endif
  print *, 'CICE: (read_restart_i2a) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 i2a data file.')
endif
end subroutine read_restart_i2a

!=================================================
subroutine read_restart_i2asum(fname, sec) !'i2a.nc', 0)

! read ice to atm coupling fields from restart file, and send to atm module

implicit none
character*(*), intent(in) :: fname
integer :: sec

integer(kind=int_kind) :: ncid
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(read_restart_i2asum) reading in i2a fields......'
  endif
  call ice_open_nc(fname, ncid)
  call ice_read_nc(ncid, 1, 'maicen1',   maicen(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'maicen2',   maicen(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'maicen3',   maicen(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'maicen4',   maicen(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'maicen5',   maicen(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown1',   msnown(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown2',   msnown(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown3',   msnown(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown4',   msnown(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'msnown5',   msnown(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn1',   mthikn(:,:,1,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn2',   mthikn(:,:,2,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn3',   mthikn(:,:,3,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn4',   mthikn(:,:,4,:),   dbug)
  call ice_read_nc(ncid, 1, 'mthikn5',   mthikn(:,:,5,:),   dbug)
  call ice_read_nc(ncid, 1, 'msst',      msst,   dbug)
  call ice_read_nc(ncid, 1, 'mssu',      mssu,   dbug)
  call ice_read_nc(ncid, 1, 'mssv',      mssv,   dbug)
  call ice_read_nc(ncid, 1, 'muvel',     muvel,  dbug)
  call ice_read_nc(ncid, 1, 'mvvel',     mvvel,  dbug)
  call ice_read_nc(ncid, 1, 'maiu',      maiu,   dbug)
  !
  !call ice_read_nc(ncid, 1, 'maice_ia', maice_ia,          dbug)

  if (my_task == master_task) then
    call ice_close_nc(ncid)
    write(il_out,*) '(read_restart_i2asum) has read in 21 i2a fields.'
  endif

else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (read_restart_i2asum) not found file *** ',fname
  endif
  print *, 'CICE: (read_restart_i2asum) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 i2a data file.')
endif
end subroutine read_restart_i2asum

!=================================================
subroutine put_restart_i2a(fname, sec)
! call this subroutine after called get_restart_oi2
! it uses ocn_sst etc to calculate average ocn fields which will be used to send 
!to atm at the 1st step of continue run, because the ocn_sst cannot be sent to ice at the end of last run. 
! average ice fields (done at end of last run) are ready by calling read_restart_i2asum() 
!
implicit none

character*(*), intent(in) :: fname
integer :: sec

  if ( file_exist('i2a.nc') ) then
    write(il_out,*)' calling read_restart_i2a at time_sec = ',sec
    call read_restart_i2a('i2a.nc', sec)  
  endif
  if ( file_exist('i2asum.nc') ) then
    write(il_out,*)' calling read_restart_i2asum at time_sec = ',sec
    call read_restart_i2asum('i2asum.nc', sec) 

    write(il_out,*)' calling ave_ocn_fields_4_i2a at time_sec = ',sec
    call time_average_ocn_fields_4_i2a  !accumulate/average ocn fields needed for IA coupling
    write(il_out,*) ' calling get_i2a_fields at time_sec =', sec
    call get_i2a_fields
  endif

end subroutine put_restart_i2a

!=================================================
subroutine get_restart_o2i(fname)

! To be called at beginning of each run trunk to read in restart o2i fields

implicit none

character*(*), intent(in) :: fname
 
integer(kind=int_kind) :: ncid_o2i
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(get_restart_o2i) reading in o2i fields......'
  endif
  call ice_open_nc(fname, ncid_o2i)
  call ice_read_nc(ncid_o2i, 1, 'sst_i',    ocn_sst,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'sss_i',    ocn_sss,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssu_i',    ocn_ssu,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssv_i',    ocn_ssv,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'sslx_i',   ocn_sslx,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssly_i',   ocn_ssly,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'pfmice_i', ocn_pfmice, dbug)
  call ice_read_nc(ncid_o2i, 1, 'co2_oi',   ocn_co2, dbug)
  call ice_read_nc(ncid_o2i, 1, 'co2fx_oi',  ocn_co2fx, dbug)
  if (my_task == master_task) then
    call ice_close_nc(ncid_o2i)
    write(il_out,*) '(get_restart_o2i) has read in 7 o2i fields.'
  endif
else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (get_restart_o2i) not found file *** ',fname
  endif
  print *, 'CICE: (get_restart_o2i) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 o2i data file.')
endif

return
end subroutine get_restart_o2i

!=================================================
subroutine get_restart_mice(fname)

! Called at beginning of the run to get 'last' IO cpl int T-M ice variables 
! which are used together with the first received i2a fields to obtain the first 
! i2o fields sent to ocn immediately as the 1st io cpl int forcing there.

implicit none

character*(*), intent(in) :: fname

integer(kind=int_kind) :: ncid_o2i
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(get_restart_mice) opening file: ', fname
  endif

  call ice_open_nc(fname, ncid_o2i)
  call ice_read_nc(ncid_o2i, 1, 'maicen1',   maicen_saved(:,:,1,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maicen2',   maicen_saved(:,:,2,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maicen3',   maicen_saved(:,:,3,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maicen4',   maicen_saved(:,:,4,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maicen5',   maicen_saved(:,:,5,:), dbug)
  call ice_read_nc(ncid_o2i, 1, 'maice',     maice,     dbug)
  call ice_read_nc(ncid_o2i, 1, 'mstrocnxT', mstrocnxT, dbug)
  call ice_read_nc(ncid_o2i, 1, 'mstrocnyT', mstrocnyT, dbug)
  call ice_read_nc(ncid_o2i, 1, 'mfresh',    mfresh,    dbug)
  call ice_read_nc(ncid_o2i, 1, 'mfsalt',    mfsalt,    dbug)
  call ice_read_nc(ncid_o2i, 1, 'mfhocn',    mfhocn,    dbug)
  call ice_read_nc(ncid_o2i, 1, 'mfswthru',  mfswthru,  dbug)
  call ice_read_nc(ncid_o2i, 1, 'msicemass', msicemass, dbug)
  write(il_out,*) '(get_restart_mice) ALL variables read in! '

  if (my_task == master_task) then
    call ice_close_nc(ncid_o2i)
    write(il_out,*) '(get_restart_mice) has read in 8 T-M variables.'
  endif
else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (get_restart_mice) not found file *** ',fname
  endif
  print *, 'CICE: (get_restart_mice) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 mice data file.')
endif

return
end subroutine get_restart_mice


!=================================================
subroutine get_lice_discharge(fname)

! Called at beginning of each run trunk to read in land ice discharge mask or iceberg
! (off Antarctica and Greenland).

implicit none

character(len=*), intent(in) :: fname
character*80 :: myvar = 'ficeberg'
integer(kind=int_kind) :: ncid_i2o, im, k, i, j
logical :: dbug = .true. 

call ice_open_nc(trim(fname), ncid_i2o)

write(il_out,*) '(get_lice_discharge) opened datafile: ', trim(fname)
write(il_out,*) '(get_lice_discharge) ncid_i2o= ', ncid_i2o

if (iceberg .lt. 1 .or. iceberg .gt. 4) then
  write(il_out,*) '(get_lice_discharge) in ESM only supports iceberg = 1,2,3,4) '
  call abort_ice('CICE stopped: ESM only supports iceberg = 1,2,3,4. Please set it to 1,2,3,4')
else
  call gather_global(gtarea, tarea, master_task, distrb_info)
  select case (iceberg)
    case (1); myvar = 'FICEBERG_AC2'
    case (2); myvar = 'FICEBERG_GC3'
    case (3); myvar = 'FICEBERG_AC2_AVE'
    case (4); myvar = 'FICEBERG_GC3_AVE'
  end select
  write(il_out,*)'(get_lice_discharge), iceberg = ', iceberg
  write(il_out,'(a,a)') '(get_lice_discharge) reading in iceberg data, myvar= ',trim(myvar)
  do im = 1, 12
    write(il_out,*) '(get_lice_discharge) reading in data, month= ',im
    call ice_read_nc(ncid_i2o, im, trim(myvar), vwork, dbug)

    ! Restrict iceberg fluxes to ocean points
    where (tmask)
      vwork = vwork
    else where
      vwork = c0
    end where

    call gather_global(gwork, vwork, master_task, distrb_info)

    if ( my_task == master_task ) then 
      gicebergfw(:,:,im) = gwork(:,:)

      ticeberg_s(im) = 0.0
      do j = 1, iceberg_je_s  !1, ny_global/2 (iceberg_je_s smaller than ny_global/2 thus saves time)
        do i = 1, nx_global
          ticeberg_s(im) = ticeberg_s(im) + gtarea(i,j) * gwork(i,j)
        enddo
      enddo
      ticeberg_n(im) = 0.0
      do j = iceberg_js_n, ny_global  !ny_global/2 + 1, ny_global !(iceberg_js_n bigger than ny_global/2 +1)
        do i = 1, nx_global
          ticeberg_n(im) = ticeberg_n(im) + gtarea(i,j) * gwork(i,j)
        enddo
      enddo

      write(il_out, *) '(get_lice_discharge) check: im, ticeberg_s, ticeberg_n = ',im, ticeberg_s(im), ticeberg_n(im)
    endif

  enddo

endif
if (my_task == master_task) then
  call ice_close_nc(ncid_i2o)
endif

return

end subroutine get_lice_discharge


!=================================================
subroutine get_restart_i2o(fname)

! To be called at beginning of each run trunk to read in restart i2o fields

implicit none

character*(*), intent(in) :: fname

integer(kind=int_kind) :: ncid_i2o, jf, jfs
logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  if (my_task==0) then
    write(il_out,*) '(get_time0_i2o_fields) reading in i2o fields......'
  endif
  call ice_open_nc(fname, ncid_i2o)
  do jf = nsend_i2a + 1, jpfldout
    call ice_read_nc(ncid_i2o, 1, cl_writ(jf) , vwork, dbug)
    select case(trim(cl_writ(jf)))
    case ('strsu_io'); io_strsu = vwork
    case ('strsv_io'); io_strsv = vwork
    case ('rain_io');  io_rain = vwork
    case ('snow_io');  io_snow = vwork
    case ('stflx_io'); io_stflx = vwork
    case ('htflx_io'); io_htflx = vwork
    case ('swflx_io'); io_swflx = vwork
    case ('qflux_io'); io_qflux = vwork
    case ('shflx_io'); io_shflx = vwork
    case ('lwflx_io'); io_lwflx = vwork
    case ('runof_io'); io_runof = vwork
    case ('press_io'); io_press = vwork
    case ('aice_io');  io_aice  = vwork
    case ('melt_io');  io_melt  = vwork
    case ('form_io');  io_form  = vwork
    case ('co2_i1');  io_co2  = vwork
    case ('wnd_i1');  io_wnd  = vwork
!2 more added 20171024:
    case ('lice_fw');  io_licefw = vwork
    case ('lice_ht');  io_liceht = vwork
    end select
  enddo
  if (my_task == master_task) then
    call ice_close_nc(ncid_i2o)
    write(il_out,*) '(get_time0_i2o_fields) has read in 19 i2o fields.'
  endif
else
  if (my_task==0) then
    write(il_out,*) 'ERROR: (get_time0_i2o_fields) not found file *** ',fname
  endif
  print *, 'CICE: (get_time0_i2o_fields_old) not found file *** ',fname
  call abort_ice('CICE stopped -- Need time0 i2o data file.')
endif

return
end subroutine get_restart_i2o

!=================================================
subroutine set_sbc_ice
!-------------------------
!This routine is NOT used! 
!-------------------------
!
! Set coupling fields (in units of GMB, from UM and MOM4) needed for CICE
!
! Adapted from "subroutine cice_sbc_in" of HadGem3 Nemo "MODULE sbcice_cice"
! for the "nsbc = 5" case.
!
! It should be called after calling "from_atm" and "from_ocn". 
!-------------------------------------------------------------------------------

implicit none

real :: r1_S0
real, dimension(nx_block,ny_block,nblocks) :: zzs

integer :: i,j,k,cat

!*** Fields from UM (all on T cell center):

!(1) windstress taux:
strax = um_taux * maice      !*tmask   

!(2) windstress tauy:
stray = um_tauy * maice      !*tmask

!(3) surface downward latent heat flux (==> multi-category)
do j = 1, ny_block
do i = 1, nx_block
  do k = 1, nblocks
    if (maice(i,j,k)==0.0) then
      do cat = 1, ncat
        flatn_f(i,j,cat,k) = 0.0
      enddo
      ! This will then be conserved in CICE (done in sfcflux_to_ocn)
      flatn_f(i,j,1,k) = um_lhflx(i,j,k) 
    else
      do cat = 1, ncat
        flatn_f(i,j,cat,k) = um_lhflx(i,j,k) * maicen(i,j,cat,k)/maice(i,j,k)
        !???: flatn_f(i,j,cat,k) = um_iceevp(i,j,cat,k) * Lsub 
        !flatn_f(i,j,cat,k) = - um_iceevp(i,j,cat,k) * Lsub
      enddo
    endif
  enddo
enddo
enddo 

! GBM conductive flux through ice:
!(4-8) top melting; (9-13) bottom belting ==> surface heatflux
do cat = 1, ncat
  fcondtopn_f(:,:,cat,:) = um_bmlt(:,:,cat,:)
  fsurfn_f   (:,:,cat,:) = um_tmlt(:,:,cat,:) + um_bmlt(:,:,cat,:)
enddo

!(14) snowfall
fsnow = max(maice * um_snow, 0.0)

!(15) rainfall
frain = max(maice * um_rain, 0.0)

!BX-20160718: "save" the ice concentration "maice" used here for scaling-up frain etc in
!ice_step for "consistency"--
!maice_saved = maice 

!*** Fields from MOM4 (SSU/V and sslx/y are on U points): 

!(1) freezing/melting potential
frzmlt = ocn_pfmice
if (limit_icemelt) then
  frzmlt(:,:,:) = max(frzmlt(:,:,:), meltlimit)
endif

!(2) SST
!make sure SST is 'all right' K==>C
sst = ocn_sst
if (maxval(sst).gt.200) then
  sst = sst -273.15
endif

!(3) SSS
sss = ocn_sss

!(4) SSU
uocn = ocn_ssu

!(5) SSV
vocn = ocn_ssv

!(6) surface slope sslx
ss_tltx = ocn_sslx

!(7) surface slope ssly
ss_tlty = ocn_ssly

!(as per S.O.) make sure Tf is properly initialized
Tf (:,:,:) = -depressT*sss(:,:,:)  ! freezing temp (C)
!
!B: May use different formula for Tf such as TEOS-10 formulation: 
!
!r1_S0 = 0.875/35.16504
!zzs(:,:,:) = sqrt(abs(sss(:,:,:)) * r1_S0)
!Tf(:,:,:) = ((((1.46873e-03 * zzs(:,:,:) - 9.64972e-03) * zzs(:,:,:) + &
!               2.28348e-02) * zzs(:,:,:) - 3.12775e-02) * zzs(:,:,:) + &
!               2.07679e-02) * zzs(:,:,:) - 5.87701e-02
!Tf(:,:,:) = Tf(:,:,:) * sss(:,:,:) ! - 7.53e-4 * 5.0 !!!5.0 is depth in meters

end subroutine set_sbc_ice

!===============================================================================
subroutine get_sbc_ice
!
! ** Purpose: set GBM coupling fields (from UM and MOM4) needed for CICE
!
! Adapted from "subroutine cice_sbc_in" of HadGem3 Nemo "MODULE sbcice_cice"
!    for the "nsbc = 5" case.
!
! It should be called after calling "from_atm" and "from_ocn". 
!-------------------------------------------------------------------------------

implicit none

real :: r1_S0
real, dimension(nx_block,ny_block,nblocks) :: zzs

integer :: i,j,k,cat

! Fields from UM (all on T cell center):

!(1) windstress taux:
strax = um_taux * aice

!(2) windstress tauy:
stray = um_tauy * aice

!(3) surface downward latent heat flux (==> multi_category)
do j = 1, ny_block
do i = 1, nx_block
  do k = 1, nblocks
    ! Notes from ACCESS-CM2:
      !BX 20160826: as in NEMO sbccpl.F90, there is no "open water field" um_lhflx involved: 
      !    qla_ice(:,:,1:jpl) = - frcv(jpr_ievp)%z3(:,:,1:jpl) * lsub 
      !-------------------------------------------------------------------------------------
    ! CM2 Uses this, as um_lhflx is not available:
    !      do cat = 1, ncat
    !        flatn_f(i,j,cat,k) = - um_iceevp(i,j,cat,k) * Lsub
    !      enddo

    if (aice(i,j,k)==0.0) then
      do cat = 1, ncat
        flatn_f(i,j,cat,k) = 0.0
      enddo
      ! This will then be conserved in CICE (done in sfcflux_to_ocn)
      flatn_f(i,j,1,k) = um_lhflx(i,j,k)
    else
      do cat = 1, ncat
        flatn_f(i,j,cat,k) = um_lhflx(i,j,k) * aicen(i,j,cat,k)/aice(i,j,k)
      enddo
    endif
  enddo
enddo
enddo

! GBM conductive flux through ice:
!(4-8) top melting; (9-13) bottom belting ==> surface heatflux
do cat = 1, ncat
  fcondtopn_f(:,:,cat,:) = um_bmlt(:,:,cat,:)
  fsurfn_f   (:,:,cat,:) = um_tmlt(:,:,cat,:) + um_bmlt(:,:,cat,:)
enddo

fsnow = max(aice * um_snow,0.0)
frain = max(aice * um_rain,0.0)  

! Fields from MOM4 (SSU/V and sslx/y are on U points): 

!(1) freezing/melting potential
frzmlt = ocn_pfmice
!20080312: set maximum melting htflux allowed from ocn, (eg, -200 W/m^2)
!          the artificial "meltlimit = -200 " is read in from input_ice.nml
!20090320: set option 'limit_icemelt' in case no limit needed if cice behaves!
if (limit_icemelt) then
  frzmlt(:,:,:) = max(frzmlt(:,:,:), meltlimit)
endif

!(2) SST
sst = ocn_sst -273.15

!(3) SSS
sss = ocn_sss

!(4) SSU
uocn = ocn_ssu

!(5) SSV
vocn = ocn_ssv
!(6) surface slope sslx
                                                             
ss_tltx = ocn_sslx

!(7) surface slope ssly
ss_tlty = ocn_ssly

! * (as per S. O'Farrel) make sure Tf if properly initialized
!----- should use eos formula to calculate Tf for "consistency" with GCx ----!
Tf (:,:,:) = -depressT*sss(:,:,:)  ! freezing temp (C)
!
!B: May use different formula for Tf such as TEOS-10 formulation: 
!
!r1_S0 = 0.875/35.16504
!zzs(:,:,:) = sqrt(abs(sss(:,:,:)) * r1_S0)
!Tf(:,:,:) = ((((1.46873e-03 * zzs(:,:,:) - 9.64972e-03) * zzs(:,:,:) + &
!               2.28348e-02) * zzs(:,:,:) - 3.12775e-02) * zzs(:,:,:) + &
!               2.07679e-02) * zzs(:,:,:) - 5.87701e-02
!Tf(:,:,:) = Tf(:,:,:) * sss(:,:,:) ! - 7.53e-4 * 5.0 !!!5.0 is depth in meters
!
end subroutine get_sbc_ice

!=================================================
subroutine save_restart_o2i(fname, nstep)

! output the last o2i forcing data received in cice by the end of the run, 
! to be read in at the beginning of next run by cice

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid
integer(kind=int_kind) :: jf, jfs, ll, ilout

if (my_task == 0) then
  call create_ncfile(fname, ncid, il_im, il_jm, ll=1, ilout=il_out)
  call write_nc_1Dtime(real(nstep), 1, 'time', ncid)
endif

do jf = nrecv_a2i + 1, jpfldin

  select case (trim(cl_read(jf)))
  case('sst_i'); vwork = ocn_sst
  case('sss_i'); vwork = ocn_sss
  case('ssu_i'); vwork = ocn_ssu
  case('ssv_i'); vwork = ocn_ssv
  case('sslx_i'); vwork = ocn_sslx
  case('ssly_i'); vwork = ocn_ssly
  case('pfmice_i'); vwork = ocn_pfmice
  case('co2_oi'); vwork = ocn_co2
  case('co2fx_oi'); vwork = ocn_co2fx
  end select

  call gather_global(gwork, vwork, master_task, distrb_info)
  if (my_task == 0) then
    call write_nc2D(ncid, cl_read(jf), gwork, 2, il_im, il_jm, 1, ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck( nf_close(ncid) )

return
end subroutine save_restart_o2i

!=================================================
subroutine save_restart_i2asum(fname, nstep)
! output the last i2a forcing data in cice at the end of the run,
! to be read in at the beginning of next run by cice and sent to atm

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid
integer(kind=int_kind) :: jf, jfs, ll, ilout

integer(kind=int_kind), parameter :: sumfldin  = 46     !21
character(len=8), dimension(sumfldin) :: sumfld

sumfld(1)='msst'
sumfld(2)='mssu'
sumfld(3)='mssv'
sumfld(4)='muvel'
sumfld(5)='mvvel'
sumfld(6)='maiu'
sumfld(7)='maicen1'
sumfld(8)='maicen2'
sumfld(9)='maicen3'
sumfld(10)='maicen4'
sumfld(11)='maicen5'
sumfld(12)='mthikn1'
sumfld(13)='mthikn2'
sumfld(14)='mthikn3'
sumfld(15)='mthikn4'
sumfld(16)='mthikn5'
sumfld(17)='msnown1'
sumfld(18)='msnown2'
sumfld(19)='msnown3'
sumfld(20)='msnown4'
sumfld(21)='msnown5'

if (my_task == 0) then
  call create_ncfile(fname, ncid, il_im, il_jm, ll=1, ilout=il_out)
endif

do jf = 1, sumfldin
    select case (trim(sumfld(jf)))
    case('msst'); vwork = msst
    case('mssu'); vwork = mssu
    case('mssv'); vwork = mssv
    case('muvel'); vwork = muvel
    case('mvvel'); vwork = mvvel
    case('maiu'); vwork = maiu
    case('maicen1'); vwork = maicen(:,:,1,:)
    case('maicen2'); vwork = maicen(:,:,2,:)
    case('maicen3'); vwork = maicen(:,:,3,:)
    case('maicen4'); vwork = maicen(:,:,4,:)
    case('maicen5'); vwork = maicen(:,:,5,:)
    case('mthikn1'); vwork = mthikn(:,:,1,:)
    case('mthikn2'); vwork = mthikn(:,:,2,:)
    case('mthikn3'); vwork = mthikn(:,:,3,:)
    case('mthikn4'); vwork = mthikn(:,:,4,:)
    case('mthikn5'); vwork = mthikn(:,:,5,:)
    case('msnown1'); vwork = msnown(:,:,1,:)
    case('msnown2'); vwork = msnown(:,:,2,:)
    case('msnown3'); vwork = msnown(:,:,3,:)
    case('msnown4'); vwork = msnown(:,:,4,:)
    case('msnown5'); vwork = msnown(:,:,5,:)
    end select
    call gather_global(gwork, vwork, master_task, distrb_info)
  if (my_task == 0) then
    call write_nc2D(ncid, sumfld(jf), gwork, 2, il_im, il_jm, 1, ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck( nf_close(ncid) )

end subroutine save_restart_i2asum

!=================================================
subroutine save_restart_mice(fname, nstep)

! output ice variable averaged over the last IO cpl int of this run, 
! cice reads in these vars at the beginning of next run, uses them with the first 
! received a2i fields to obtain the first i2o fields to be sent to ocn  

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid
integer(kind=int_kind) :: jf, jfs, ll, ilout

if (my_task == 0) then
  call create_ncfile(fname, ncid, il_im, il_jm, ll=1, ilout=il_out)
  call write_nc_1Dtime(real(nstep), 1, 'time', ncid)
endif

! maicen_saved appears to be the same as maicen in CICE5-UM7.3 
!B: 20170825 ==> add maicen_saved for atm_icefluxes_back2GBM calculation!
!        note maicen_saved is the last ia interval mean.  
vwork(:,:,:) = maicen_saved(:,:,1,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen1', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork(:,:,:) = maicen_saved(:,:,2,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen2', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork(:,:,:) = maicen_saved(:,:,3,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen3', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork(:,:,:) = maicen_saved(:,:,4,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen4', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork(:,:,:) = maicen_saved(:,:,5,:)
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maicen5', gwork, 2, il_im, il_jm, 1, ilout=il_out)
!b.

!The following fields are actually the ice state of last timestep 
!(no time-averaging is required in each timestep io coupling, see time_average_fields_4_i2o)
vwork = maice
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maice', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mstrocnxT
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mstrocnxT', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mstrocnyT
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mstrocnyT', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mfresh
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mfresh', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mfsalt
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mfsalt', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mfhocn
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mfhocn', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = mfswthru
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mfswthru', gwork, 2, il_im, il_jm, 1, ilout=il_out)
vwork = msicemass
call gather_global(gwork, vwork, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'msicemass', gwork, 2, il_im, il_jm, 1, ilout=il_out)

if (my_task == 0) call ncheck( nf_close(ncid) )

return
end subroutine save_restart_mice

!=================================================
subroutine get_i2a_fields

! all fields (except for vector) obtained here are all on T cell center

!(1) ocean surface temperature 
ia_sst(:,:,:) = msst(:,:,:)

!(2-3) ice/ocn combined surface velocity 
ia_uvel(:,:,:) = mssu(:,:,:) * (1. - maiu(:,:,:)) + muvel(:,:,:) * maiu(:,:,:) 
ia_vvel(:,:,:) = mssv(:,:,:) * (1. - maiu(:,:,:)) + mvvel(:,:,:) * maiu(:,:,:)

!(4-8) ice concentration
ia_aicen(:,:,:,:) = maicen(:,:,:,:)
!BX: save it for use in atm_icefluxes_back2GBM ---
maicen_saved = maicen

!(9-13) ice thickness
ia_thikn(:,:,:,:) = mthikn(:,:,:,:)

!(14-18) snow thickness
ia_snown(:,:,:,:) = msnown(:,:,:,:)

!(19-20) co2 flux stuff
ia_co2 = mco2
ia_co2fx = mco2fx

return
end subroutine get_i2a_fields

!=================================================
subroutine get_i2o_fields

! All fluxes should be in GBM units before passing into coupler.
! e.g.,  io_htflx(:,:,:) = fhocn_gbm(:,:,:)	

implicit none
real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: pice
integer (kind=int_kind) :: i,j,k

real (kind=dbl_kind) :: &
        trunoff_s  = c0, &
        trunoff_n  = c0, &
        r_s = c1, &
        r_n = c1, &
        r_runoff= c1        !=(1-min(r_max_iceberg. r_s(or r_n))

! Fields obtained here are all at T cell center. before being sent to MOM4, vector 
! (Taux, Tauy) should be shifted on to U point as required
!------------------------------------------------------------------------------- 

!(1-2) air/ice - sea stress TAUX/TAUY
!    Caution: in nemo, "strocnx/y" are NOT weighted by aice here, 'cos strocnx/y
!    have already been 'weighted' using aice (when calculated in "evp_finish". 
!    But, this weight has been removed in strocnx/yT (see "evp_finish"), therfore 
!    we need put it on again here. 
io_strsu = um_taux * (c1 - maice) - mstrocnxT * maice
io_strsv = um_tauy * (c1 - maice) - mstrocnyT * maice

!(3) freshwater flux to ocean: rainfall (+ ice melting water flux ?)
io_rain = um_rain * (c1 - maice)
!202412: fixing watermass loss from ocean by adding a small, constant fwflux into rain--
io_rain = io_rain + add_lprec

!(4) freshwater flux to ocean: snowfall
io_snow = um_snow * (c1 - maice)

!(5) salt flux to ocean
io_stflx = mfsalt

!(6) ice/snow melting heatflux into ocean 
io_htflx = mfhocn

!(7) short wave radiation 
!(CH: the (1-aice) weight should not be here 'cos all fluxes passed in from
!     UM have already been aice-weighted when they are calculated there!!!) 
!io_swflx = um_swflx * (c1 - maice) + mfswthru
io_swflx = um_swflx + mfswthru

!(8) latent heat flux (positive out of ocean as required by MOM4)
io_qflux = um_evap * Lvap        !Note it's already weighted in UM for open sea.

!(9) sensible heat flux (positive out of ocean as required by MOM4)
io_shflx = um_shflx

!(10) net long wave radiation positive down
io_lwflx = um_lwflx

!(11) runoff 
!*** mask off "extra/useless" runoff on dry points ***
where (tmask)
  io_runof = um_runoff
else where
  um_runoff = c0
end where

call gather_global(grunoff, io_runof, master_task, distrb_info)

if (my_task == master_task) then

  trunoff_s = c0
  do j = 1, runoff_je_s
    do i = 1, nx_global
      trunoff_s = trunoff_s + gtarea(i,j) * grunoff(i,j)
      grunoff(i,j) = grunoff(i,j) * (c1 - iceberg_rate_s)  !do deduction
    enddo
  enddo
  trunoff_n = c0
  do j = runoff_js_n, runoff_je_n
    do i = runoff_is_n, runoff_ie_n
      trunoff_n = trunoff_n + gtarea(i,j) * grunoff(i,j)
      grunoff(i,j) = grunoff(i,j) * (c1 - iceberg_rate_n)  !do deduction
    enddo
  enddo
  !Now global runoff has been "updated" (deduction done for iceberg).
endif

!distributed the resultant runoff and iceberg fluxes: 
call scatter_global(vwork, grunoff, master_task, distrb_info, &
                                  field_loc_center, field_type_scalar)
io_runof(:,:,:) = vwork(:,:,:)

!2 new flux items associated with the iceberg discharged into ocean
!(18) water flux due to land ice melt off Antarctica and Greenland (kg/m^2/s)
!(19) heat  flux due to land ice melt off Antarctica and Greenland

!XXXXXX 
IF (my_task == master_task) THEN

  gwork(:,:) = c0
  do i = 1, nx_global
    do j = 1, iceberg_je_s
      gwork(i, j) = gicebergfw(i, j, month) * iceberg_rate_s * trunoff_s / ticeberg_s(month)
    enddo
    do j = iceberg_js_n, ny_global
      gwork(i, j) = gicebergfw(i, j, month) * iceberg_rate_n * trunoff_n / ticeberg_n(month)
    enddo
  enddo
  !Now global iceberg has been defined (using the deduction from runoff)

ENDIF
!XXXXXX

call scatter_global(vwork, gwork, master_task, distrb_info, &
                    field_loc_center, field_type_scalar)
io_licefw(:,:,:) = vwork(:,:,:)         !i2o field No 18.

!Also count in the latent heat carried with the runoff part, as done below, thus allowing
!for (rough) consistency of energy exchange no matter what iceberg_rate_s/n are used.
!Warning: the follow approach would lose all the runoff LH, if runoff_lh=.false., in no-iceberg case 

if (my_task == master_task) then
  if ( runoff_lh ) then

    do i = 1, nx_global
      do j = 1, runoff_je_s
        gwork(i,j) = gwork(i,j) + grunoff(i,j)
      enddo
    enddo
    do i = runoff_is_n, runoff_ie_n
      do j = runoff_js_n, runoff_je_n
        gwork(i,j) = gwork(i,j) + grunoff(i,j)
      enddo
    enddo
  else
  !If runoff with latent heat flux crashes the model in no-iceberg case, due probably to too big LH(?)
  !all the LH carried by runoff is applied to iceberg areas.
    do i = 1, nx_global
      do j = 1, iceberg_je_s
        gwork(i,j) = gwork(i,j)/max(iceberg_rate_s, 0.0001_dbl_kind) !get the whole runoff LH onto iceberg
      enddo
      do j = iceberg_js_n, ny_global
        gwork(i,j) = gwork(i,j)/max(iceberg_rate_n, 0.0001_dbl_kind)
      enddo
    enddo
  endif
endif

call scatter_global(vwork, gwork, master_task, distrb_info, &
                    field_loc_center, field_type_scalar)
io_liceht = - vwork * Lfresh * iceberg_lh       !FW converted into LH flux (W/m^2).
                                                !!i2o field No 19.
!(12) pressure
pice = gravit * msicemass
!----------------------------------------------------------------------------
!sicemass = rho_ice x hi + rho_snow x hs (in m)
if (ice_pressure_on) then
  io_press = pice * maice
endif
if (air_pressure_on) then
   !as GFDL SIS, we use patm anormaly, i.e., taking off 1.e5 Pa !
   io_press(:,:,:) = io_press(:,:,:) + um_press(:,:,:) - 1.0e5_dbl_kind 
endif
!(13) ice concentration
io_aice = maice
!(14) ice melt fwflux 
io_melt = max(c0,mfresh(:,:,:)) 
!(15) ice form fwflux
io_form = min(c0,mfresh(:,:,:))

!(16) CO2
io_co2 = um_co2
!(17) 10m winnspeed
io_wnd = um_wnd

return
end subroutine get_i2o_fields

!=================================================
subroutine initialize_mice_fields_4_i2o

implicit none

maice = 0.
mstrocnxT = 0.
mstrocnyT = 0.
mfresh = 0.
mfsalt = 0.
mfhocn = 0. 
mfswthru = 0.
msicemass = 0.

return
end subroutine initialize_mice_fields_4_i2o

!=================================================
subroutine initialize_mice_fields_4_i2a

implicit none

muvel = 0.
mvvel = 0.

maiu = 0.
maicen = 0.
mthikn = 0.
msnown = 0. 

return
end subroutine initialize_mice_fields_4_i2a

!=================================================
subroutine initialize_mocn_fields_4_i2a

implicit none

msst = 0.
mssu = 0.
mssv = 0.
mco2 = 0.
mco2fx = 0.

return
end subroutine initialize_mocn_fields_4_i2a

!=================================================
subroutine time_average_ocn_fields_4_i2a

implicit none

msst(:,:,:) = msst(:,:,:) + ocn_sst(:,:,:) * coef_ai
mssu(:,:,:) = mssu(:,:,:) + ocn_ssu(:,:,:) * coef_ai
mssv(:,:,:) = mssv(:,:,:) + ocn_ssv(:,:,:) * coef_ai
mco2(:,:,:) = mco2(:,:,:) + ocn_co2(:,:,:) * coef_ai
mco2fx(:,:,:) = mco2fx(:,:,:) + ocn_co2fx(:,:,:) * coef_ai

return
end subroutine time_average_ocn_fields_4_i2a

!=================================================
subroutine time_average_fields_4_i2o
!now for each timestep io coupling, so no time-averaging is required. 
implicit none

maice(:,:,:)     = aice(:,:,:)
mstrocnxT(:,:,:) = strocnxT(:,:,:) 
mstrocnyT(:,:,:) = strocnyT(:,:,:) 
mfresh(:,:,:)    = fresh(:,:,:) 
mfsalt(:,:,:)    = fsalt(:,:,:) 
mfhocn(:,:,:)    = fhocn(:,:,:) 
mfswthru(:,:,:)  = fswthru(:,:,:)
msicemass(:,:,:) = sicemass(:,:,:)

return
end subroutine time_average_fields_4_i2o

!=================================================
subroutine time_average_fields_4_i2a

implicit none

! ice fields:
muvel(:,:,:) = muvel(:,:,:) + uvel(:,:,:) * coef_ai
mvvel(:,:,:) = mvvel(:,:,:) + vvel(:,:,:) * coef_ai
maicen(:,:,:,:) = maicen(:,:,:,:) + aicen(:,:,:,:) * coef_ai  !T cat. ice concentration
mthikn(:,:,:,:) = mthikn(:,:,:,:) + vicen(:,:,:,:) * coef_ai  !T cat. ice thickness
msnown(:,:,:,:) = msnown(:,:,:,:) + vsnon(:,:,:,:) * coef_ai  !T cat. snow thickness

call to_ugrid(aice, aiiu)
maiu(:,:,:)  = maiu(:,:,:) + aiiu(:,:,:) * coef_ai            !U cell ice concentraction

!ocn fields:
!must be done after calling from_ocn so as to get the most recently updated ocn fields,
!therefore a separate call to "time_average_ocn_fields_4_i2a" is done for this purpose.

return
end subroutine time_average_fields_4_i2a

!=================================================
subroutine check_i2a_fields(nstep)

implicit none

integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ilout, ll, jf
integer(kind=int_kind), save :: ncid,currstep 
data currstep/0/ 

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist('fields_i2a_in_ice.nc') ) then
  call create_ncfile('fields_i2a_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening file fields_i2a_in_ice.nc at nstep = ', nstep
  call ncheck( nf_open('fields_i2a_in_ice.nc',nf_write,ncid) )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid) 
end if

do jf = 1, nsend_i2a
   
  select case(trim(cl_writ(jf)))
    case('isst_ia');  vwork = ia_sst
    case('icecon01'); vwork(:,:,:) = ia_aicen(:,:,1,:) 
    case('icecon02'); vwork(:,:,:) = ia_aicen(:,:,2,:) 
    case('icecon03'); vwork(:,:,:) = ia_aicen(:,:,3,:) 
    case('icecon04'); vwork(:,:,:) = ia_aicen(:,:,4,:) 
    case('icecon05'); vwork(:,:,:) = ia_aicen(:,:,5,:) 
    case('snwthk01'); vwork(:,:,:) = ia_snown(:,:,1,:)
    case('snwthk02'); vwork(:,:,:) = ia_snown(:,:,2,:)  
    case('snwthk03'); vwork(:,:,:) = ia_snown(:,:,3,:)
    case('snwthk04'); vwork(:,:,:) = ia_snown(:,:,4,:)
    case('snwthk05'); vwork(:,:,:) = ia_snown(:,:,5,:)
    case('icethk01'); vwork(:,:,:) = ia_thikn(:,:,1,:)
    case('icethk02'); vwork(:,:,:) = ia_thikn(:,:,2,:)
    case('icethk03'); vwork(:,:,:) = ia_thikn(:,:,3,:)
    case('icethk04'); vwork(:,:,:) = ia_thikn(:,:,5,:)
    case('icethk05'); vwork(:,:,:) = ia_thikn(:,:,5,:)
    case('uvel_ia');  vwork = ia_uvel 
    case('vvel_ia');  vwork = ia_vvel 
    case('co2_i2');  vwork = ia_co2
    case('co2fx_i2');  vwork = ia_co2fx
  end select

  call gather_global(gwork, vwork, master_task, distrb_info)
    
  if (my_task == 0 ) then
    call write_nc2D(ncid, trim(cl_writ(jf)), gwork, 1, il_im,il_jm,currstep,ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_i2a_fields

!=================================================
subroutine check_a2i_fields(nstep)

implicit none

integer(kind=int_kind), intent(in) :: nstep
character*80 :: ncfile='fields_a2i_in_ice_2.nc'
integer(kind=int_kind) :: ncid, currstep, ll, ilout, jf
data currstep/0/
save currstep

currstep=currstep+1

if ( my_task == 0 .and. .not. file_exist(trim(ncfile)) ) then
  call create_ncfile(trim(ncfile),ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening file ', trim(ncfile), ' at nstep = ', nstep
  call ncheck( nf_open(trim(ncfile),nf_write,ncid) )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

do jf = 1, nrecv_a2i

  select case (trim(cl_read(jf)))
    case ('thflx_i');  vwork = um_thflx
    case ('pswflx_i'); vwork = um_pswflx
    case ('runoff_i'); vwork = um_runoff
    case ('wme_i');    vwork = um_wme
    case ('rain_i');   vwork = um_rain
    case ('snow_i');   vwork = um_snow
    case ('evap_i');   vwork = um_evap
    case ('lhflx_i');  vwork = um_lhflx
    case ('tmlt01'); vwork(:,:,:) = um_tmlt(:,:,1,:)
    case ('tmlt02'); vwork(:,:,:) = um_tmlt(:,:,2,:)
    case ('tmlt03'); vwork(:,:,:) = um_tmlt(:,:,3,:)
    case ('tmlt04'); vwork(:,:,:) = um_tmlt(:,:,4,:)
    case ('tmlt05'); vwork(:,:,:) = um_tmlt(:,:,5,:)
    case ('bmlt01'); vwork(:,:,:) = um_tmlt(:,:,1,:)
    case ('bmlt02'); vwork(:,:,:) = um_tmlt(:,:,2,:)
    case ('bmlt03'); vwork(:,:,:) = um_tmlt(:,:,3,:)
    case ('bmlt04'); vwork(:,:,:) = um_tmlt(:,:,4,:)
    case ('bmlt05'); vwork(:,:,:) = um_tmlt(:,:,5,:)
    case ('taux_i'); vwork = um_taux
    case ('tauy_i'); vwork = um_tauy
    case ('swflx_i'); vwork = um_swflx
    case ('lwflx_i'); vwork = um_lwflx
    case ('shflx_i'); vwork = um_shflx
    case ('press_i'); vwork = um_press
    case ('co2_ai'); vwork = um_co2
    case ('wnd_ai'); vwork = um_wnd
  end select 

  call gather_global(gwork, vwork, master_task, distrb_info)

  if (my_task == 0) then
    call write_nc2D(ncid, trim(cl_read(jf)), gwork, 1, il_im,il_jm,currstep,ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_a2i_fields

!=================================================
subroutine check_i2o_fields(nstep, scale)

implicit none

integer(kind=int_kind), intent(in) :: nstep
real, intent(in) :: scale
integer(kind=int_kind) :: ncid, currstep, ll, ilout, jf
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist('fields_i2o_in_ice.nc') ) then
  call create_ncfile('fields_i2o_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening file fields_i2o_in_ice.nc at nstep = ', nstep
  call ncheck( nf_open('fields_i2o_in_ice.nc',nf_write,ncid) )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

do jf = nsend_i2a + 1, jpfldout

  select case(trim(cl_writ(jf)))
    case('strsu_io')
      vwork = scale * io_strsu
    case('strsv_io')
      vwork = scale * io_strsv
    case('rain_io')
      vwork = scale * io_rain
    case('snow_io')
      vwork = scale * io_snow
    case('stflx_io')
      vwork = scale * io_stflx
    case('htflx_io')
      vwork = scale * io_htflx
    case('swflx_io')
      vwork = scale * io_swflx
    case('qflux_io')
      vwork = scale * io_qflux
    case('shflx_io')
      vwork = scale * io_shflx
    case('lwflx_io')
      vwork = scale * io_lwflx
    case('runof_io')
      vwork = scale * io_runof
    case('press_io')
      vwork = scale * io_press
    case('aice_io')
      vwork = scale * io_aice
    case('form_io')
      vwork = scale * io_form
    case('melt_io')
      vwork = scale * io_melt
    case('co2_i1')
      vwork = scale * io_co2
    case('wnd_i1')
      vwork = scale * io_wnd
    !202407: 2 more fields added:
    case('lice_fw')
      vwork = scale * io_licefw
    case('lice_ht')
      vwork = scale * io_liceht      
  end select

  call gather_global(gwork, vwork, master_task, distrb_info)

  if (my_task == 0 ) then
    call write_nc2D(ncid, trim(cl_writ(jf)), gwork, 1, il_im,il_jm,currstep,ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_i2o_fields

!=================================================
subroutine check_o2i_fields(nstep)

implicit none

integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid, currstep, ilout, ll, jf
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist('fields_o2i_in_ice.nc') ) then
  call create_ncfile('fields_o2i_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening file fields_o2i_in_ice.nc at nstep = ', nstep
  call ncheck( nf_open('fields_o2i_in_ice.nc',nf_write,ncid) )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

do jf = nrecv_a2i + 1, jpfldin

  select case (trim(cl_read(jf)))
    case ('sst_i'); vwork = ocn_sst
    case ('sss_i'); vwork = ocn_sss
    case ('ssu_i'); vwork = ocn_ssu
    case ('ssv_i'); vwork = ocn_ssv
    case ('sslx_i'); vwork = ocn_sslx
    case ('ssly_i'); vwork = ocn_ssly
    case ('pfmice_i'); vwork = ocn_pfmice
    case ('co2_oi'); vwork = ocn_co2
    case ('co2fx_oi'); vwork = ocn_co2fx
  end select

  call gather_global(gwork, vwork, master_task, distrb_info)

  if (my_task == 0) then
    call write_nc2D(ncid, trim(cl_read(jf)), gwork, 1, il_im,il_jm,currstep,ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_o2i_fields

!============================================================================
subroutine check_frzmlt_sst(ncfilenm)

!this is (mainly) used to check cice solo run frzmlt and sst !
! (for comparison against a coupled run forcing into cice)

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, sst, master_task, distrb_info) 
if (my_task == 0) call write_nc2D(ncid, 'sst', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, frzmlt, master_task, distrb_info) 
if (my_task == 0) call write_nc2D(ncid, 'frzmlt', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_frzmlt_sst

!=================================================
subroutine check_i2o_uvfluxes(ncfilenm)

!this is temporarily used to check i2o fields (uflux, vflux and maice) 
!for debug purpose

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, maice, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'maice', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, um_taux, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'um_taux', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, um_tauy, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'um_tauy', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, mstrocnxT, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mstrocnxT', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, mstrocnyT, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'mstrocnyT', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_i2o_uvfluxes

!=================================================
subroutine check_ice_fields(ncfilenm)
!this is temporarily used to check ice fields immediately after ice_step 
!for debug purpose

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, aice, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aice', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,1,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,2,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,3,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,4,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, aicen(:,:,5,:), master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aicen5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_ice_fields

!=================================================
subroutine check_ice_sbc_fields(ncfilenm)

!this is temporarily used to check ice_sbc fields got from get_ice_sbc 
!for debug purpose

implicit none

real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: v3d

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

v3d = 0.0

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, aice, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'aice', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, strax, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'strax', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, stray, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'stray', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

v3d(:,:,:) = flatn_f(:,:,1,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = flatn_f(:,:,2,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = flatn_f(:,:,3,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = flatn_f(:,:,4,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = flatn_f(:,:,5,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'flatn_f5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

v3d(:,:,:) = fcondtopn_f(:,:,1,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fcondtopn_f(:,:,2,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fcondtopn_f(:,:,3,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fcondtopn_f(:,:,4,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fcondtopn_f(:,:,5,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fcondtopn_f5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

v3d(:,:,:) = fsurfn_f(:,:,1,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f1', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fsurfn_f(:,:,2,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f2', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fsurfn_f(:,:,3,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f3', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fsurfn_f(:,:,4,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f4', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
v3d(:,:,:) = fsurfn_f(:,:,5,:)
call gather_global(gwork, v3d, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsurfn_f5', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

call gather_global(gwork, fsnow, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'fsnow', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, frain, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'frain', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

!from ocen:
call gather_global(gwork, frzmlt, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'frzmlt', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, sst, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sst', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, sss, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sss', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, uocn, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'uocn', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, vocn, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'vocn', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

call gather_global(gwork, ss_tltx, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ss_tltx', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, ss_tlty, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ss_tlty', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

call gather_global(gwork, Tf, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'Tf', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_ice_sbc_fields

!=================================================
subroutine check_sstsss(ncfilenm)

!this is used to check cice sst/sss : temporary use (20091019)

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, sst, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sst', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, sss, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sss', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return
end subroutine check_sstsss

!=================================================
subroutine check_iceberg_fields(ncfilenm)

!this is used to check land ice fields

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, io_licefw, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'io_licefw', gwork, 1, il_im,il_jm,currstep,ilout=il_out)
call gather_global(gwork, io_liceht, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'io_liceht', gwork, 1, il_im,il_jm,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid))

return

end subroutine check_iceberg_fields

!=================================================
subroutine check_landice_fields_1(ncfilenm)

!this is used to check land ice fields

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

if (my_task == 0) call ncheck(nf_close(ncid))

return

end subroutine check_landice_fields_1

!=================================================
subroutine check_landice_fields_2(ncfilenm)

!this is used to check land ice fields

implicit none

character*(*), intent(in) :: ncfilenm
integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(ncfilenm) ) then
  call create_ncfile(ncfilenm,ncid,il_im,il_jm,ll=1,ilout=il_out)
endif

if (my_task == 0) then
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
  call ncheck( nf_open(ncfilenm, nf_write,ncid) )
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

if (my_task == 0) call ncheck(nf_close(ncid))

return

end subroutine check_landice_fields_2

!=================================================
function file_exist (file_name)
!
character(len=*), intent(in) :: file_name
logical  file_exist

file_exist = .false.
if (len_trim(file_name) == 0) return
if (file_name(1:1) == ' ')    return

inquire (file=trim(file_name), exist=file_exist)

end function file_exist

!=================================================

end module cpl_forcing_handler
