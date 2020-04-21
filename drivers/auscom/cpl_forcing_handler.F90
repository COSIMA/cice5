module cpl_forcing_handler
!
!==============================================================================
!
    use ice_blocks
    use ice_forcing
    use ice_read_write
    use ice_domain_size
    use ice_domain,    only : distrb_info, nblocks 
    use ice_flux                           !forcing data definition (Tair, Qa, uocn, etc.)
    use ice_state,     only : aice, trcr   !ice concentration and tracers
    use ice_gather_scatter           !gather and scatter routines...
!        use ice_constants, only : gravit, Tocnfrz
    use ice_constants
    use ice_grid,      only : tmask
    use ice_communicate, only : my_task, master_task
    use ice_ocean,     only : cprho
!20091118
    use ice_shortwave, only : ocn_albedo2D
    use ice_exit, only: abort_ice

    use cpl_parameters
    use cpl_netcdf_setup
    use cpl_arrays_setup
    use ice_calendar, only: dt

implicit none

contains

!===============================================================================
subroutine nullify_i2o_fluxes()

    iostrsu(:,:,:) = 0.0
    iostrsv(:,:,:) = 0.0
    iorain(:,:,:) = 0.0
    iosnow(:,:,:) = 0.0
    iostflx(:,:,:) = 0.0
    iohtflx(:,:,:) = 0.0
    ioswflx(:,:,:) = 0.0
    ioqflux(:,:,:) = 0.0
    ioshflx(:,:,:) = 0.0
    iolwflx(:,:,:) = 0.0
    iopress(:,:,:) = 0.0
    iorunof(:,:,:) = 0.0
    ioaice(:,:,:)  = 0.0
    iomelt(:,:,:)  = 0.0
    ioform(:,:,:)  = 0.0
    iolicefw(:,:,:)  = 0.0
    iolicefh(:,:,:)  = 0.0

end subroutine nullify_i2o_fluxes

!===============================================================================
subroutine tavg_i2o_fluxes

    iostrsu(:,:,:) = iostrsu(:,:,:) + tiostrsu(:,:,:)*coef_ic
    iostrsv(:,:,:) = iostrsv(:,:,:) + tiostrsv(:,:,:)*coef_ic
    iorain (:,:,:) = iorain (:,:,:) + tiorain (:,:,:)*coef_ic
    iosnow (:,:,:) = iosnow (:,:,:) + tiosnow (:,:,:)*coef_ic
    iostflx(:,:,:) = iostflx(:,:,:) + tiostflx(:,:,:)*coef_ic
    iohtflx(:,:,:) = iohtflx(:,:,:) + tiohtflx(:,:,:)*coef_ic
    ioswflx(:,:,:) = ioswflx(:,:,:) + tioswflx(:,:,:)*coef_ic
    ioqflux(:,:,:) = ioqflux(:,:,:) + tioqflux(:,:,:)*coef_ic
    ioshflx(:,:,:) = ioshflx(:,:,:) + tioshflx(:,:,:)*coef_ic
    iolwflx(:,:,:) = iolwflx(:,:,:) + tiolwflx(:,:,:)*coef_ic
    iorunof(:,:,:) = iorunof(:,:,:) + tiorunof(:,:,:)*coef_ic
    iopress(:,:,:) = iopress(:,:,:) + tiopress(:,:,:)*coef_ic
    ioaice(:,:,:)  = ioaice(:,:,:)  + tioaice(:,:,:)*coef_ic
!!!
    iomelt (:,:,:) = iomelt (:,:,:) + tiomelt (:,:,:)*coef_ic
    ioform (:,:,:) = ioform (:,:,:) + tioform (:,:,:)*coef_ic
    iolicefw (:,:,:) = iolicefw (:,:,:) + tiolicefw (:,:,:)*coef_ic
    iolicefh (:,:,:) = iolicefh (:,:,:) + tiolicefh (:,:,:)*coef_ic

return
end subroutine tavg_i2o_fluxes

!===============================================================================
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
  if (my_task == 0) print *, '(get_core_runoff) file doesnt exist: ', fname
  stop 'CICE stopped: core runoff (remapped) file not found.'
endif

return
end subroutine get_core_runoff

!===============================================================================
subroutine get_time0_sstsss(fname, nmonth)

implicit none

character*(*), intent(in) :: fname
    integer(kind=int_kind), intent(in) :: nmonth
    logical :: dbug
    integer(kind=int_kind) :: ncid

dbug = .true.

if ( file_exist(fname) ) then
  call ice_open_nc(fname, ncid)
  call ice_read_global_nc(ncid, nmonth, 'TEMP', gwork, dbug)
  call scatter_global(sst, gwork, master_task, distrb_info, &
                      field_loc_center, field_type_scalar)
  call ice_read_global_nc(ncid, nmonth, 'SALT', gwork, dbug)
  call scatter_global(sss, gwork, master_task, distrb_info, &
                      field_loc_center, field_type_scalar)
  if (my_task == master_task) call ice_close_nc(ncid)
else
  if (my_task == 0) print *, '(get_time0_sstsss) file doesnt exist: ', fname
  stop 'CICE stopped--initial SST and SSS ncfile not found.'
endif

return
end subroutine get_time0_sstsss

!===============================================================================
subroutine get_time0_sstsss_OLD(fname, nmonth)

implicit none

character*(*), intent(in) :: fname
    integer(kind=int_kind), intent(in) :: nmonth
    logical :: dbug
    integer(kind=int_kind) :: ncid

dbug = .true.
!dbug = .false.
if ( file_exist(fname) ) then
  call ice_open_nc(fname, ncid)
  call ice_read_nc(ncid, nmonth, 'TEMP', sst, dbug)
  call gather_global(gwork, sst, master_task, distrb_info)
  call ice_read_nc(ncid, nmonth, 'SALT', sss, dbug)
  call gather_global(gwork, sss, master_task, distrb_info)
  if (my_task == master_task) call ice_close_nc(ncid)
else
  if (my_task == 0) print*, '(get_time0_sstsss) file doesnt exist: ', fname
  stop 'CICE stopped--initial SST and SSS ncfile not found.'
endif

return
end subroutine get_time0_sstsss_OLD

!===============================================================================
subroutine get_time0_o2i_fields(fname)

implicit none

character*(*), intent(in) :: fname
 
    integer(kind=int_kind) :: ncid_o2i
    logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  call ice_open_nc(fname, ncid_o2i)
  call ice_read_nc(ncid_o2i, 1, 'sst_i',    ssto,   dbug, field_loc_center, field_type_scalar)
  call ice_read_nc(ncid_o2i, 1, 'sss_i',    ssso,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssu_i',    ssuo,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssv_i',    ssvo,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'sslx_i',   sslx,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'ssly_i',   ssly,   dbug)
  call ice_read_nc(ncid_o2i, 1, 'pfmice_i', pfmice, dbug)
  if (my_task == master_task) call ice_close_nc(ncid_o2i)
else
  print *, 'CICE: (get_time0_o2i_fields_old) not found file *** ',fname
  stop 'CICE stopped -- Need time0 o2i data file.'
endif

return
end subroutine get_time0_o2i_fields

!===============================================================================
subroutine get_time0_i2o_fields(fname)

implicit none

    character*(*), intent(in) :: fname

    integer(kind=int_kind) :: ncid_i2o, i

    if ( file_exist(fname) ) then
#if defined(DEBUG)
        if (my_task == master_task) write(il_out,*) '(get_time0_i2o_fields) reading in i2o fields......'
#endif
        call ice_open_nc(fname, ncid_i2o) 
        do i=1, num_fields_to_ocn
            vwork(:, :, :) = 0.0
            call ice_read_nc(ncid_i2o, 1, trim(fields_to_ocn(i)) , vwork, .true.)

            if (trim(fields_to_ocn(i)) == 'strsu_io') then
                iostrsu = vwork
            elseif (trim(fields_to_ocn(i)) == 'strsv_io') then
                iostrsv = vwork
            elseif (trim(fields_to_ocn(i)) == 'rain_io') then
                iorain = vwork
            elseif (trim(fields_to_ocn(i)) == 'snow_io') then
                iosnow = vwork
            elseif (trim(fields_to_ocn(i)) == 'stflx_io') then
                iostflx = vwork
            elseif (trim(fields_to_ocn(i)) == 'htflx_io') then
                iohtflx = vwork
            elseif (trim(fields_to_ocn(i)) == 'swflx_io') then
                ioswflx = vwork
            elseif (trim(fields_to_ocn(i)) == 'qflux_io') then
                ioqflux = vwork
            elseif (trim(fields_to_ocn(i)) == 'shflx_io') then
                ioshflx = vwork
            elseif (trim(fields_to_ocn(i)) == 'lwflx_io') then
                iolwflx = vwork
            elseif (trim(fields_to_ocn(i)) == 'runof_io') then
                iorunof = vwork
            elseif (trim(fields_to_ocn(i)) == 'press_io') then
                iopress = vwork
            elseif (trim(fields_to_ocn(i)) == 'aice_io') then
                ioaice = vwork
            elseif (trim(fields_to_ocn(i)) == 'melt_io') then
                iomelt = vwork
            elseif (trim(fields_to_ocn(i)) == 'form_io') then
                ioform = vwork
            elseif (trim(fields_to_ocn(i)) == 'licefw_io') then
                iolicefw = vwork
            elseif (trim(fields_to_ocn(i)) == 'licefh_io') then
                iolicefh = vwork
            else
                call abort_ice('ice: bad initialization array name '//trim(fields_to_ocn(i)))
            endif
      enddo
      if (my_task == master_task) call ice_close_nc(ncid_i2o)

#if defined(DEBUG)
      if (my_task == master_task) write(il_out,*) '(get_time0_i2o_fields) has read in 11 i2o fields.'
#endif
    else
        print *, 'CICE: (get_time0_i2o_fields_old) not found file *** ',fname
        stop 'CICE stopped -- Need time0 i2o data file.'
    endif

end subroutine get_time0_i2o_fields

!===============================================================================
subroutine get_u_star(fname)

implicit none

character*(*), intent(in) :: fname

    integer(kind=int_kind) :: ncid
    logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
#if defined(DEBUG)
  if (my_task == master_task) write(il_out,*) '(get_u_star) reading in initial u_star field ......'
#endif
  call ice_open_nc(fname, ncid)
  call ice_read_nc(ncid, 1, 'u_star', vwork, dbug)
  u_star0 = vwork
  if (my_task == master_task) call ice_close_nc(ncid)
#if defined(DEBUG)
  if (my_task == master_task) write(il_out,*) '(get_u_star) has read in u_star field.'
#endif
else
  print *, 'CICE: (get_u_star) not found file *** ',fname
  stop 'CICE stopped -- Need u_star restart datafile.'
endif

return
end subroutine get_u_star
!===============================================================================
subroutine get_sicemass(fname)

implicit none

character*(*), intent(in) :: fname

    integer(kind=int_kind) :: ncid
    logical :: dbug

dbug = .true.
if ( file_exist(fname) ) then
  call ice_open_nc(fname, ncid)
  call ice_read_nc(ncid, 1, 'sicemass', vwork, dbug)
  sicemass = vwork
  if (my_task == master_task) then
      call ice_close_nc(ncid)
#if defined(DEBUG)
      write(il_out,*) '(sicemass) reading in initial sicemass field ......'
      write(il_out,*) '(get_sicemass) has read in sicemass field.'
#endif
  endif
else
  if (my_task == master_task) then
#if defined(DEBUG)
      write(il_out,*) 'WARNING: (get_sicemass) not found file *** ',fname
      write(il_out,*) 'CICE: (get_sicemass) not found file *** ',fname
      write(il_out,*) 'CICE: Will estimate sicemass as iopress/(ioaice*gravit) *** '
#endif
  endif
  ! below crashes for some reason when cold start
  !where (ioaice>0.0) 
  !    sicemass = iopress/(ioaice*gravit)
  !elsewhere
  !    sicemass = 0.0
  !end where
  sicemass = 0.0
endif

return
end subroutine get_sicemass
 
!===============================================================================
subroutine newt_forcing_raw 

implicit none

! --- from atm: 
tair(:,:,:)  =  tair0(:,:,:)
fsw(:,:,:)   =  swflx0(:,:,:)
flw(:,:,:)   =  lwflx0(:,:,:)
uatm(:,:,:)  =  uwnd0(:,:,:)
vatm(:,:,:)  =  vwnd0(:,:,:)
qa(:,:,:)    =  qair0(:,:,:)
frain(:,:,:) =  rain0(:,:,:)
fsnow(:,:,:) =  snow0(:,:,:) 
press(:,:,:) =  press0(:,:,:)
runof(:,:,:) =  runof0(:,:,:)
calv(:,:,:) =  calv0(:,:,:)

! --- from ocean:
uocn(:,:,:) = ssuo(:,:,:)
vocn(:,:,:) = ssvo(:,:,:)
sst(:,:,:)  = ssto(:,:,:)
!make sure SST is 'all right' K==>C
if (maxval(sst).gt.200) then
    sst = sst - 273.15
endif
sss(:,:,:)     = ssso(:,:,:)
ss_tltx(:,:,:) = sslx(:,:,:)
ss_tlty(:,:,:) = ssly(:,:,:)
!frzmlt(:,:,:)  = pfmice(:,:,:) * frazil_factor / dt_cpl_io   !W/m^2 as required by cice.
frzmlt(:,:,:)  = pfmice(:,:,:)    
!frazil_factor and unit conversion done in ocean model (16/07/2008)

!12/03/2008: set max melt htflux allowed from ocn, eg, meltlimit = -200 (W/m^2)
!            meltlimit is read in from input_ice.nml
!20090320: set option 'limit_icemelt' in case no limit is needed when cice 'behaves'!
if (limit_icemelt) then
  frzmlt(:,:,:) = max(frzmlt(:,:,:), meltlimit)  
endif

end subroutine newt_forcing_raw

!===============================================================================
! output the last i2o forcing data, to be read in at the beginning of next run 
! by cice and sent to ocn immediately

subroutine save_time0_i2o_fields(fname, nstep)

implicit none

    character*(*), intent(in) :: fname
    integer(kind=int_kind), intent(in) :: nstep
    integer(kind=int_kind) :: ncid
    integer(kind=int_kind) :: i, ll, ilout

    if (my_task == 0) then 
        call create_ncfile(fname, ncid, nx_global, ny_global, ll=1, ilout=il_out)
        call write_nc_1Dtime(real(nstep), 1, 'time', ncid)
    endif

    do i=1, num_fields_to_ocn
        if (trim(fields_to_ocn(i)) == 'strsu_io') then
            vwork = iostrsu
        elseif (trim(fields_to_ocn(i)) == 'strsv_io') then
            vwork = iostrsv
        elseif (trim(fields_to_ocn(i)) == 'rain_io') then
            vwork = iorain
        elseif (trim(fields_to_ocn(i)) == 'snow_io') then
            vwork = iosnow
        elseif (trim(fields_to_ocn(i)) == 'stflx_io') then
            vwork = iostflx
        elseif (trim(fields_to_ocn(i)) == 'htflx_io') then
            vwork = iohtflx
        elseif (trim(fields_to_ocn(i)) == 'swflx_io') then
            vwork = ioswflx
        elseif (trim(fields_to_ocn(i)) == 'qflux_io') then
            vwork = ioqflux
        elseif (trim(fields_to_ocn(i)) == 'shflx_io') then
            vwork = ioshflx
        elseif (trim(fields_to_ocn(i)) == 'lwflx_io') then
            vwork = iolwflx
        elseif (trim(fields_to_ocn(i)) == 'runof_io') then
            vwork = iorunof
        elseif (trim(fields_to_ocn(i)) == 'press_io') then
            vwork = iopress
        elseif (trim(fields_to_ocn(i)) == 'aice_io') then
            vwork = ioaice
        elseif (trim(fields_to_ocn(i)) == 'melt_io') then
            vwork = iomelt
        elseif (trim(fields_to_ocn(i)) == 'form_io') then
            vwork = ioform
        elseif (trim(fields_to_ocn(i)) == 'licefw_io') then
            vwork = iolicefw
        elseif (trim(fields_to_ocn(i)) == 'licefh_io') then
            vwork = iolicefh
        else
            call abort_ice('ice: bad save array name '//trim(fields_to_ocn(i)))
        endif

        call gather_global(gwork, vwork, master_task, distrb_info)
        if (my_task == 0) then 
          call write_nc2D(ncid, fields_to_ocn(i), gwork, 2, nx_global, ny_global, 1, ilout=il_out)
        endif
    enddo

    if (my_task == 0) then
        call ncheck(nf_close(ncid), 'save_time0_i2o_fields: nf_close')
    endif

end subroutine save_time0_i2o_fields

!===============================================================================
subroutine save_u_star(fname, nstep)

! --- output the last u_star fields, to be read in at the beginning of next run 
!     by cice for roughness calculation if the gfdl surface flux scheme is used

implicit none

character*(*), intent(in) :: fname
    integer(kind=int_kind), intent(in) :: nstep
    integer(kind=int_kind) :: ncid
    integer(kind=int_kind) :: ll, ilout

if (my_task == 0) then
  call create_ncfile(fname, ncid, nx_global, ny_global, ll=1, ilout=il_out)
  call write_nc_1Dtime(real(nstep), 1, 'time', ncid)
endif

vwork = u_star0
call gather_global(gwork, vwork, master_task, distrb_info)

if (my_task == 0) then
  call write_nc2D(ncid, 'u_star', gwork, 2, nx_global, ny_global, 1, ilout=il_out)
  call ncheck(nf_close(ncid), 'save_u_star: nf_close')
endif

return
end subroutine save_u_star
!===============================================================================
subroutine save_sicemass(fname, nstep)

! --- output the last sicemass field, to be read in at the beginning of next run 
!     by cice for iopress calculation

implicit none

character*(*), intent(in) :: fname
    integer(kind=int_kind), intent(in) :: nstep
    integer(kind=int_kind) :: ncid
    integer(kind=int_kind) :: ll, ilout

if (my_task == 0) then
  call create_ncfile(fname, ncid, nx_global, ny_global, ll=1, ilout=il_out)
  call write_nc_1Dtime(real(nstep), 1, 'time', ncid)
endif

vwork = sicemass
call gather_global(gwork, vwork, master_task, distrb_info)

if (my_task == 0) then
  call write_nc2D(ncid, 'sicemass', gwork, 2, nx_global, ny_global, 1, ilout=il_out)
  call ncheck(nf_close(ncid), 'save_sicemass: nf_close')
endif

return
end subroutine save_sicemass

!===============================================================================

subroutine get_i2o_fluxes

    use ice_atmo, only : atmo_boundary_layer, atmbndy, atmo_boundary_const
 
implicit none

    real (kind=dbl_kind) :: &
   dtice       ! time step 

    integer(kind=int_kind) :: &
   i, j           , & ! horizontal indices
   iblk               !, & ! block index
   !ilo,ihi,jlo,jhi    ! beginning and end of physical domain

!ars599: 26032014 new code add in (CODE:Cdn)
    real (kind=dbl_kind), dimension(nx_block,ny_block) :: &
         delt    ,&         ! potential temperature difference   (K)
         delq    ,&         ! specific humidity difference   (kg/kg)
         shcoef  ,&         ! transfer coefficient for sensible heat
         lhcoef  ,&         ! transfer coefficient for latent heat
         Cdn_atm ,&         ! neutral drag coefficient
         Cdn_atm_ocn_n      ! ratio drag coeff / neutral drag coeff

!
    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
   sst_B, frzmlt_B  
   !B: we don't want CICE to update sst and frzmlt, which suppose to be from
   !   ocn model. *** May need revisit this part... and, rethink how the ocn 
   !   model 'top layer' (with free surface) changes its properties upon
   !   receiving the i2o forcing... e.g., the Heatfluxes only works on the  
   !   top layer no matter how thin it is? -- it could ca    use instaneous 
   !   extreme SST for some points, I assume......          Feb, 2008    ***
    real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
   swabs_ocn, pice
!

    integer(kind=int_kind) :: &
   ij        ! combined ij index

    integer(kind=int_kind), save :: &
   icells    ! number of ocean cells

    integer(kind=int_kind), dimension(nx_block*ny_block), save :: &
   indxi, indxj    ! compressed indices for ocean cells

!type (block) :: &
!   this_block           ! block information for current block

!
dtice = dt
swabs_ocn = 0.0
pice      = 0.0    
!initialization of array is critical -- otherwise it would be assigned HUGE value
!and lead to '* 253 Invalid operation ' (Fatal) error whenever it is involved in
!any operation such as a simple maxval(swabs_ocn), or even var = pice + 0.0 ....! 


do iblk = 1, nblocks

!-----------------------------------------------------------------
! Identify ocean cells.
! Set fluxes to zero in land cells.
!-----------------------------------------------------------------

   icells = 0
   do j = 1, ny_block
   do i = 1, nx_block
      if (tmask(i,j,iblk)) then
         icells = icells + 1
         indxi(icells) = i
         indxj(icells) = j
      else
         sst       (i,j,iblk) = c0
         frzmlt    (i,j,iblk) = c0
         flwout_ocn(i,j,iblk) = c0
         fsens_ocn (i,j,iblk) = c0
         flat_ocn  (i,j,iblk) = c0
         evap_ocn  (i,j,iblk) = c0
      endif
   enddo                  ! i
   enddo                  ! j
!-----------------------------------------------------------------
! Compute boundary layer quantities
!-----------------------------------------------------------------
   if (trim(atmbndy) == 'constant') then
!ars599: 26032014 new code add in (CODE:Cdn)
      call atmo_boundary_const (nx_block,  ny_block,   &
                                'ice',     icells,     &
                                indxi,     indxj,      &
                                uatm       (:,:,iblk), &   
                                vatm       (:,:,iblk), &   
                                wind       (:,:,iblk), &   
                                rhoa       (:,:,iblk), &
                                strairx_ocn(:,:,iblk), & 
                                strairy_ocn(:,:,iblk), & 
                                sst        (:,:,iblk), &
				potT       (:,:,iblk), &
                                Qa         (:,:,iblk), &
                                delt       (:,:),      &
				delq       (:,:),      &
                                lhcoef     (:,:),      &
                                shcoef     (:,:),      &
                                Cdn_atm    (:,:))  
   else ! default
      call atmo_boundary_layer (nx_block,  ny_block,   &
                                'ocn',     icells,     &
                                indxi,     indxj,      &
                                sst        (:,:,iblk), &    
                                potT       (:,:,iblk), &
                                uatm       (:,:,iblk), &   
                                vatm       (:,:,iblk), &   
                                wind       (:,:,iblk), &   
                                zlvl       (:,:,iblk), &   
                                Qa         (:,:,iblk), &     
                                rhoa       (:,:,iblk), &
                                strairx_ocn(:,:,iblk), & 
                                strairy_ocn(:,:,iblk), & 
                                Tref_ocn   (:,:,iblk), & 
                                Qref_ocn   (:,:,iblk), & 
                                delt       (:,:),      &    
                                delq       (:,:),      &
                                lhcoef     (:,:),      &
                                shcoef     (:,:),      &
                                Cdn_atm    (:,:),      &
				Cdn_atm_ocn_n     (:,:))     
!ars599: 26032014 new code add in (CODE:Cdn)
   endif

!-----------------------------------------------------------------
! Ocean  albedo
! For now, assume albedo = albocn in each spectral band.
!-----------------------------------------------------------------

!   alvdr_ocn(:,:,iblk) = albocn
!   alidr_ocn(:,:,iblk) = albocn
!   alvdf_ocn(:,:,iblk) = albocn
!   alidf_ocn(:,:,iblk) = albocn
!B: 20091118
   alvdr_ocn(:,:,iblk) = ocn_albedo2D(:,:,iblk)
   alidr_ocn(:,:,iblk) = ocn_albedo2D(:,:,iblk)
   alvdf_ocn(:,:,iblk) = ocn_albedo2D(:,:,iblk)
   alidf_ocn(:,:,iblk) = ocn_albedo2D(:,:,iblk)

!-----------------------------------------------------------------
! Compute ocean fluxes *BUT NOT update SST (let ocn model do it!)*
!     (see comments above)
!-----------------------------------------------------------------

   sst_B   (:,:,iblk) = sst   (:,:,iblk)
   frzmlt_B(:,:,iblk) = frzmlt(:,:,iblk) !actually no need to do so!

   call ocean_energy_budget_B (nx_block, ny_block,   &
                             dtice,    icells,     &
                             indxi,    indxj,      &
                             delt      (:,:),      &   
                             delq      (:,:),      &
                             lhcoef    (:,:),      &
                             shcoef    (:,:),      &
                             aice      (:,:,iblk), &
                             Tf        (:,:,iblk), &
                             swvdr     (:,:,iblk), &
                             swidr     (:,:,iblk), &
                             swvdf     (:,:,iblk), &
                             swidf     (:,:,iblk), &    
                             alvdr_ocn (:,:,iblk), &
                             alidr_ocn (:,:,iblk), &
                             alvdf_ocn (:,:,iblk), &
                             alidf_ocn (:,:,iblk), &
                             flw       (:,:,iblk), &
                             qdp       (:,:,iblk), &
                             hmix      (:,:,iblk), &
                             flwout_ocn(:,:,iblk), &
                             fsens_ocn (:,:,iblk), &
                             flat_ocn  (:,:,iblk), &
                             evap_ocn  (:,:,iblk), &
                             sst_B     (:,:,iblk), &
                             frzmlt_B  (:,:,iblk), &
                             swabs_ocn (:,:,iblk))

enddo                     ! iblk

!tioqflux = flat_ocn
!tioshflx = fsens_ocn
!tiolwflx = flw + flwout_ocn   !net lw flux (down) into ocean 
!tioswflx = swabs_ocn

if (gfdl_surface_flux) then

  !!call gather_global(gwork, u_star0, master_task, distrb_info)
  !!if (my_task == master_task) write(52,'(10e12.4)')gwork
  !overwrite the surface fluxes (sh, lh, lw and windstress) with the GFDL calculation
  call gfdl_ocean_fluxes ( fsens_ocn, flat_ocn, flwout_ocn, strairx_ocn, strairy_ocn )

endif

tioqflux = flat_ocn
tioshflx = fsens_ocn
tiolwflx = flw + flwout_ocn   !net lw flux (down) into ocean 
tioswflx = swabs_ocn

! === double checked with Siobhan the following mergeing under sea ice:  ===
!     Note: for those i2o fluxes which are already weighted with ice catagory
!     fractions in routine 'merge_fluxes', they must NOT be weighted here
!     AGAIN using the 'total' ice fraction aice!                     30/07/07
!     (eg, using 'fresh' instead of 'aice * fresh', 
!                'fswthru' instead of 'aice * fswthru' and so on. 
!     However, ice-ocean stresses are different: they are calculated (see 
!     ice_flux_in.F for details) with no aice(n) weighing. they thus must be
!     wighted by " * aice " here before being passed to ocean model. 08/08/07  
! -----------------------------------------------------------------------------

! 1) air- or ice-sea interface stress x
!  tiostrsu(:,:,:) = strairx_ocn(:,:,:) * (1. - aice(:,:,:)) + strocnxT(:,:,:) * aice(:,:,:)
! 2) air- or ice-sea interface stress y
!  tiostrsv(:,:,:) = strairy_ocn(:,:,:) * (1. - aice(:,:,:)) + strocnyT(:,:,:) * aice(:,:,:)
!
!!! *** BUG found here: strocnx/yT MUST change sign to '-'!		 11/03/2009 *** !!!
!
  tiostrsu(:,:,:) = strairx_ocn(:,:,:) * (1. - aice(:,:,:)) - strocnxT(:,:,:) * aice(:,:,:)
  tiostrsv(:,:,:) = strairy_ocn(:,:,:) * (1. - aice(:,:,:)) - strocnyT(:,:,:) * aice(:,:,:)

! 3) rainfall to ocean
  tiorain(:,:,:) = frain(:,:,:) * (1. - aice(:,:,:)) 
  
!!!  if (ice_fwflux) tiorain(:,:,:) = tiorain(:,:,:) + fresh(:,:,:)  !!!CH confirmed--ice_fwflux=.t.!!!
!!!
!!! Now ice_fwflux (melt/form) is sent to ocn seperately (see below)
!!!

! 4) snowfall to ocean
  tiosnow(:,:,:) = fsnow(:,:,:) * (1. - aice(:,:,:))

! 5) salt flux to ocean
  tiostflx(:,:,:) = fsalt(:,:,:) 

! 6) ice melting heatflux into ocean (fhnet renamed to fhocn)
  tiohtflx(:,:,:) = fhocn(:,:,:)                    

! 7) short wave radiation into ocean 
  tioswflx(:,:,:) = tioswflx(:,:,:) * (1. - aice(:,:,:)) + fswthru(:,:,:)

! 8) latent heat flux: 
!    *** Note sign changed here, to become positive out of ocean as required by mom4!
  tioqflux(:,:,:) = - tioqflux(:,:,:) * (1. - aice(:,:,:)) 

! 9) sensible heat flux: 
!    *** Note sign changed here, to become positive out of ocean as required by mom4!
  tioshflx(:,:,:) = - tioshflx(:,:,:) * (1. - aice(:,:,:))

!109)long wave radiation
  tiolwflx(:,:,:) = tiolwflx(:,:,:) * (1. - aice(:,:,:))

!11)runoff: relocated onto coastal grid points (pre-processed by Steve Phipps)
  tiorunof(:,:,:) = runof(:,:,:) + calv(:,:,:)
  tiolicefw(:,:,:) = 0.0
  tiolicefh(:,:,:) = 0.0

!12)pressure
!  if (my_task == 0)  write(il_out,*)'size of pice,     ',&
!                  size(pice(:,1,1)),    size(pice(1,:,1)),    size(pice(1,1,:))
!  if (my_task == 0) write(il_out,*)'size of sicemass: ',&
!                  size(sicemass(:,1,1)),size(sicemass(1,:,1)),size(sicemass(1,1,:)) 
  pice(:,:,:) = gravit * sicemass(:,:,:)  !sicemass = rho_ice x hi + rho_snow x hs (in m)
  !----------------------------------------------------------------------------
  ! Should we set limit to the ovelying ice pressure as suggested in MOM4 code?
  !(see ocean_sbc.F90) if yes, we may     use following 
  !pice(i,j) = min(pice(i,j), gravit*rhow*max_ice_thickness) 
  ! (note  rhow = 1026 kg/m^3 here, but mom4 instead uses rho0 = 1035 kg/m^3)
  ! No, let mom4 handle it (see ocean_sbc.F90)
  !----------------------------------------------------------------------------
  tiopress = press - 1.0e5   !as GFDL SIS, we     use patm anormaly!   
  if (ice_pressure_on) then
!sjm 20101118
!   tiopress(:,:,:) = press(:,:,:) + pice(:,:,:) * aice(:,:,:)
    tiopress(:,:,:) = pice(:,:,:) * aice(:,:,:)
  endif
  !----------------------------------------------------------------------------
  !as GFDL SIS, we     use patm anormaly and then add in the ice/snow pressure! 
  !29/11/2007
  !---------------------------------------------------------------------------- 
!13) ice coverage
  tioaice(:,:,:) = aice(:,:,:)

!--------------------------
!14) ice melt waterflux:
  tiomelt(:,:,:) = max(0.0,fresh(:,:,:))
!15) ice form waterflux:
  tioform(:,:,:) = min(0.0,fresh(:,:,:))

return
end subroutine get_i2o_fluxes

!=======================================================================
subroutine ocean_energy_budget_B ( nx_block,   ny_block,  &
                                dtice,      icells,    &
                                indxi,      indxj,     &
                                delt,       delq,      &
                                lhcoef,     shcoef,    &
                                Baice,       Tf,        &
                                swvdr,      swidr,     &
                                swvdf,      swidf,     &
                                alvdr_ocn,  alidr_ocn, &
                                alvdf_ocn,  alidf_ocn, &
                                Bflw,                   &
                                qdp,        hmix,      &
                                flwout_ocn, fsens_ocn, &
                                flat_ocn,   evap_ocn,  &
                                sst_B,      frzmlt_B,  &
                                swabs_ocn )

! Compute ocean energy budget. 
! 
! HISTORY: adapted from ocean_energy_budget, Feb 2008 by D. Bi for i2o_fluxes.
!
! !INPUT/OUTPUT PARAMETERS:
!
    integer(kind=int_kind), intent(in) :: &
   nx_block, ny_block, & ! block dimensions
   icells                ! number of cells that require atmo fluxes

    integer(kind=int_kind), dimension(nx_block*ny_block), &
   intent(in) :: &
   indxi, indxj    ! compressed i and j indices

    real (kind=dbl_kind), intent(in) :: &
   dtice           ! time step

    real (kind=dbl_kind), dimension(nx_block,ny_block), &
   intent(in) :: &
   delt         , & ! potential T difference   (K)
   delq         , & ! humidity difference      (kg/kg)
   shcoef       , & ! transfer coefficient for sensible heat
   lhcoef       , & ! transfer coefficient for latent heat
   Baice         , & ! fractional ice area
   Tf           , & ! ocean freezing temperature (C)
   swvdr        , & ! incoming shortwave, visible direct (W/m^2)
   swidr        , & ! incoming shortwave, near IR direct (W/m^2)
   swvdf        , & ! incoming shortwave, visible diff    use (W/m^2)
   swidf        , & ! incoming shortwave, near IR diff    use (W/m^2)
   alvdr_ocn    , & ! visible albedo, direct   (fraction)
   alidr_ocn    , & ! near-ir albedo, direct   (fraction)
   alvdf_ocn    , & ! visible albedo, diff    use  (fraction)
   alidf_ocn    , & ! near-ir albedo, diff    use  (fraction)
   Bflw          , & ! incoming longwave (W/m^2)
   hmix             ! ocean mixed layer depth (m)

    real (kind=dbl_kind), dimension(nx_block,ny_block), &
   intent(inout) :: &
   sst_B        , & ! sea surface temperature (C)
   qdp          , & ! deep ocean heat flux (W/m^2)
   frzmlt_B     , & ! freeze-melt potential (W/m^2)
   fsens_ocn    , & ! sensible heat flux (W/m^2)
   flat_ocn     , & ! latent heat flux   (W/m^2)
   flwout_ocn   , & ! outgoing longwave  (W/m^2)
   evap_ocn         ! evaporative vapor flux (kg/m^2/s)

    real (kind=dbl_kind), dimension(nx_block,ny_block), intent(out) :: &
   swabs_ocn

    real (kind=dbl_kind) :: &
   TsfK , & ! surface temperature (K)
   swabs    ! surface absorbed shortwave heat flux (W/m^2)

    real (kind=dbl_kind), parameter :: &
   frzmlt_max = c1000   ! max magnitude of frzmlt (W/m^2)

    integer(kind=int_kind) :: i, j, ij

do ij = 1, icells
   i = indxi(ij)
   j = indxj(ij)

   ! shortwave radiative flux
   swabs = (c1-alvdr_ocn(i,j)) * swvdr(i,j) &
         + (c1-alidr_ocn(i,j)) * swidr(i,j) &
         + (c1-alvdf_ocn(i,j)) * swvdf(i,j) &
         + (c1-alidf_ocn(i,j)) * swidf(i,j) 

   swabs_ocn(i,j) = swabs

   ! ocean surface temperature in Kelvin
   TsfK = sst_B(i,j) + Tffresh

   ! longwave radiative flux
   flwout_ocn(i,j) = -stefan_boltzmann * TsfK**4

   ! downward latent and sensible heat fluxes
   fsens_ocn(i,j) =  shcoef(i,j) * delt(i,j)
   flat_ocn (i,j) =  lhcoef(i,j) * delq(i,j)
   evap_ocn (i,j) = -flat_ocn(i,j) / Lvap

   !B: sst and frzmlt remain unchanged for a coupling period, thus the 
   !   following part is not needed:

   if (.false.) then

   ! Compute sst change due to exchange with atm/ice above
   ! Note: fhnet, fswthru are added in ice_therm_vertical.F
   sst_B(i,j) = sst_B(i,j) + &
        (fsens_ocn(i,j) + flat_ocn(i,j) + flwout_ocn(i,j) &
       + Bflw(i,j) + swabs) * (c1-Baice(i,j)) * dtice &
       / (cprho*hmix(i,j))

   ! adjust qdp if cooling of mixed layer would occur when sst <= Tf
   if (sst_B(i,j) <= Tf(i,j) .and. qdp(i,j) > c0) qdp(i,j) = c0

   ! computed T change due to exchange with deep layers:
   sst_B(i,j) = sst_B(i,j) - qdp(i,j)*dtice/(cprho*hmix(i,j))

   ! compute potential to freeze or melt ice
   frzmlt_B(i,j) = (Tf(i,j)-sst_B(i,j))*cprho*hmix(i,j)/dtice
   frzmlt_B(i,j) = min(max(frzmlt_B(i,j),-frzmlt_max),frzmlt_max)

   ! if sst is below freezing, reset sst to Tf
   if (sst_B(i,j) <= Tf(i,j)) sst_B(i,j) = Tf(i,j)

   endif

enddo                     ! ij

return
end subroutine ocean_energy_budget_B

!=======================================================================
subroutine gfdl_ocean_fluxes(sh,lh,lwo,taox,taoy)

    use ice_fileunits
    use ice_exit, only : abort_ice
    use ocean_rough_mod
    use surface_flux_mod

implicit none

!    real (kind=dbl_kind), intent(out), dimension(:,:,:) :: sh,lh,lwo,taox,taoy
    real (kind=dbl_kind), intent(out), dimension(nx_block, ny_block, nblocks) :: &
   sh,lh,lwo,taox,taoy

!real, dimension(size(sh(:,1,1)),size(sh(1,:,1))) :: &
real, dimension(nx_block, ny_block) :: &
   t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
   p_surf,    t_surf,     t_ca,      q_surf,                         &
   u_surf,    v_surf,                                                &
   rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
   flux_t,    flux_q,     flux_r,    flux_u,    flux_v,              &
   cd_m, cd_t, cd_q, w_atm, u_star, b_star, q_star, &
   dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
   dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm

   real :: dt_xx = 0.0   !not used!

!    logical, dimension(size(sh(:,1,1)),size(sh(1,:,1))) :: land, seawater, avail 
    logical, dimension(nx_block, ny_block) :: land, seawater, avail

!real, dimension(size(sh(:,1,1)),size(sh(1,:,1))) :: &
real, dimension(nx_block, ny_block) :: tv_atm, d_atm

real :: gust0 = 1.0

integer :: iblk

real, parameter :: d622   = rdgas/rvgas
real, parameter :: d378   = 1.-d622
real            :: d608   = d378/d622

z_atm = 10.0      !zlvl 
t_ca  = 222.222   !not really used (below we set avail = seawater)

!
!call gather_global(gwork, u_star0, master_task, distrb_info)
!if (my_task == master_task) write(53,'(10e12.4)')gwork
!

do iblk = 1, nblocks

   land(:,:) = .not. tmask(:,:,iblk)
   seawater(:,:) = tmask(:,:,iblk)
   avail(:,:) = seawater(:,:)

   u_star(:,:) = u_star0(:,:,iblk)

   !
   ! bi003:     use .../mom4p1/src/ice_param/ocean_rough.F90 code to calculate roughness.
   !---------------------------------------------------------------------------------------
   ! * note the surface friction velocity u_star required here is actually dependent 
   ! on roughness (see mo_drag), so we have to     use 'previous' step u_star to calculate 
   ! the current step roughness here. In this case u_star should be save in a restart 
   ! file for     use at the beginning in next run. (probably not crucial though, as per J McG)
   !--------------------------------------------------------------------------------------- 
   ! 
   call compute_ocean_roughness ( avail, u_star,  &
                                  rough_mom, rough_heat, rough_moist )

   !============= NOT SURE WHY 1.0 ===========!
   rough_scale = 1.             ! ??????????? !
   !rough_scale = rough_scale0  ! --- ??? --- !
   !==========================================!   

   gust = gust0

   where (avail)
 
   t_atm(:,:) = tair0(:,:,iblk)
   q_atm_in(:,:) = qair0(:,:,iblk)
   u_atm(:,:) = uwnd0(:,:,iblk)
   v_atm(:,:) = vwnd0(:,:,iblk)
   
   !'rough' estimate of virtual temperature 
   tv_atm(:,:) = tair0(:,:,iblk) * (1.0 + d608*qair0(:,:,iblk))
   !air density 10 meter height, 'rough', 'cos of using p_surf
   d_atm(:,:) = press0(:,:,iblk) / (rdgas * tv_atm(:,:))
   !estimate of 10 meter (zlvl) pressure
   p_atm(:,:) = press0(:,:,iblk) - d_atm(:,:) * gravit * 10.0

   p_surf(:,:) = press0(:,:,iblk)
   t_surf(:,:) = ssto(:,:,iblk)        !in Kelvin
   where (t_surf < 250)
      t_surf = t_surf + 273.15
   endwhere

   !q_surf is not really used for seawater!
   q_surf(:,:) = q_atm_in(:,:) 

   u_surf(:,:) = ssuo(:,:,iblk)
   v_surf(:,:) = ssvo(:,:,iblk)

   endwhere

   call surface_flux ( &
     t_atm,     q_atm_in,   u_atm,     v_atm,     p_atm,     z_atm,    &
     p_surf,    t_surf,     t_ca,      q_surf,                         &
     u_surf,    v_surf,                                                &
     rough_mom, rough_heat, rough_moist, rough_scale, gust,            &
     flux_t,    flux_q,     flux_r,    flux_u,    flux_v,              &
     cd_m,      cd_t,       cd_q,                                      &
     w_atm,     u_star,     b_star,     q_star,                        &
     dhdt_surf, dedt_surf,  dedq_surf,  drdt_surf,                     &
     dhdt_atm,  dedq_atm,   dtaudu_atm, dtaudv_atm,                    &
     dt_xx,     land,       seawater,  avail  )

! change sign for all fields as required by the mom4 ocean model! 
   sh(:,:,iblk) =  - flux_t(:,:)
   lh(:,:,iblk) =  - flux_q(:,:) * Lvap
   lwo(:,:,iblk) = - flux_r(:,:)
   taox(:,:,iblk) = - flux_u(:,:)
   taoy(:,:,iblk) = - flux_v(:,:)

   !u_star0 will be saved for restarting the next run  
   u_star0(:,:,iblk) = u_star(:,:)
   rough_mom0(:,:,iblk) = rough_mom(:,:)
   rough_heat0(:,:,iblk) = rough_heat(:,:)
   rough_moist0(:,:,iblk) = rough_moist(:,:)

enddo
return

end subroutine gfdl_ocean_fluxes

!===============================================================================
subroutine check_roughness(ncfile,nstep)

implicit none

character*(*), intent(in) :: ncfile
    integer(kind=int_kind), intent(in) :: nstep
    integer(kind=int_kind) :: ilout, ll
    integer(kind=int_kind) :: ncid,currstep
data currstep/0/
save currstep

currstep=currstep+1

if ( my_task == 0 .and. .not. file_exist(trim(ncfile)) )  then
  call create_ncfile(trim(ncfile),ncid,nx_global,ny_global,ll=1,ilout=il_out)
endif

if (my_task == 0) then
#if defined(DEBUG)
  write(il_out,*) 'opening file ',trim(ncfile), ' at nstep = ', nstep
#endif
  call ncheck(nf_open(trim(ncfile),nf_write,ncid), 'check_roughness: nf_open')
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

call gather_global(gwork, u_star0,   master_task, distrb_info)
!if (my_task == 0) write(61,'(10e12.4)') gwork
if (my_task == 0) call write_nc2D(ncid, 'u_star', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, rough_mom0,   master_task, distrb_info)
!if (my_task == 0) write(62,'(10e12.4)') gwork
if (my_task == 0) call write_nc2D(ncid, 'rough_mom', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, rough_heat0,   master_task, distrb_info)
!if (my_task == 0) write(63,'(10e12.4)') gwork
if (my_task == 0) call write_nc2D(ncid, 'rough_heat', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, rough_moist0,   master_task, distrb_info)
!if (my_task == 0) write(64,'(10e12.4)') gwork
if (my_task == 0) call write_nc2D(ncid, 'rough_moist', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid), 'check_roughness: nf_close')

return
end subroutine check_roughness

!============================================================================
subroutine check_a2i_fields(ncfile,nstep)

implicit none

character*(*), intent(in) :: ncfile
    integer(kind=int_kind), intent(in) :: nstep
    integer(kind=int_kind) :: ncid,currstep,ll,ilout
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(trim(ncfile)) ) then
  call create_ncfile(trim(ncfile),ncid,nx_global,ny_global,ll=1,ilout=il_out)
endif

if (my_task == 0) then
#if defined(DEBUG)
  write(il_out,*) 'opening file ',trim(ncfile),' at nstep = ', nstep
#endif
  call ncheck(nf_open(trim(ncfile),nf_write,ncid), 'check_a2i_fields: nf_open' )
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

call gather_global(gwork, tair0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'tair0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, swflx0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'swflx0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, lwflx0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'lwflx0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, uwnd0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'uwnd0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, vwnd0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'vwnd0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, qair0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'qair0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, rain0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'rain0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, snow0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'snow0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, runof0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'runof0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, calv0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'calv0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, press0, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'press0', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid), 'check_a2i_fields: nf_close')

return
end subroutine check_a2i_fields

!============================================================================
subroutine check_i2o_fields(ncfile,nstep, scale)

implicit none

character*(*), intent(in) :: ncfile
    integer(kind=int_kind), intent(in) :: nstep
real, intent(in) :: scale
    integer(kind=int_kind) :: ncid,currstep, ll, ilout
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(trim(ncfile)) ) then
  call create_ncfile(trim(ncfile),ncid,nx_global,ny_global,ll=1,ilout=il_out)
endif

if (my_task == 0) then
#if defined(DEBUG)
  write(il_out,*) 'opening file ', trim(ncfile), ' at nstep = ', nstep
#endif
  call ncheck(nf_open('fields_i2o_in_ice.nc',nf_write,ncid), &
              'check_i2o_fields: nf_open')
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

call gather_global(gwork, scale*iostrsu, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iostrsu', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*iostrsv, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iostrsv', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*iorain, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iorain', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*iosnow, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iosnow', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*iostflx, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iostflx', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*iohtflx, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iohtflx', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*ioswflx, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ioswflx', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*ioqflux, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ioqflux', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*ioshflx, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ioshflx', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*iolwflx, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iolwflx', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*iorunof, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iorunof', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*iopress, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iopress', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*ioaice,  master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ioaice',  gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
!!!
call gather_global(gwork, scale*iomelt,  master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'iomelt',  gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, scale*ioform,  master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ioform',  gwork, 2, nx_global,ny_global,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid), 'check_i2o_fields: nf_close')

return
end subroutine check_i2o_fields

!============================================================================
subroutine check_o2i_fields(ncfile,nstep)

implicit none

character*(*), intent(in) :: ncfile
    integer(kind=int_kind), intent(in) :: nstep
    integer(kind=int_kind) :: ncid,currstep, ilout, ll
data currstep/0/
save currstep

currstep=currstep+1

if (my_task == 0 .and. .not. file_exist(trim(ncfile)) ) then
  call create_ncfile(trim(ncfile),ncid,nx_global,ny_global,ll=1,ilout=il_out)
endif

if (my_task == 0) then
#if defined(DEBUG)
  write(il_out,*) 'opening file ',trim(ncfile),' at nstep = ', nstep
#endif
  call ncheck(nf_open(trim(ncfile),nf_write,ncid), 'check_o2i_fields: nf_open')
  call write_nc_1Dtime(real(nstep),currstep,'time',ncid)
end if

call gather_global(gwork, ssto, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ssto', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, ssso, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ssso', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, ssuo, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ssuo', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, ssvo, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ssvo', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, sslx, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'sslx', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, ssly, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'ssly', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, pfmice, master_task, distrb_info)
if (my_task == 0) call write_nc2D(ncid, 'pfmice', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid), 'check_o2i_fields: nf_close')

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
  call create_ncfile(ncfilenm,ncid,nx_global,ny_global,ll=1,ilout=il_out)
endif

if (my_task == 0) then
#if defined(DEBUG)
  write(il_out,*) 'opening ncfile at nstep ', ncfilenm,  currstep
#endif
  call ncheck(nf_open(ncfilenm, nf_write,ncid), 'check_frzmlt_sst: nf_open')
  call write_nc_1Dtime(real(currstep),currstep,'time',ncid)
end if

call gather_global(gwork, sst, master_task, distrb_info) 
if (my_task == 0) call write_nc2D(ncid, 'sst', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)
call gather_global(gwork, frzmlt, master_task, distrb_info) 
if (my_task == 0) call write_nc2D(ncid, 'frzmlt', gwork, 2, nx_global,ny_global,currstep,ilout=il_out)

if (my_task == 0) call ncheck(nf_close(ncid), 'check_frzmlt_sst: nf_close')

return
end subroutine check_frzmlt_sst

!============================================================================
subroutine new_freezingT
 
!ars599: 26032014 temperally removed
!	new code has no such variable
!	so remove the if
!        if (trim(Tfrzpt) == 'constant') then
           Tf (:,:,:) = Tocnfrz ! deg C
!        else ! default:  Tfrzpt = 'linear_S'
           Tf    (:,:,:) = -depressT*sss(:,:,:)  ! freezing temp (C)
!        endif
 
end subroutine new_freezingT

!============================================================================
function file_exist (file_name)
!
character(len=*), intent(in) :: file_name
    logical  file_exist

file_exist = .false.
if (len_trim(file_name) == 0) return
if (file_name(1:1) == ' ')    return

inquire (file=trim(file_name), exist=file_exist)

end function file_exist

!============================================================================

end module cpl_forcing_handler
