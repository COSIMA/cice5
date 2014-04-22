!  SVN:$Id: CICE_InitMod.F90 746 2013-09-28 22:47:56Z eclare $
!=======================================================================
!
!  This module contains the CICE initialization routine that sets model
!  parameters and initializes the grid and CICE state variables.
!
!  authors Elizabeth C. Hunke, LANL
!          William H. Lipscomb, LANL
!          Philip W. Jones, LANL
!
! 2006: Converted to free form source (F90) by Elizabeth Hunke
! 2008: E. Hunke moved ESMF code to its own driver

      module CICE_InitMod

      use ice_kinds_mod

!ars599: 11042014: remove details of forcing handler, interface, arrays setup, and gather scatter
#ifdef AusCOM
      use cpl_parameters
      use cpl_forcing_handler
      use cpl_interface
!ars599: 27032014: defind my_task
      use ice_communicate, only: my_task
#endif

      implicit none
      private
      public :: CICE_Initialize, cice_init
      save

#ifdef AusCOM
      integer :: nrec
#endif

!=======================================================================

      contains

!=======================================================================

!  Initialize the basic state, grid and all necessary parameters for
!  running the CICE model.  Return the initial state in routine
!  export state.
!  Note: This initialization driver is designed for standalone and
!        CCSM-coupled applications.  For other
!        applications (e.g., standalone CAM), this driver would be
!        replaced by a different driver that calls subroutine cice_init,
!        where most of the work is done.

      subroutine CICE_Initialize

   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------

      call cice_init

      end subroutine CICE_Initialize

!=======================================================================
!
!  Initialize CICE model.

      subroutine cice_init

      use ice_aerosol, only: faero_default
      use ice_algae, only: get_forcing_bgc
      use ice_calendar, only: dt, dt_dyn, time, istep, istep1, write_ic, &
          init_calendar, calendar, idate, month
!ars599: 27032014
      use ice_communicate, only: MPI_COMM_ICE
      use ice_communicate, only: init_communicate
      use ice_diagnostics, only: init_diags
      use ice_domain, only: init_domain_blocks
      use ice_dyn_eap, only: init_eap
      use ice_dyn_shared, only: kdyn, init_evp
      use ice_fileunits, only: init_fileunits
      use ice_flux, only: init_coupler_flux, init_history_therm, &
          init_history_dyn, init_flux_atm, init_flux_ocn
      use ice_forcing, only: init_forcing_ocn, init_forcing_atmo, &
          get_forcing_atmo, get_forcing_ocn
      use ice_grid, only: init_grid1, init_grid2
      use ice_history, only: init_hist, accum_hist
      use ice_restart_shared, only: restart, runid, runtype
      use ice_init, only: input_data, init_state
      use ice_itd, only: init_itd
      use ice_kinds_mod
      use ice_restoring, only: ice_HaloRestore_init
      use ice_shortwave, only: init_shortwave
      use ice_state, only: tr_aero
      use ice_therm_vertical, only: init_thermo_vertical
      use ice_timers, only: timer_total, init_ice_timers, ice_timer_start
      use ice_transport_driver, only: init_transport
      use ice_zbgc, only: init_zbgc
      use ice_zbgc_shared, only: skl_bgc
#ifdef popcice
      use drv_forcing, only: sst_sss
#endif

#ifdef AusCOM
      integer(kind=int_kind) :: idate_save
#endif

      call init_communicate     ! initial setup for message passing
#ifdef AusCOM 
      call prism_init		! called in init_communicate	
      MPI_COMM_ICE = il_commlocal
!      call init_cpl     ! initialize message passing
      call get_cpl_timecontrol
      write(il_out,*)' CICE (cice_init) 1    jobnum = ',jobnum
      write(il_out,*)' CICE (cice_init) 1   inidate = ',inidate
      write(il_out,*)' CICE (cice_init) 1 init_date = ',init_date
      write(il_out,*)' CICE (cice_init) 1  runtime0 = ',runtime0
      write(il_out,*)' CICE (cice_init) 1   runtime = ',runtime
      write(il_out,*)' CICE (cice_init) 1     idate = ',my_task, idate
      !write(il_out,*)' CICE (cice_init) 1   runtype = ',runtype
#endif
      call init_fileunits       ! unit numbers

      call input_data           ! namelist variables
      if (trim(runid) == 'bering') call check_finished_file
      call init_zbgc            ! vertical biogeochemistry namelist

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution
#ifdef AusCOM
      call init_cpl     ! initialize message passing
#endif
      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables

      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file

      if (kdyn == 2) then
         call init_eap (dt_dyn) ! define eap dynamics parameters, variables
      else                      ! for both kdyn = 0 or 1
         call init_evp (dt_dyn) ! define evp dynamics parameters, variables
      endif

      call init_coupler_flux    ! initialize fluxes exchanged with coupler
#ifdef popcice
      call sst_sss              ! POP data for CICE initialization
#endif 
      call init_thermo_vertical ! initialize vertical thermodynamics
      call init_itd             ! initialize ice thickness distribution
      call calendar(time)       ! determine the initial date

!ars599: 11042014: remove most of the lines based on cice4.1_fm
!	however try to keep idate_save
#ifdef AusCOM
      idate_save = idate  !save for late re-set in case 'restart' is used for jobnum=1
                          !and mess up the calendar idate for this exp...!
      if (jobnum==1) then 
        nrec = month - 1        !month is from calendar
        if (nrec == 0) nrec = 12 
        call get_time0_sstsss(trim(inputdir)//'/monthly_sstsss.nc', nrec)
      endif
      !the read in sst/sss determines the initial ice state (in init_state)
#else
      call init_forcing_ocn(dt) ! initialize sss and sst from data
#endif
      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport
!      call ice_HaloRestore_init ! restored boundary conditions
!ars599: 11042014: note:
!	the "if (runtype == 'continue') then ! start from core restart file"
!	from dhb599 might add in the new code "init_restart together with
!	tracer part
!	however transport and Halorestore might keep in OM but not in CM??

      call init_restart         ! initialize restart variables

#ifdef AusCOM
      write(il_out,*) 'CICE (cice_init) 2      time = ', my_task, time
      write(il_out,*) 'CICE (cice_init) 2  runtime0 = ', my_task, runtime0
      write(il_out,*) 'CICE (cice_init) 2     idate = ', my_task, idate
 
      if (jobnum == 1 ) then
        time = 0.0            !NOTE, the first job must be set back to 0 and 
        idate = idate_save    !idate back to the 'initial' value, in any case
      endif 
#endif

      call init_diags           ! initialize diagnostic output points
      call init_history_therm   ! initialize thermo history variables
      call init_history_dyn     ! initialize dynamic history variables

      ! Initialize shortwave components using swdn from previous timestep 
      ! if restarting. These components will be scaled to current forcing 
      ! in prep_radiation.
      if (trim(runtype) == 'continue' .or. restart) &
         call init_shortwave    ! initialize radiative transfer

         istep  = istep  + 1    ! update time step counters
         istep1 = istep1 + 1
         time = time + dt       ! determine the time and date
#ifndef AusCOM
         call calendar(time)    ! at the end of the first timestep
#else
!ars599: 26032014 original code
!         call calendar(time)    ! at the end of the first timestep
         call calendar(time-runtime0)
      write(il_out,*) 'CICE (cice_init) 3     time = ', my_task, time
      write(il_out,*) 'CICE (cice_init) 3 runtime0 = ', my_task, runtime0
      write(il_out,*) 'CICE (cice_init) 3    idate = ', my_task, idate
#endif

   !--------------------------------------------------------------------
   ! coupler communication or forcing data initialization
   !--------------------------------------------------------------------

      call init_forcing_atmo    ! initialize atmospheric forcing (standalone)

#ifndef coupled
      call get_forcing_atmo     ! atmospheric forcing from data
      call get_forcing_ocn(dt)  ! ocean forcing from data
!      if (tr_aero) call faero_data          ! aerosols
      if (tr_aero) call faero_default ! aerosols
      if (skl_bgc) call get_forcing_bgc
#endif

      if (runtype == 'initial' .and. .not. restart) &
         call init_shortwave    ! initialize radiative transfer using current swdn

!ars599: 09042014: add in
!	based on cice4.1_fm
!20091020
!#ifndef AusCOM
      call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
      call init_flux_ocn        ! initialize ocean fluxes sent to coupler
!#endif

!ars599: 11042014: ice_write_hist is no longer there now change to accum_hist
!	so wrapup this line n use the new code
!dhb599 20111128: this call is moved here from 'downstair', because it *re-initilaise*
!     the o2i fileds read in from o2i.nc (e.g., sst=>-1.84 sss=>34.) !!!  
!      call ice_write_hist(dt)   ! write initial conditions if write_ic = T
      if (write_ic) call accum_hist(dt) ! write initial conditions 

#ifdef AusCOM
      write(il_out,*)' calling init_mocn_fields_4_i2a at time_sec = ',0
      !call initialize_mice_fields_4_i2a
      call initialize_mocn_fields_4_i2a

      ! for continue runs, need restart o2i forcing fields and time-averaged ice 
      ! variables ('mice')saved at the end of last run from ice models; 
      ! for initial run, pre-processed o2i (and maybe mice) fields are required.
!      call get_restart_o2i('o2i.nc')
      call get_restart_o2i(trim(restartdir)//'/o2i.nc')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !if no lag for ice to atm coupling, then cice has to read restart file i2a.nc and 
      !put the data to atm. the call is not needed if there is lag for ice2atm coupling
      !must call after get_restart_o2i(), by which the ocn_sst ect are read in and re-used by put_restart_i2a()  
!      call put_restart_i2a('i2a.nc', 0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      if ( file_exist('CICE_restart/mice.nc') ) then
      if ( file_exist(trim(restartdir)//'/mice.nc') ) then
        !for continue runs, mice data MUST be available.
!        call get_restart_mice('CICE_restart/mice.nc')
        call get_restart_mice(trim(restartdir)//'/mice.nc')
      else
write(6,*)'*** CICE WARNING: No initial mice.nc data available here! **'
write(6,*)'*** CICE WARNING: ALL mice variables will be set to ZERO! **'
write(6,*)'*** CICE WARNING: This is allowed for the init run ONLY ! **' 
      endif
      if (use_core_runoff) then
!         call get_core_runoff('CICE_input/core_runoff_regrid.nc',&
         call get_core_runoff(trim(inputdir)//'/core_runoff_regrid.nc',&
                              'runoff',1)
      endif

        write(il_out,*)' calling ave_ocn_fields_4_i2a at time_sec = ',0 !time_sec
        call time_average_ocn_fields_4_i2a  !accumulate/average ocn fields needed for IA coupling

      !get a2i fields and then set up initial SBC for ice
      !call from_atm(0)
      !call get_sbc_ice
      !now "most of" the IC of this run are 'proper' for "call ice_write_hist"
#endif

!ars599: 11042014: ice_write_hist is no longer there now change to accum_hist
!	so wrapup this line n markout
!dhb599: 20111128: the following call is moved 'upstair'
!      if (write_ic) call accum_hist(dt) ! write initial conditions 

      end subroutine cice_init

!=======================================================================

      subroutine init_restart

      use ice_aerosol, only: init_aerosol
      use ice_age, only: init_age, restart_age, read_restart_age
      use ice_blocks, only: nx_block, ny_block
      use ice_brine, only: init_hbrine
      use ice_calendar, only: time, calendar
      use ice_domain, only: nblocks
      use ice_domain_size, only: ncat
      use ice_dyn_eap, only: read_restart_eap
      use ice_dyn_shared, only: kdyn
      use ice_firstyear, only: init_fy, restart_FY, read_restart_FY
      use ice_flux, only: sss
      use ice_init, only: ice_ic
      use ice_lvl, only: init_lvl, restart_lvl, read_restart_lvl
      use ice_meltpond_cesm, only: init_meltponds_cesm, &
          restart_pond_cesm, read_restart_pond_cesm
      use ice_meltpond_lvl, only: init_meltponds_lvl, &
          restart_pond_lvl, read_restart_pond_lvl, dhsn
      use ice_meltpond_topo, only: init_meltponds_topo, &
          restart_pond_topo, read_restart_pond_topo
      use ice_restart_shared, only: runtype, restart
      use ice_restart_driver, only: restartfile, restartfile_v4
      use ice_state, only: tr_iage, tr_FY, tr_lvl, tr_pond_cesm, &
          tr_pond_lvl, tr_pond_topo, tr_aero, trcrn, &
          nt_iage, nt_FY, nt_alvl, nt_vlvl, nt_apnd, nt_hpnd, nt_ipnd, tr_brine
      use ice_zbgc, only: init_bgc
      use ice_zbgc_shared, only: skl_bgc

      integer(kind=int_kind) :: iblk

      if (trim(runtype) == 'continue') then 
         ! start from core restart file
         call restartfile()           ! given by pointer in ice_in
!ars599: 11042014: markout call calendar
!	according to dhb599 initmod at cice4.1_fm
         call calendar(time)          ! update time parameters
         if (kdyn == 2) call read_restart_eap ! EAP
      else if (restart) then          ! ice_ic = core restart file
         call restartfile (ice_ic)    !  or 'default' or 'none'
         !!! uncomment to create netcdf
         ! call restartfile_v4 (ice_ic)  ! CICE v4.1 binary restart file
         !!! uncomment if EAP restart data exists
         ! if (kdyn == 2) call read_restart_eap
      endif         

      ! tracers
      ! ice age tracer   
      if (tr_iage) then 
         if (trim(runtype) == 'continue') &
              restart_age = .true.
         if (restart_age) then
            call read_restart_age
         else
            do iblk = 1, nblocks 
               call init_age(nx_block, ny_block, ncat, trcrn(:,:,nt_iage,:,iblk))
            enddo ! iblk
         endif
      endif
      ! first-year area tracer
      if (tr_FY) then
         if (trim(runtype) == 'continue') restart_FY = .true.
         if (restart_FY) then
            call read_restart_FY
         else
            do iblk = 1, nblocks 
               call init_FY(nx_block, ny_block, ncat, trcrn(:,:,nt_FY,:,iblk))
            enddo ! iblk
         endif
      endif
      ! level ice tracer
      if (tr_lvl) then
         if (trim(runtype) == 'continue') restart_lvl = .true.
         if (restart_lvl) then
            call read_restart_lvl
         else
            do iblk = 1, nblocks 
!ars599: 11042014: markout call init_lvl
!	according to dhb599 initmod at cice4.1_fm
!               call init_lvl(nx_block, ny_block, ncat, &
!                    trcrn(:,:,nt_alvl,:,iblk), trcrn(:,:,nt_vlvl,:,iblk))
            enddo ! iblk
         endif
      endif
      ! CESM melt ponds
      if (tr_pond_cesm) then
         if (trim(runtype) == 'continue') &
              restart_pond_cesm = .true.
         if (restart_pond_cesm) then
            call read_restart_pond_cesm
         else
            do iblk = 1, nblocks 
               call init_meltponds_cesm(nx_block, ny_block, ncat, &
                    trcrn(:,:,nt_apnd,:,iblk), trcrn(:,:,nt_hpnd,:,iblk))
            enddo ! iblk
         endif
      endif
      ! level-ice melt ponds
      if (tr_pond_lvl) then
         if (trim(runtype) == 'continue') &
              restart_pond_lvl = .true.
         if (restart_pond_lvl) then
            call read_restart_pond_lvl
         else
            do iblk = 1, nblocks 
               call init_meltponds_lvl(nx_block, ny_block, ncat, &
                    trcrn(:,:,nt_apnd,:,iblk), trcrn(:,:,nt_hpnd,:,iblk), &
                    trcrn(:,:,nt_ipnd,:,iblk), dhsn(:,:,:,iblk))
            enddo ! iblk
         endif
      endif
      ! topographic melt ponds
      if (tr_pond_topo) then
         if (trim(runtype) == 'continue') &
              restart_pond_topo = .true.
         if (restart_pond_topo) then
            call read_restart_pond_topo
         else
            do iblk = 1, nblocks 
               call init_meltponds_topo(nx_block, ny_block, ncat, &
                    trcrn(:,:,nt_apnd,:,iblk), trcrn(:,:,nt_hpnd,:,iblk), &
                    trcrn(:,:,nt_ipnd,:,iblk))
            enddo ! iblk
         endif ! .not restart_pond
      endif
      if (tr_aero)  call init_aerosol ! ice aerosol
      if (tr_brine) call init_hbrine  ! brine height tracer
      if (skl_bgc)  call init_bgc     ! biogeochemistry

      end subroutine init_restart

!=======================================================================
!
! Check whether a file indicating that the previous run finished cleanly
! If so, then do not continue the current restart.  This is needed only 
! for runs on machine 'bering' (set using runid = 'bering').
!
!  author: Adrian Turner, LANL

      subroutine check_finished_file()

      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_restart_shared, only: restart_dir

      character(len=char_len_long) :: filename
      logical :: lexist = .false.

      if (my_task == master_task) then
           
         filename = trim(restart_dir)//"finished"
         inquire(file=filename, exist=lexist)
         if (lexist) then
            call abort_ice("Found already finished file - quitting")
         end if

      endif

      end subroutine check_finished_file

!=======================================================================

      end module CICE_InitMod

!=======================================================================
