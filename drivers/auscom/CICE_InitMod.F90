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

#ifdef AusCOM
      use accessom2_mod, only : accessom2_type => accessom2
      use cpl_parameters
      use cpl_parameters, only : read_namelist_parameters, accessom2_config_dir
      use cpl_forcing_handler, only : get_time0_sstsss, get_u_star
      use cpl_interface , only : prism_init, init_cpl, il_commlocal
      use cpl_interface, only: coupler
      use cpl_arrays_setup, only : gwork, u_star0
      use ice_gather_scatter

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

      subroutine CICE_Initialize(accessom2)
        type(accessom2_type), intent(out) :: accessom2

   !--------------------------------------------------------------------
   ! model initialization
   !--------------------------------------------------------------------

      call cice_init(accessom2)

      end subroutine CICE_Initialize

!=======================================================================
!
!  Initialize CICE model.

      subroutine cice_init(accessom2)

      use ice_aerosol, only: faero_default
      use ice_algae, only: get_forcing_bgc
      use ice_calendar, only: dt, npt, dt_dyn, time, istep, istep1, write_ic, &
          init_calendar, calendar, idate, month
!ars599: 27032014
      use ice_communicate, only: MPI_COMM_ICE
      use ice_communicate, only: init_communicate
      use ice_communicate, only: my_task, master_task
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
      use ice_restart_shared, only: restart_dir, input_dir
#ifdef popcice
      use drv_forcing, only: sst_sss
#endif
      use version_mod, only: CICE_COMMIT_HASH

      type(accessom2_type), intent(inout) :: accessom2

      integer(kind=int_kind) :: idate_save


      call read_namelist_parameters()

      ! initial setup for message passing
      call init_communicate()

      call prism_init(trim(accessom2_config_dir))
      MPI_COMM_ICE = il_commlocal

      call init_fileunits       ! unit numbers

      ! Initialise libaccessom2
      call accessom2%init('cicexx', config_dir=trim(accessom2_config_dir))

      ! Tell libaccessom2 about any global configs/state
      call accessom2%set_cpl_field_counts(num_atm_to_ice_fields=num_fields_from_atm, &
                                          num_ice_to_ocean_fields=num_fields_to_ocn, &
                                          num_ocean_to_ice_fields=num_fields_from_ocn)

      ! Synchronise accessom2 configuration between all models and PEs
      call accessom2%sync_config(coupler)

      ! Use accessom2 configuration
      call input_data(accessom2%get_forcing_start_date_array(), &
                      accessom2%get_cur_exp_date_array(), &
                      accessom2%get_seconds_since_cur_exp_year(), &
                      accessom2%get_total_runtime_in_seconds(), &
                      accessom2%get_ice_ocean_timestep(), &
                      accessom2%get_calendar_type())

      if (trim(runid) == 'bering') call check_finished_file
      call init_zbgc            ! vertical biogeochemistry namelist

      call init_domain_blocks   ! set up block decomposition
      call init_grid1           ! domain distribution

      call init_ice_timers      ! initialize all timers
      call ice_timer_start(timer_total)   ! start timing entire run
      call init_grid2           ! grid variables

#ifdef AusCOM
     ! initialize message passing, pass in total runtime in seconds and field
     ! coupling timesteps for oasis.
      call init_cpl(int(npt*dt), accessom2%get_coupling_field_timesteps())
#endif
      call init_calendar        ! initialize some calendar stuff
      call init_hist (dt)       ! initialize output history file

      call get_cpl_timecontrol(accessom2%get_atm_ice_timestep(), &
                               accessom2%get_ice_ocean_timestep())

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

#ifdef AusCOM
#if defined(DEBUG)
      write(il_out,*)' CICE: calendar called!'
#endif

      if (runtype == 'initial') then
        nrec = month - 1            !month is from calendar
        if (nrec == 0) nrec = 12 
        call get_time0_sstsss(trim(input_dir)//'monthly_sstsss.nc', nrec)
#if defined(DEBUG)
        write(il_out,*) 'CICE called  get_time0_sstsss. my_task = ',my_task
#endif
      endif
      !the read in sst/sss determines the initial ice state (in init_state)
      !which is overwritten by call to restartfile if restart=.t.
      !
      !20100111: get the surface friction velocity for gfdl surface flux calculation
      !          (roughness calculation requires last time step u_star...)
      if (gfdl_surface_flux) then
         call get_u_star(trim(input_dir)//'u_star.nc')
      endif  

#else
      call init_forcing_ocn(dt) ! initialize sss and sst from data
#endif

      call init_state           ! initialize the ice state
      call init_transport       ! initialize horizontal transport
      call ice_HaloRestore_init ! restored boundary conditions

      call init_restart         ! initialize restart variables

#ifdef AusCOM
#if defined(DEBUG)
      write(il_out,*) 'CICE (cice_init) 2      time = ', my_task, time
      write(il_out,*) 'CICE (cice_init) 2     idate = ', my_task, idate
#endif
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
         call calendar(time)    ! at the end of the first timestep

#if defined(DEBUG)
      write(il_out,*) 'CICE (cice_init) 3     time = ', my_task, time
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

#ifndef AusCOM
      call init_flux_atm        ! initialize atmosphere fluxes sent to coupler
      call init_flux_ocn        ! initialize ocean fluxes sent to coupler
#endif

      ! Print out my version
      if (my_task == master_task) then
          print*, CICE_COMMIT_HASH
          call accessom2%print_version_info()
      endif

      if (write_ic) call accum_hist(dt) ! write initial conditions 

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
               call init_lvl(nx_block, ny_block, ncat, &
                    trcrn(:,:,nt_alvl,:,iblk), trcrn(:,:,nt_vlvl,:,iblk))
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
