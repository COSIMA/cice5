!  SVN:$Id: CICE_RunMod.F90 746 2013-09-28 22:47:56Z eclare $
!=======================================================================
!
!  Main driver for time stepping of CICE.
!
!  authors Elizabeth C. Hunke, LANL
!          Philip W. Jones, LANL
!          William H. Lipscomb, LANL
!
! 2006 ECH: moved exit timeLoop to prevent execution of unnecessary timestep
! 2006 ECH: Streamlined for efficiency 
! 2006 ECH: Converted to free source form (F90)
! 2007 BPB: Modified Delta-Eddington shortwave interface
! 2008 ECH: moved ESMF code to its own driver

      module CICE_RunMod

      use ice_kinds_mod

#ifdef AusCOM 
      !For stuff in this AusCOM's own driver the "#ifdef AusCOM" is NOT needed!
      !but for consistency with the code in other places, we keep it anyway ...
      !...to "indentify" the modification to the original code, easier 
      !...to locate future code update                              Aug. 2008  
      
      use cpl_parameters
      use cpl_arrays_setup
      use cpl_interface
      use cpl_forcing_handler
      use cpl_interface, only : write_boundary_checksums
      use accessom2_mod, only : accessom2_type => accessom2
      use simple_timer_mod, only : simple_timer_type => simple_timer
      use logger_mod, only : logger_type => logger
#endif

      implicit none
      private
      public :: CICE_Run, ice_step
      save

!=======================================================================

      contains

!=======================================================================
!
!  This is the main driver routine for advancing CICE forward in time.
!
!  author Elizabeth C. Hunke, LANL
!         Philip W. Jones, LANL
!         William H. Lipscomb, LANL

      subroutine CICE_Run(accessom2, logger)

      use ice_aerosol, only: faero_default
      use ice_algae, only: get_forcing_bgc
      use ice_calendar, only: istep, istep1, time, dt, npt, stop_now, calendar
      use ice_communicate, only : my_task, master_task
#ifdef AusCOM
!ars599: 27032014 add in
      use ice_calendar, only: month, mday
      use ice_forcing, only: get_forcing_atmo, get_forcing_ocn, &
          get_forcing_atmo_ready
#endif
      use ice_flux, only: init_flux_atm, init_flux_ocn
      use ice_state, only: tr_aero
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_couple, timer_step
      use ice_zbgc_shared, only: skl_bgc
      use ice_restart_shared, only: restart_dir, input_dir

#ifdef AusCOM 
!ars599: 27032014 add in
      use ice_timers, only: timer_into_ocn
      use ice_grid, only: t2ugrid_vector, u2tgrid_vector

      type(accessom2_type), intent(inout) :: accessom2
      type(logger_type), intent(in) :: logger

      integer (kind=int_kind) :: time_sec, itap, icpl_ai, icpl_io
      integer (kind=int_kind) :: stimestamp_ai
      integer (kind=int_kind) :: rtimestamp_io, stimestamp_io
      !receive and send timestamps (seconds)
      integer (kind=int_kind) :: imon
      logical :: first_ice_step, first_ocean_wait

      ! Keep some stats about ice_step performance and ocean wait times.
      type(simple_timer_type) :: ice_step_timer, ocean_wait_timer
      type(simple_timer_type) :: coupling_step_timer
#endif

   !--------------------------------------------------------------------
   !  initialize error code and step timer
   !--------------------------------------------------------------------

      call ice_timer_start(timer_step)   ! start timing entire run

      ! Initialise simple timers
      !call ice_step_timer%init('ice_step', logger)
      !call ocean_wait_timer%init('ocean_wait', logger)
      !call coupling_step_timer%init('coupling_step', logger)

   !--------------------------------------------------------------------
   ! timestep loop
   !--------------------------------------------------------------------

      ! Input 2-timelevel a2i data for better diurnal cycle forcing (critical to the 
      ! ice model). interpolation will be done to get the 'right' forcing in between 
      ! these two timelevels.
      ! To make this possible, we need pre-process the a2i data, let the model read 
      ! in for the time 0 (eg, 00h) at the beginning of each run segment; then we let 
      ! the coupler send the next coupling point (06h) a2i data, ie, one-cpl-interval 
      ! ahead of the real coupling point, which must be done properly in the data atm 
      ! model (matm).
      !
!XXX      call get_time0_a2i_fields('CICE_input/A2I_time0.nc') 
     
      ! restart runs need 'initial' o2i and i2o forcing fields saved at the end of 
      !    last run from ocn and ice model;
      ! initial run needs the pre-processed o2i and i2o fields.

      call get_time0_o2i_fields(trim(input_dir)//'o2i.nc')

      call get_time0_i2o_fields(trim(input_dir)//'i2o.nc')
      call get_sicemass(trim(input_dir)//'sicemass.nc')

      if (use_core_nyf_runoff) then
         stop "Don't do this"
        call get_core_runoff('INPUT/core_runoff_regrid.nc','runoff',1)
      endif

      time_sec = 0
      imon = 0

      call from_atm(time_sec)
      call update_halos_from_atm(time_sec)
      ! Shift windstress/ice-ocean stress from T onto U grid before sending into ocn
      call t2ugrid_vector(iostrsu)
      call t2ugrid_vector(iostrsv)

      DO icpl_ai = 1, num_cpl_ai   !begin I <==> A coupling iterations

! In case of CORE-IAF RUNOFF:
      if (use_core_iaf_runoff) then
         stop "Don't do this either"
        call calendar(time)
        if (imon /= month ) then
          imon = month
#if defined(DEBUG)
          print *, "use_core_iaf_runoff: icpl_ai, month, mday", icpl_ai, month, mday
          write(1001,*)'icpl_ai, month, mday = ',icpl_ai, month, mday
#endif
          call get_core_runoff('INPUT/core_runoff_regrid.nc','RUNOFF', month)
        endif
      endif

      Do icpl_io = 1, num_cpl_io   !begin I <==> O coupling iterations
        !call coupling_step_timer%start()

          if (my_task == master_task) then
            print*, 'master beginning of coupling loop', time_sec
          endif
 

        stimestamp_io = time_sec

        ! ---temp check for roughness etc.---
        if (chk_gfdl_roughness) then
           !
           !call gather_global(gwork, u_star0, master_task, distrb_info)
           !if (my_task == master_task) write(54,'(10e12.4)')gwork
           !
           call check_roughness(trim(input_dir)//'fields_roughness.nc',stimestamp_io)
        endif
        ! ----------------------------------- 

        if (debug_output) then
            call write_boundary_checksums(time_sec)
        endif

          if (my_task == master_task) then
            print*, 'master here 0', time_sec
          endif
 
        call ice_timer_start(timer_into_ocn)  ! atm/ocn coupling
        call into_ocn(stimestamp_io, 1.0)
        call ice_timer_stop(timer_into_ocn)  ! atm/ocn coupling
        !set i2o fields back to 0 for next i2o coupling period 'sum-up'
        call nullify_i2o_fluxes() 

          if (my_task == master_task) then
            print*, 'master here 1', time_sec
          endif
 
        ! Communication with atmosphere and ocean has completed. Update halos
        ! ready for ice timestep.
        call update_halos_from_ocn(time_sec)

        sss=ssso
        call new_freezingT

          if (my_task == master_task) then
            print*, 'master here 2', time_sec
          endif
 
        do itap = 1, num_ice_io    !ice time loop within each i2o cpl interval


          !put in place all (atm and ocn) 'raw' forcing fields:
          call newt_forcing_raw

          !convert the 'raw' atm forcing into that required by cice
          call get_forcing_atmo_ready

          !call ice_step_timer%start()

          if (my_task == master_task) then
            print*, 'master calling ice_step at time: ', time_sec
          endif
          call ice_step()
          if (my_task == master_task) then
            print*, 'master done calling ice_step'
          endif
          !call ice_step_timer%stop()

          istep  = istep  + 1    ! update time step counters
          istep1 = istep1 + 1
          time = time + dt       ! determine the time and date

          time_sec = time_sec + dt
          if (my_task == master_task) then
            print*, 'master done calling ice_step at time: ', time_sec
          endif

          if (my_task == master_task) then
            print*, 'master calling progress_date', time_sec
            call accessom2%progress_date(int(dt))
            print*, 'master done calling progress_date', time_sec
          endif
 
          call calendar(time)

          !initialize fluxes sent to coupler (WHY should still need do this?)
          call init_flux_atm
          call init_flux_ocn
  
        end do     ! itap

        ! ---temp check for roughness etc.---
        !if (chk_gfdl_roughness) then
        !   !
        !   call gather_global(gwork, u_star0, master_task, distrb_info)
        !   if (my_task == master_task) write(54,'(10e12.4)')gwork
        !   !
        !   call check_roughness(time_sec)
        !endif
        ! ----------------------------------- 
          if (my_task == master_task) then
            print*, 'master calling from_atm', time_sec
          endif
 
        if (icpl_io == num_cpl_io .and. icpl_ai < num_cpl_ai) then
          call from_atm(time_sec)
          call update_halos_from_atm(time_sec)
          ! Shift windstress/ice-ocean stress from T onto U grid before sending into ocn
          call t2ugrid_vector(iostrsu)
          call t2ugrid_vector(iostrsv)
        endif

          if (my_task == master_task) then
            print*, 'master done calling from_atm', time_sec
          endif
 
        rtimestamp_io = time_sec
        if (rtimestamp_io < (dt*npt)) then
          !call ocean_wait_timer%start()
          call from_ocn(rtimestamp_io)
          !call ocean_wait_timer%stop()
        endif

          if (my_task == master_task) then
            print*, 'master done calling from_ocn', time_sec
          endif
 
        !call coupling_step_timer%stop()

        print *, 'CICE: in coupling loop PE, time_sec ', my_task, time_sec
      End Do      !icpl_io

      END DO        !icpl_ai

          if (my_task == master_task) then
            print*, 'master after coupling loop'
          endif
 
        print *, 'CICE: after coupling loop PE ', my_task

      ! final update of the stimestamp_io, ie., put back the last dt_ice:
      stimestamp_io = stimestamp_io + dt

!XXX      call save_time0_a2i_fields('CICE_input/A2I_time1.nc', stimestamp_io)

      call save_time0_i2o_fields(trim(restart_dir)//'i2o.nc', stimestamp_io) 

      call save_u_star(trim(restart_dir)//'u_star.nc',stimestamp_io)    

      call save_sicemass(trim(restart_dir)//'sicemass.nc',stimestamp_io)    

        print *, 'CICE: after save_sicemass PE ', my_task

   !--------------------------------------------------------------------
   ! end of timestep loop
   !--------------------------------------------------------------------

      call ice_timer_stop(timer_step)   ! end timestepping loop timer     
      !call ice_step_timer%write_stats()
      !call ocean_wait_timer%write_stats()
      !call coupling_step_timer%write_stats()

        print *, 'CICE: cice_run finished PE ', my_task

      end subroutine CICE_Run

!=======================================================================
!
!  Calls drivers for physics components, some initialization, and output
!
!  author Elizabeth C. Hunke, LANL
!         William H. Lipscomb, LANL

      subroutine ice_step

      use ice_age, only: write_restart_age
      use ice_aerosol, only: write_restart_aero
      use ice_boundary, only: ice_HaloUpdate
      use ice_brine, only: hbrine_diags, write_restart_hbrine
      use ice_calendar, only: dt, dt_dyn, ndtd, diagfreq, write_restart, istep
      use ice_constants, only: field_loc_center, field_type_scalar
      use ice_diagnostics, only: init_mass_diags, runtime_diags
      use ice_domain, only: halo_info, nblocks
      use ice_domain_size, only: nslyr
      use ice_dyn_eap, only: write_restart_eap
      use ice_dyn_shared, only: kdyn
      use ice_firstyear, only: write_restart_FY
      use ice_flux, only: scale_factor, init_history_therm
      use ice_history, only: accum_hist
      use ice_lvl, only: write_restart_lvl
      use ice_restart, only: final_restart
      use ice_restart_driver, only: dumpfile
      use ice_meltpond_cesm, only: write_restart_pond_cesm
      use ice_meltpond_lvl, only: write_restart_pond_lvl
      use ice_meltpond_topo, only: write_restart_pond_topo
      use ice_restoring, only: restore_ice, ice_HaloRestore
      use ice_state, only: nt_qsno, trcrn, tr_iage, tr_FY, tr_lvl, &
          tr_pond_cesm, tr_pond_lvl, tr_pond_topo, tr_brine, tr_aero
      use ice_step_mod, only: prep_radiation, step_therm1, step_therm2, &
          post_thermo, step_dynamics, step_radiation
      use ice_therm_shared, only: calc_Tsfc
      use ice_timers, only: ice_timer_start, ice_timer_stop, &
          timer_diags, timer_column, timer_thermo, timer_bound, &
          timer_hist, timer_readwrite
      use ice_algae, only: bgc_diags, write_restart_bgc
      use ice_zbgc, only: init_history_bgc, biogeochemistry
      use ice_zbgc_shared, only: skl_bgc

      integer (kind=int_kind) :: &
         iblk        , & ! block index 
         k               ! dynamics supercycling index

      !-----------------------------------------------------------------
      ! restoring on grid boundaries
      !-----------------------------------------------------------------

         if (restore_ice) call ice_HaloRestore

      !-----------------------------------------------------------------
      ! initialize diagnostics
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics/history
         call init_mass_diags   ! diagnostics per timestep
         call init_history_therm
         call init_history_bgc
         call ice_timer_stop(timer_diags)   ! diagnostics/history

         call ice_timer_start(timer_column)  ! column physics
         call ice_timer_start(timer_thermo)  ! thermodynamics

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

#ifndef AusCOM
      !-----------------------------------------------------------------
      ! Scale radiation fields
      !-----------------------------------------------------------------

            if (calc_Tsfc) call prep_radiation (dt, iblk)
#endif

      !-----------------------------------------------------------------
      ! thermodynamics
      !-----------------------------------------------------------------

            call step_therm1     (dt, iblk) ! vertical thermodynamics
         enddo ! iblk

        !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks
            call biogeochemistry (dt, iblk) ! biogeochemistry
            call step_therm2     (dt, iblk) ! ice thickness distribution thermo
         enddo ! iblk
         !$OMP END PARALLEL DO

         call post_thermo (dt)             ! finalize thermo update

         call ice_timer_stop(timer_thermo) ! thermodynamics
         call ice_timer_stop(timer_column) ! column physics

      !-----------------------------------------------------------------
      ! dynamics, transport, ridging
      !-----------------------------------------------------------------

         do k = 1, ndtd
            call step_dynamics (dt_dyn, ndtd)
         enddo

      !-----------------------------------------------------------------
      ! albedo, shortwave radiation
      !-----------------------------------------------------------------

         call ice_timer_start(timer_column)  ! column physics
         call ice_timer_start(timer_thermo)  ! thermodynamics

         !$OMP PARALLEL DO PRIVATE(iblk)
         do iblk = 1, nblocks

            call step_radiation (dt, iblk)

      !-----------------------------------------------------------------
      ! get ready for coupling and the next time step
      !-----------------------------------------------------------------

            call coupling_prep (iblk)

         enddo ! iblk
         !$OMP END PARALLEL DO

         ! Calculate/merge i2o fields for each ice time step
         call get_i2o_fluxes

         ! Do time-weighted sum-up for the i2o fields
         call tavg_i2o_fluxes

         call ice_timer_start(timer_bound)
         call ice_HaloUpdate (scale_factor,     halo_info, &
                              field_loc_center, field_type_scalar)
         call ice_timer_stop(timer_bound)

         call ice_timer_stop(timer_thermo) ! thermodynamics
         call ice_timer_stop(timer_column) ! column physics

      !-----------------------------------------------------------------
      ! write data
      !-----------------------------------------------------------------

         call ice_timer_start(timer_diags)  ! diagnostics
         if (mod(istep,diagfreq) == 0) then
            call runtime_diags(dt)          ! log file
            if (skl_bgc)  call bgc_diags (dt)
            if (tr_brine) call hbrine_diags (dt)
         endif
         call ice_timer_stop(timer_diags)   ! diagnostics

         call ice_timer_start(timer_hist)   ! history
         call accum_hist (dt)               ! history file
         call ice_timer_stop(timer_hist)    ! history

         call ice_timer_start(timer_readwrite)  ! reading/writing
         if (write_restart == 1) then
            call dumpfile     ! core variables for restarting
            if (tr_iage)      call write_restart_age
            if (tr_FY)        call write_restart_FY
            if (tr_lvl)       call write_restart_lvl
            if (tr_pond_cesm) call write_restart_pond_cesm
            if (tr_pond_lvl)  call write_restart_pond_lvl
            if (tr_pond_topo) call write_restart_pond_topo
            if (tr_aero)      call write_restart_aero
            if (skl_bgc)      call write_restart_bgc  
            if (tr_brine)     call write_restart_hbrine
            if (kdyn == 2)    call write_restart_eap
            call final_restart
         endif

         call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine ice_step
    
!=======================================================================
!
! Prepare for coupling
!
! authors: Elizabeth C. Hunke, LANL

      subroutine coupling_prep (iblk)

      use ice_blocks, only: block, nx_block, ny_block
      use ice_calendar, only: dt, nstreams
      use ice_constants, only: c0, c1, puny, rhofresh
      use ice_domain_size, only: ncat
      use ice_flux, only: alvdf, alidf, alvdr, alidr, albice, albsno, &
          albpnd, albcnt, apeff_ai, coszen, fpond, fresh, &
          alvdf_ai, alidf_ai, alvdr_ai, alidr_ai, fhocn_ai, &
          fresh_ai, fsalt_ai, fsalt, &
          fswthru_ai, fhocn, fswthru, scale_factor, &
          swvdr, swidr, swvdf, swidf, Tf, Tair, Qa, strairxT, strairyt, &
          fsens, flat, fswabs, flwout, evap, Tref, Qref, faero_ocn, &
          fsurfn_f, flatn_f, scale_fluxes, frzmlt_init, frzmlt
      use ice_grid, only: tmask
      use ice_ocean, only: oceanmixed_ice, ocean_mixed_layer
      use ice_shortwave, only: alvdfn, alidfn, alvdrn, alidrn, &
                               albicen, albsnon, albpndn, apeffn
      use ice_state, only: aicen, aice, aice_init, nbtrcr
      use ice_therm_shared, only: calc_Tsfc
      use ice_timers, only: timer_couple, ice_timer_start, ice_timer_stop
      use ice_zbgc_shared, only: flux_bio, flux_bio_ai

      integer (kind=int_kind), intent(in) :: & 
         iblk            ! block index 

      ! local variables

      integer (kind=int_kind) :: & 
         n           , & ! thickness category index
         i,j         , & ! horizontal indices
         k               ! tracer index

      real (kind=dbl_kind) :: cszn ! counter for history averaging

      !-----------------------------------------------------------------
      ! Save current value of frzmlt for diagnostics.
      ! Update mixed layer with heat and radiation from ice.
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            frzmlt_init  (i,j,iblk) = frzmlt(i,j,iblk)
         enddo
         enddo

         call ice_timer_start(timer_couple)   ! atm/ocn coupling

         if (oceanmixed_ice) &
         call ocean_mixed_layer (dt,iblk) ! ocean surface fluxes and sst

#ifdef AusCOM
      if (chk_frzmlt_sst) call check_frzmlt_sst('frzmlt_sst1.nc')
#endif

      !-----------------------------------------------------------------
      ! Aggregate albedos
      !-----------------------------------------------------------------

         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = c0
            alidf(i,j,iblk) = c0
            alvdr(i,j,iblk) = c0
            alidr(i,j,iblk) = c0

            albice(i,j,iblk) = c0
            albsno(i,j,iblk) = c0
            albpnd(i,j,iblk) = c0
            apeff_ai(i,j,iblk) = c0

            ! for history averaging
            cszn = c0
            if (coszen(i,j,iblk) > puny) cszn = c1
            do n = 1, nstreams
               albcnt(i,j,iblk,n) = albcnt(i,j,iblk,n) + cszn
            enddo
         enddo
         enddo
         do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            alvdf(i,j,iblk) = alvdf(i,j,iblk) &
               + alvdfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidf(i,j,iblk) = alidf(i,j,iblk) &
               + alidfn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alvdr(i,j,iblk) = alvdr(i,j,iblk) &
               + alvdrn(i,j,n,iblk)*aicen(i,j,n,iblk)
            alidr(i,j,iblk) = alidr(i,j,iblk) &
               + alidrn(i,j,n,iblk)*aicen(i,j,n,iblk)

            if (coszen(i,j,iblk) > puny) then ! sun above horizon
            albice(i,j,iblk) = albice(i,j,iblk) &
               + albicen(i,j,n,iblk)*aicen(i,j,n,iblk)
            albsno(i,j,iblk) = albsno(i,j,iblk) &
               + albsnon(i,j,n,iblk)*aicen(i,j,n,iblk)
            albpnd(i,j,iblk) = albpnd(i,j,iblk) &
               + albpndn(i,j,n,iblk)*aicen(i,j,n,iblk)
            endif

            apeff_ai(i,j,iblk) = apeff_ai(i,j,iblk) &       ! for history
               + apeffn(i,j,n,iblk)*aicen(i,j,n,iblk)
         enddo
         enddo
         enddo

         do j = 1, ny_block
         do i = 1, nx_block

      !-----------------------------------------------------------------
      ! reduce fresh by fpond for coupling
      !-----------------------------------------------------------------

            fpond(i,j,iblk) = fpond(i,j,iblk) * rhofresh/dt
            fresh(i,j,iblk) = fresh(i,j,iblk) - fpond(i,j,iblk)

      !----------------------------------------------------------------
      ! Store grid box mean albedos and fluxes before scaling by aice
      !----------------------------------------------------------------

            alvdf_ai  (i,j,iblk) = alvdf  (i,j,iblk)
            alidf_ai  (i,j,iblk) = alidf  (i,j,iblk)
            alvdr_ai  (i,j,iblk) = alvdr  (i,j,iblk)
            alidr_ai  (i,j,iblk) = alidr  (i,j,iblk)
            fresh_ai  (i,j,iblk) = fresh  (i,j,iblk)
            fsalt_ai  (i,j,iblk) = fsalt  (i,j,iblk)
            fhocn_ai  (i,j,iblk) = fhocn  (i,j,iblk)
            fswthru_ai(i,j,iblk) = fswthru(i,j,iblk)

            if (nbtrcr > 0) then
            do k = 1, nbtrcr
              flux_bio_ai  (i,j,k,iblk) = flux_bio  (i,j,k,iblk)
            enddo
            endif

      !-----------------------------------------------------------------
      ! Save net shortwave for scaling factor in scale_factor
      !-----------------------------------------------------------------
            scale_factor(i,j,iblk) = &
                       swvdr(i,j,iblk)*(c1 - alvdr_ai(i,j,iblk)) &
                     + swvdf(i,j,iblk)*(c1 - alvdf_ai(i,j,iblk)) &
                     + swidr(i,j,iblk)*(c1 - alidr_ai(i,j,iblk)) &
                     + swidf(i,j,iblk)*(c1 - alidf_ai(i,j,iblk))

         enddo
         enddo

#ifndef AusCOM
      !B: Note this 'Scaling' operation is NOT needed for the AusCOM system
      ! 'cos the i2o fields are all properly weighted before being sent to ocn.
      ! (if mistakenly done, the ice model would send ridiculously large fluxes
      ! into ocn when aice is very small, and cause the mom4 to stop. eg- fhocn
      ! could be extremely large neg number after scaling which would cool down
      ! the 1st layer temp out of the allowed temp range......      17/03/2008)
      !
      !            ************* revisit this part later **************
      !            "per unit ice area" --- tricky --- 
      !            why hadgem3 uses this scaling?????  (check again...) 
      !

      !-----------------------------------------------------------------
      ! Divide fluxes by ice area 
      !  - the CCSM coupler assumes fluxes are per unit ice area
      !  - also needed for global budget in diagnostics
      !-----------------------------------------------------------------

         call scale_fluxes (nx_block,            ny_block,           &
                            tmask    (:,:,iblk), nbtrcr,             &
                            aice     (:,:,iblk), Tf      (:,:,iblk), &
                            Tair     (:,:,iblk), Qa      (:,:,iblk), &
                            strairxT (:,:,iblk), strairyT(:,:,iblk), &
                            fsens    (:,:,iblk), flat    (:,:,iblk), &
                            fswabs   (:,:,iblk), flwout  (:,:,iblk), &
                            evap     (:,:,iblk),                     &
                            Tref     (:,:,iblk), Qref    (:,:,iblk), &
                            fresh    (:,:,iblk), fsalt   (:,:,iblk), &
                            fhocn    (:,:,iblk), fswthru (:,:,iblk), &
                            faero_ocn(:,:,:,iblk),                   &
                            alvdr    (:,:,iblk), alidr   (:,:,iblk), &
                            alvdf    (:,:,iblk), alidf   (:,:,iblk), &
                            flux_bio(:,:,1:nbtrcr,iblk))
 
#endif

!echmod - comment this out for efficiency, if .not. calc_Tsfc
         if (.not. calc_Tsfc) then

       !---------------------------------------------------------------
       ! If surface fluxes were provided, conserve these fluxes at ice 
       ! free points by passing to ocean. 
       !---------------------------------------------------------------

            call sfcflux_to_ocn & 
                         (nx_block,              ny_block,             &
                          tmask   (:,:,iblk),    aice_init(:,:,iblk),  &
                          fsurfn_f (:,:,:,iblk), flatn_f(:,:,:,iblk),  &
                          fresh    (:,:,iblk),   fhocn    (:,:,iblk))
         endif                 
!echmod

         call ice_timer_stop(timer_couple)   ! atm/ocn coupling

      end subroutine coupling_prep

!=======================================================================
!
! If surface heat fluxes are provided to CICE instead of CICE calculating
! them internally (i.e. .not. calc_Tsfc), then these heat fluxes can 
! be provided at points which do not have ice.  (This is could be due to
! the heat fluxes being calculated on a lower resolution grid or the
! heat fluxes not recalculated at every CICE timestep.)  At ice free points, 
! conserve energy and water by passing these fluxes to the ocean.
!
! author: A. McLaren, Met Office

      subroutine sfcflux_to_ocn(nx_block,   ny_block,     &
                                tmask,      aice,         &
                                fsurfn_f,   flatn_f,      &
                                fresh,      fhocn)

      use ice_domain_size, only: ncat

      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block  ! block dimensions

      logical (kind=log_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          tmask       ! land/boundary mask, thickness (T-cell)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: &
          aice        ! initial ice concentration

      real (kind=dbl_kind), dimension(nx_block,ny_block,ncat), &
          intent(in) :: &
          fsurfn_f, & ! net surface heat flux (provided as forcing)
          flatn_f     ! latent heat flux (provided as forcing)

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: &
          fresh        , & ! fresh water flux to ocean         (kg/m2/s)
          fhocn            ! actual ocn/ice heat flx           (W/m**2)

#ifdef CICE_IN_NEMO

      ! local variables
      integer (kind=int_kind) :: &
          i, j, n    ! horizontal indices
      
      real (kind=dbl_kind)    :: &
          rLsub            ! 1/Lsub

      rLsub = c1 / Lsub

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j) .and. aice(i,j) <= puny) then
               fhocn(i,j)      = fhocn(i,j)              &
                            + fsurfn_f(i,j,n) + flatn_f(i,j,n)
               fresh(i,j)      = fresh(i,j)              &
                                 + flatn_f(i,j,n) * rLsub
            endif
         enddo   ! i
         enddo   ! j
      enddo      ! n

#endif 

      end subroutine sfcflux_to_ocn

!=======================================================================

      end module CICE_RunMod

!=======================================================================
