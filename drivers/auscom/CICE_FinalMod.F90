!  SVN:$Id: CICE_FinalMod.F90 744 2013-09-27 22:53:24Z eclare $
!=======================================================================
!
!  This module contains routines for the final exit of the CICE model,
!  including final output and clean exit from any message passing
!  environments and frameworks.
!
!  authors: Philip W. Jones, LANL
!  2006: Converted to free source form (F90) by Elizabeth Hunke
!  2008: E. Hunke moved ESMF code to its own driver

      module CICE_FinalMod

      use ice_kinds_mod

      use accessom2_mod, only : accessom2_type => accessom2
      use cpl_interface, only : coupler_termination

      implicit none
      private
      public :: CICE_Finalize
      save

!=======================================================================

      contains

!=======================================================================
!
!  This routine shuts down CICE by exiting all relevent environments.

      subroutine CICE_Finalize(accessom2)

      use ice_exit, only: end_run
      use ice_fileunits, only: nu_diag, release_all_fileunits
      use ice_restart_shared, only: runid
      use ice_timers, only: ice_timer_stop, ice_timer_print_all, timer_total
      use ice_communicate, only: my_task, master_task
      use ice_calendar, only : year_init, nyr, month, mday, hour, sec
      use ice_calendar, only : calendar, time, dt

      type(accessom2_type), intent(inout) :: accessom2

      integer, dimension(6) :: date_array
   !-------------------------------------------------------------------
   ! stop timers and print timer info
   !-------------------------------------------------------------------

      call ice_timer_stop(timer_total)        ! stop timing entire run
#ifdef AusCOM
      call ice_timer_print_all(stats=.true.) ! print timing information
#else
      call ice_timer_print_all(stats=.false.) ! print timing information
#endif

!echmod      if (nu_diag /= 6) close (nu_diag) ! diagnostic output
      call release_all_fileunits

   !-------------------------------------------------------------------
   ! write 'finished' file if needed
   !-------------------------------------------------------------------

      if (runid == 'bering') call writeout_finished_file()

   !-------------------------------------------------------------------
   ! quit MPI
   !-------------------------------------------------------------------

    ! Allow libaccessom2 to check that datetime of all models is synchronised at
    ! the end of the run.
    call calendar(time-dt)
    date_array(1) = nyr + year_init - 1
    date_array(2) = month
    date_array(3) = mday
    date_array(4) = int(sec / 3600)
    date_array(5) = int(mod(sec, 3600) / 60)
    date_array(6) = mod(sec, 60)
    call accessom2%deinit(cur_date_array=date_array)

   call coupler_termination

end subroutine CICE_Finalize

!=======================================================================
!
! Write a file indicating that this run finished cleanly.  This is
! needed only for runs on machine 'bering' (set using runid = 'bering').
!
!  author: Adrian Turner, LANL

      subroutine writeout_finished_file()
      
      use ice_restart_shared, only: restart_dir
      use ice_communicate, only: my_task, master_task

      character(len=char_len_long) :: filename

      if (my_task == master_task) then
           
         filename = trim(restart_dir)//"finished"
         open(11,file=filename)
         write(11,*) "finished"
         close(11)

      endif

      end subroutine writeout_finished_file

!=======================================================================

      end module CICE_FinalMod

!=======================================================================
