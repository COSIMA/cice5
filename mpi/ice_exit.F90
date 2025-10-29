!  SVN:$Id: ice_exit.F90 925 2015-03-04 00:34:27Z eclare $
!=======================================================================
!
! Exit the model. 
! authors William H. Lipscomb (LANL)
!         Elizabeth C. Hunke (LANL)
! 2006 ECH: separated serial and mpi functionality

      module ice_exit

      use ice_kinds_mod

      implicit none
      public

!=======================================================================

      contains

!=======================================================================

      subroutine abort_ice(error_message)

!  This routine aborts the ice model and prints an error message.

#if (defined CCSMCOUPLED)
      use shr_sys_mod
#else
      use ice_fileunits, only: nu_diag, ice_stderr, ice_stdout, &
                               flush_fileunit
      include 'mpif.h'   ! MPI Fortran include file
#endif

      character (len=*), intent(in) :: error_message

      ! local variables

#ifndef CCSMCOUPLED
      ! MPI error flag, 
      integer (int_kind) :: errorcode
      ! MPI return value
      integer (int_kind) :: ierr
#endif

#if (defined CCSMCOUPLED)
      call flush_fileunit(nu_diag)
      write (nu_diag,*) error_message
      call flush_fileunit(nu_diag)
      call shr_sys_abort(error_message)
#else
      call flush_fileunit(nu_diag)

      write (ice_stdout,*) error_message
      call flush_fileunit(ice_stdout)
      write (ice_stderr,*) error_message
      call flush_fileunit(ice_stderr)

      error_code = 1, !default to non-zero error
      call MPI_ABORT(MPI_COMM_WORLD, error_code, ierr)
      stop
#endif

      end subroutine abort_ice

!=======================================================================

      subroutine end_run

! Ends run by calling MPI_FINALIZE.

      integer (int_kind) :: ierr ! MPI error flag

      call MPI_FINALIZE(ierr)

      end subroutine end_run

!=======================================================================

      end module ice_exit

!=======================================================================
