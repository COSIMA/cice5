
#if (MXBLCKS != 1)
#error The code assumes that max_blocks == 1
#endif

!============================================================================
  module cpl_interface
!============================================================================
! coupling interface between CICE and the oasis3_25 coupler (via MPI2) using 
! the PRISM System Model Interface (PSMILe).
!----------------------------------------------------------------------------

  !prism stuff
  use mpi
  use mod_prism

  !cice stuff
  use ice_kinds_mod
  use ice_communicate, only : my_task, master_task, MPI_COMM_ICE
  use ice_blocks,      only : nx_block, ny_block, nghost
  use ice_grid,        only : grid_file, kmt_file
  use ice_read_write,  only : ice_read_global_nc, ice_open_nc, ice_close_nc
  use ice_domain_size
  use ice_exit, only: abort_ice
!ars599: 26032014 add distrb
  use ice_distribution, only : nprocsX, nprocsY, distrb
  use ice_gather_scatter
  use ice_constants
  use ice_boundary, only : ice_HaloUpdate
  use ice_domain, only : distrb_info,ew_boundary_type,ns_boundary_type, halo_info

  !cpl stuff
  use cpl_parameters
  use cpl_netcdf_setup
  use cpl_arrays_setup
  use cpl_forcing_handler

  use ice_timers, only: ice_timer_start, ice_timer_stop
  use ice_timers, only: timer_from_atm_halos, timer_from_ocn_halos
  use ice_timers, only: timer_from_atm, timer_waiting_atm, timer_waiting_ocn
  use ice_timers, only: timer_from_ocn, timer_runoff_remap

!ars599: 27032014 add distrb
  !mpi stuff
  use ice_broadcast, only :  broadcast_array

  implicit none

  public :: prism_init, init_cpl, coupler_termination, get_time0_sstsss, &
            from_atm, into_ocn, from_ocn, il_commlocal
  public :: update_halos_from_ocn, update_halos_from_atm

  private

  integer(kind=int_kind), dimension(jpfldout) :: il_var_id_out ! ID for fields sent 
  integer(kind=int_kind), dimension(jpfldin)  :: il_var_id_in  ! ID for fields rcvd

  character(len=6), parameter :: cp_modnam='cicexx' ! Component model name

  integer(kind=int_kind) :: il_commlocal  ! Component internal communicator 
  integer(kind=int_kind) :: ierror
  integer(kind=int_kind) :: il_comp_id    ! Component ID
  integer(kind=int_kind) :: il_nbtotproc   ! Total number of processes
  integer(kind=int_kind) :: il_nbcplproc   ! Number of processes involved in the coupling
  integer(kind=int_kind) :: l_ilo, l_ihi, l_jlo, l_jhi !local partition

  integer(kind=int_kind) :: il_commatm, my_commatm_task

  integer :: sendsubarray, recvsubarray , resizedrecvsubarray
  integer, dimension(:), allocatable :: counts, disps
!  integer :: subgw_type

  integer(kind=int_kind) :: il_bufsize
  real(kind=dbl_kind), dimension(:,:), allocatable :: rla_array
  real(kind=dbl_kind), dimension(:),   allocatable :: rla_bufsend
  real(kind=dbl_kind), dimension(:,:), allocatable :: vwork2d

  integer(kind=int_kind), dimension(:,:), allocatable :: land_points
  integer(int_kind) :: num_land_points, num_ocean_points

  contains

!======================================================================
  subroutine prism_init
!-----------------------!

  include 'mpif.h'

  integer(kind=int_kind) :: io_size, ii, integer_byte_size, integer_io_size
  integer(kind=int_kind) :: il_real, il_bufsizebyt

  logical :: mpiflag

  character(len=12) :: chiceout
  character(len=6) :: chout

  !-----------------------------------
  ! 'define' the model global domain: 
  !-----------------------------------
  nt_cells = nx_global * ny_global

  !allocate rla_array to be used below
  allocate (rla_array(nx_global,ny_global) )

  !-------------------
  ! Initialize PSMILe.
  !-------------------

  ! Initialise MPI
  mpiflag = .FALSE.
  call MPI_Initialized (mpiflag, ierror)

  if ( .not. mpiflag ) then
    call MPI_INIT(ierror)
  endif

  call MPI_Initialized (mpiflag, ierror)

  print *, 'CICE: (prism_init) calling prism_init_comp_proto...'

  call prism_init_comp_proto (il_comp_id, cp_modnam, ierror)

  if (ierror /= PRISM_Ok) then 
      call prism_abort_proto(il_comp_id, 'cice prism_init','STOP 1')
  else
      print *, 'CICE: (prism_init) called prism_init_comp_proto !'
  endif

  !B: the following part may not be really needed(?)
  !
  ! Let's suppose the model attaches to a MPI buffer for its own use
  !
  ! ! Sophisticated way to determine buffer size needed (without "kind")
  ! ! Here one message containing rla_array

  integer_byte_size = BIT_SIZE(ii)/8
  inquire (iolength=io_size) ii
  integer_io_size = io_size
  inquire (iolength=io_size) rla_array(1,1)
  il_real = io_size/integer_io_size*integer_byte_size
  il_bufsize = nt_cells + MPI_BSEND_OVERHEAD/il_real + 1
  allocate (rla_bufsend(il_bufsize), stat = ierror)
  il_bufsizebyt = il_bufsize * il_real
  call MPI_Buffer_Attach(rla_bufsend, il_bufsizebyt, ierror)

  if (ierror /= PRISM_Ok) then
      print *, 'CICE: (prism_init) Error in MPI_Buffer_Attach.'
      call prism_abort_proto(il_comp_id, 'cice prism_init','STOP 2')
  else
      print *, 'CICE: (prism_init) MPI_Buffer_Attach ok!'
  endif
  !
  ! PSMILe attribution of local communicator.
  ! 
  !   Either MPI_COMM_WORLD if MPI2 is used,
  !   or a local communicator created by Oasis if MPI1 is used.
  !
  call prism_get_localcomm_proto(il_commlocal, ierror)
  !
  if (ierror /= PRISM_Ok) then
      print *, 'CICE: Error in prism_get_localcomm_proto'
      call prism_abort_proto(il_comp_id, 'cice prism_init','STOP 3')
  else
      print *, 'CICE: _get_localcomm_ OK! il_commlocal= ',il_commlocal
  endif
  !
  ! Inquire if model is parallel or not and open the process log file
  !
  print *, '* CICE: Entering init_cpl.....'

  print *, '* CICE4 (init_cpl) calling MPI_Comm_Size ...'
  call MPI_Comm_Size(il_commlocal, il_nbtotproc, ierror)
  print *, '* CICE4 (init_cpl) calling MPI_Comm_Rank ...'
  call MPI_Comm_Rank(il_commlocal, my_task, ierror)

  print *, '* CICE4 (init_cpl) il_commlocal, il_nbtotproc, my_task = ',&
                             il_commlocal, il_nbtotproc, my_task

  il_nbcplproc = il_nbtotproc   !multi-process coupling

  ! Open the process log file
  !if (my_task == 0 .or. ll_comparal) then
    il_out = 85 + my_task
    write(chout,'(I6.6)')il_out
    chiceout='iceout'//trim(chout)
    open(il_out,file=chiceout,form='formatted')
  
    write(il_out,*) 'Number of processes:', il_nbtotproc
    write(il_out,*) 'Local process number:', my_task
    write(il_out,*) 'Local communicator is : ',il_commlocal
    write(il_out,*) 'Grid layout: nx_global,ny_global= ',nx_global,ny_global
    write(il_out,*) 'Grid decomposition: nx_block,ny_block,max_blocks= ',&
                     nx_block,ny_block,max_blocks
  !endif

  print *, '* CICE: prism_init called OK!'  

  call prism_get_intercomm(il_commatm, 'matmxx', ierror)
  call mpi_comm_rank(il_commatm, my_commatm_task, ierror)

  end subroutine prism_init

!=======================================================================
  subroutine init_cpl

  !use mpi
  include 'mpif.h'
!--------------------!

  integer(kind=int_kind) :: jf

!  logical :: ll_comparal                   ! paralell or mono-cpl coupling
!  integer(kind=int_kind) :: il_nbtotproc   ! Total number of processes
!  integer(kind=int_kind) :: il_nbcplproc   ! Number of processes involved in the coupling
  integer(kind=int_kind), dimension(2) :: il_var_nodims ! see below
  integer(kind=int_kind), dimension(4) :: il_var_shape  ! see below
  integer(kind=int_kind) :: il_part_id     ! Local partition ID
  integer(kind=int_kind) :: il_length      ! Size of partial field for each process

!  integer,parameter :: nprocsX=3, nprocsY=2

  integer, dimension(2) :: starts,sizes,subsizes
  integer(kind=mpi_address_kind) :: start, extent
  integer (int_kind),dimension(:), allocatable :: vilo, vjlo
  real(kind=dbl_kind) :: realvalue
  integer (int_kind) :: nprocs

!-------------------------------------------------------------------------

  integer(kind=int_kind) :: ilo,ihi,jlo,jhi,iblk,i,j, n

 type (block) ::  this_block           ! block information for current block

  ! Send ice grid details to atmosphere. This is used to regrid runoff.
  call send_grid_to_atm()

!calculate partition using nprocsX and nprocsX
  write(il_out,*) 'nprocsX and nprocsY:', nprocsX, nprocsY
  l_ilo=mod(my_task,nprocsX)*nx_global/nprocsX+1
  l_ihi=l_ilo + nx_global/nprocsX -1
  l_jlo=int(my_task/nprocsX) * ny_global/nprocsY+1
  l_jhi=l_jlo+ny_global/nprocsY - 1

  write(il_out,*) '  2local partion, ilo, ihi, jlo, jhi=', l_ilo, l_ihi, l_jlo, l_jhi
  write(il_out,*) '  2partition x,y sizes:', l_ihi-l_ilo+1, l_jhi-l_jlo+1

!#ifdef _SLOW___
  nprocs = il_nbtotproc
  allocate(vilo(nprocs))
  allocate(vjlo(nprocs))

!  call mpi_type_vector(l_ihi-l_ilo+1, l_jhi-l_jlo+1, nx_global, mpi_real8, subgw_type, ierror)
!  call mpi_type_commit(subgw_type, ierror)
  call mpi_gather(l_ilo, 1, mpi_integer, vilo, 1, mpi_integer, 0, MPI_COMM_ICE, ierror)
  call broadcast_array(vilo, 0)
  call mpi_gather(l_jlo, 1, mpi_integer, vjlo, 1, mpi_integer, 0, MPI_COMM_ICE, ierror)
  call broadcast_array(vjlo, 0)

!create subarray of this rank
    sizes(1)=l_ihi-l_ilo+1; sizes(2)=l_jhi-l_jlo+1
    subsizes(1)=l_ihi-l_ilo+1; subsizes(2)=l_jhi-l_jlo+1
    starts(1)=0; starts(2)=0
    call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, &
                                mpi_real8, sendsubarray, ierror)
    call mpi_type_commit(sendsubarray,ierror)
    if (my_task == 0) then ! create recv buffer in main cpu
     sizes(1)=nx_global; sizes(2)=ny_global
     subsizes(1)=l_ihi-l_ilo+1; subsizes(2)=l_jhi-l_jlo+1
     starts(1)=0; starts(2)=0
     call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, &
                                   mpi_real8, recvsubarray, ierror)
     call mpi_type_commit(recvsubarray, ierror)
     extent = sizeof(realvalue)
     start = 0
     call mpi_type_create_resized(recvsubarray, start, extent, resizedrecvsubarray, ierror)
     call mpi_type_commit(resizedrecvsubarray,ierror)
    end if
    allocate(counts(nprocs),disps(nprocs))
    forall (jf=1:nprocs) counts(jf) = 1
    do jf=1, nprocs
      disps(jf) = ((vjlo(jf)-1)*nx_global + (vilo(jf)-1))
      !disps(n) = ((vilo(n)-1)*ny_global + (vjlo(n)-1))
    end do

    write(il_out,*) ' vilo ', vilo
    write(il_out,*) ' vjlo ', vjlo
    write(il_out,*) ' counts ', counts
    write(il_out,*) ' disps ', disps
!#endif

    !
    ! The following steps need to be done:
    ! -> by the process if cice is monoprocess;
    ! -> only by the master process, if cice is parallel and only 
    !    master process is involved in the coupling;
    ! -> by all processes, if cice is parallel and all processes 
    ! are involved in the coupling.
    
    call decomp_def (il_part_id, il_length, nt_cells, &
         my_task, il_nbcplproc, .true., il_out)

    write(il_out,*)'(init_cpl) called decomp_def, my_task, ierror = ',my_task, ierror

    !
    ! PSMILe coupling fields declaration
    !

    il_var_nodims(1) = 2 ! rank of coupling field
    il_var_nodims(2) = 1 ! number of bundles in coupling field (always 1)
    il_var_shape(1)= 1 !l_ilo ! min index for the coupling field local dimension
    il_var_shape(2)= l_ihi-l_ilo+1 ! max index for the coupling field local dim
    il_var_shape(3)= 1 !l_jlo ! min index for the coupling field local dim
    il_var_shape(4)= l_jhi-l_jlo+1 ! max index for the coupling field local dim

    !
    ! Define name (as in namcouple) and declare each field sent by ice 
    !

    !ice ==> atm
    cl_writ(1)='isst_ia'
    !ice ==> ocn
    cl_writ(n_i2a+1 )='strsu_io'
    cl_writ(n_i2a+2 )='strsv_io'
    cl_writ(n_i2a+3 )='rain_io'
    cl_writ(n_i2a+4 )='snow_io'
    cl_writ(n_i2a+5 )='stflx_io'
    cl_writ(n_i2a+6 )='htflx_io'
    cl_writ(n_i2a+7 )='swflx_io'
    cl_writ(n_i2a+8 )='qflux_io'
    cl_writ(n_i2a+9 )='shflx_io'
    cl_writ(n_i2a+10)='lwflx_io'
    cl_writ(n_i2a+11)='runof_io'
    cl_writ(n_i2a+12)='press_io'
    cl_writ(n_i2a+13)='aice_io'
    cl_writ(n_i2a+14)='melt_io'
    cl_writ(n_i2a+15)='form_io'
    !cl_writ(n_i2a+16)='co2_i1'
    !cl_writ(n_i2a+17)='wnd_i1'

    do jf=1, jpfldout
      call prism_def_var_proto (il_var_id_out(jf),cl_writ(jf), il_part_id, &
         il_var_nodims, PRISM_Out, il_var_shape, PRISM_Real, ierror)
    enddo 
    !
    ! Define name (as in namcouple) and declare each field received by ice
    !

    !atm ==> ice
    cl_read(1) ='swfld_i'
    cl_read(2) ='lwfld_i'
    cl_read(3) ='rain_i'
    cl_read(4) ='snow_i'
    cl_read(5) ='press_i'
    cl_read(6) ='runof_i'
    cl_read(7) ='tair_i'
    cl_read(8) ='qair_i'
    cl_read(9) ='uwnd_i'
    cl_read(10)='vwnd_i'
    !ocn ==> ice
    cl_read(n_a2i+1)='sst_i'
    cl_read(n_a2i+2)='sss_i'
    cl_read(n_a2i+3)='ssu_i'
    cl_read(n_a2i+4)='ssv_i'
    cl_read(n_a2i+5)='sslx_i'
    cl_read(n_a2i+6)='ssly_i'
    cl_read(n_a2i+7)='pfmice_i'
    !cl_read(n_a2i+8)='co2_oi'
    !cl_read(n_a2i+9)='co2fx_oi'
    !
    do jf=1, jpfldin
      call prism_def_var_proto (il_var_id_in(jf), cl_read(jf), il_part_id, &
         il_var_nodims, PRISM_In, il_var_shape, PRISM_Real, ierror)
    enddo 
    !
    ! 7- PSMILe end of declaration phase 
    !
    call prism_enddef_proto (ierror)

  !
  ! Allocate the 'coupling' fields (to be used) for EACH PROCESS:! 
  !

  ! fields in: (local domain)
  allocate ( tair0(nx_block,ny_block,max_blocks));  tair0(:,:,:) = 0
  allocate (swflx0(nx_block,ny_block,max_blocks)); swflx0(:,:,:) = 0
  allocate (lwflx0(nx_block,ny_block,max_blocks)); lwflx0(:,:,:) = 0
  allocate ( uwnd0(nx_block,ny_block,max_blocks));  uwnd0(:,:,:) = 0
  allocate ( vwnd0(nx_block,ny_block,max_blocks));  vwnd0(:,:,:) = 0
  allocate ( qair0(nx_block,ny_block,max_blocks));  qair0(:,:,:) = 0
  allocate ( rain0(nx_block,ny_block,max_blocks));  rain0(:,:,:) = 0
  allocate ( snow0(nx_block,ny_block,max_blocks));  snow0(:,:,:) = 0
  allocate ( runof0(nx_block,ny_block,max_blocks)); runof0(:,:,:) = 0
  allocate ( press0(nx_block,ny_block,max_blocks)); press0(:,:,:) = 0

  allocate ( runof(nx_block,ny_block,max_blocks)); runof(:,:,:) = 0
  allocate ( press(nx_block,ny_block,max_blocks)); press(:,:,:) = 0

  !
  allocate ( core_runoff(nx_block,ny_block,max_blocks));  core_runoff(:,:,:) = 0.
  !

  allocate (ssto(nx_block,ny_block,max_blocks));  ssto(:,:,:) = 0
  allocate (ssso(nx_block,ny_block,max_blocks));  ssso(:,:,:) = 0
  allocate (ssuo(nx_block,ny_block,max_blocks));  ssuo(:,:,:) = 0
  allocate (ssvo(nx_block,ny_block,max_blocks));  ssvo(:,:,:) = 0
  allocate (sslx(nx_block,ny_block,max_blocks));  sslx(:,:,:) = 0
  allocate (ssly(nx_block,ny_block,max_blocks));  ssly(:,:,:) = 0
  allocate (pfmice(nx_block,ny_block,max_blocks));  pfmice(:,:,:) = 0

  !
  allocate (iostrsu(nx_block,ny_block,max_blocks)); iostrsu(:,:,:) = 0
  allocate (iostrsv(nx_block,ny_block,max_blocks)); iostrsv(:,:,:) = 0
  allocate (iorain(nx_block,ny_block,max_blocks));  iorain(:,:,:) = 0
  allocate (iosnow(nx_block,ny_block,max_blocks));  iosnow(:,:,:) = 0
  allocate (iostflx(nx_block,ny_block,max_blocks)); iostflx(:,:,:) = 0
  allocate (iohtflx(nx_block,ny_block,max_blocks)); iohtflx(:,:,:) = 0
  allocate (ioswflx(nx_block,ny_block,max_blocks)); ioswflx(:,:,:) = 0
  allocate (ioqflux(nx_block,ny_block,max_blocks)); ioqflux(:,:,:) = 0
  allocate (iolwflx(nx_block,ny_block,max_blocks)); iolwflx(:,:,:) = 0
  allocate (ioshflx(nx_block,ny_block,max_blocks)); ioshflx(:,:,:) = 0
  allocate (iorunof(nx_block,ny_block,max_blocks)); iorunof(:,:,:) = 0
  allocate (iopress(nx_block,ny_block,max_blocks)); iopress(:,:,:) = 0
  allocate (ioaice (nx_block,ny_block,max_blocks)); ioaice(:,:,:) = 0
  !!!
  allocate (iomelt (nx_block,ny_block,max_blocks)); iomelt(:,:,:) = 0
  allocate (ioform (nx_block,ny_block,max_blocks)); ioform(:,:,:) = 0

  allocate (tiostrsu(nx_block,ny_block,max_blocks)); tiostrsu(:,:,:) = 0
  allocate (tiostrsv(nx_block,ny_block,max_blocks)); tiostrsv(:,:,:) = 0
  allocate (tiorain(nx_block,ny_block,max_blocks));  tiorain(:,:,:) = 0
  allocate (tiosnow(nx_block,ny_block,max_blocks));  tiosnow(:,:,:) = 0
  allocate (tiostflx(nx_block,ny_block,max_blocks)); tiostflx(:,:,:) = 0
  allocate (tiohtflx(nx_block,ny_block,max_blocks)); tiohtflx(:,:,:) = 0
  allocate (tioswflx(nx_block,ny_block,max_blocks)); tioswflx(:,:,:) = 0
  allocate (tioqflux(nx_block,ny_block,max_blocks)); tioqflux(:,:,:) = 0
  allocate (tiolwflx(nx_block,ny_block,max_blocks)); tiolwflx(:,:,:) = 0
  allocate (tioshflx(nx_block,ny_block,max_blocks)); tioshflx(:,:,:) = 0
  allocate (tiorunof(nx_block,ny_block,max_blocks)); tiorunof(:,:,:) = 0
  allocate (tiopress(nx_block,ny_block,max_blocks)); tiopress(:,:,:) = 0
  allocate (tioaice(nx_block,ny_block,max_blocks));  tioaice(:,:,:) = 0
  !!!
  allocate (tiomelt(nx_block,ny_block,max_blocks));  tiomelt(:,:,:) = 0
  allocate (tioform(nx_block,ny_block,max_blocks));  tioform(:,:,:) = 0

  allocate (vwork(nx_block,ny_block,max_blocks)); vwork(:,:,:) = 0
  allocate (gwork(nx_global,ny_global)); gwork(:,:) = 0

  !
  allocate (sicemass(nx_block,ny_block,max_blocks)); sicemass(:,:,:) = 0.
  allocate (u_star0(nx_block,ny_block,max_blocks)); u_star0(:,:,:) = 0.
  allocate (rough_mom0(nx_block,ny_block,max_blocks)); rough_mom0(:,:,:) = 0.
  allocate (rough_heat0(nx_block,ny_block,max_blocks)); rough_heat0(:,:,:) = 0.
  allocate (rough_moist0(nx_block,ny_block,max_blocks)); rough_moist0(:,:,:) = 0.
  !
  allocate (vwork2d(l_ilo:l_ihi, l_jlo:l_jhi)); vwork2d(:,:) = 0.

end subroutine init_cpl

subroutine send_grid_to_atm()

  integer(kind=int_kind) :: tag, buf_int(2)
  real(kind=dbl_kind), dimension(:), allocatable :: buf_real

  integer :: fid
  real(kind=dbl_kind), dimension(:, :), allocatable :: tlat_global, tlon_global
  real(kind=dbl_kind), dimension(:, :), allocatable :: mask_global

  if (my_task == master_task) then
    write(il_out,*) 'Before send_grid_to_atm'

    call ice_open_nc(grid_file, fid)

    allocate(tlat_global(nx_global, ny_global))
    allocate(tlon_global(nx_global, ny_global))
    call ice_open_nc(grid_file, fid)
    call ice_read_global_nc(fid , 1, 'tlat' , tlat_global, .false.)
    call ice_read_global_nc(fid , 1, 'tlon' , tlon_global, .false.)
    call ice_close_nc(fid)

    allocate(mask_global(nx_global, ny_global))
    call ice_open_nc(kmt_file, fid)
    call ice_read_global_nc(fid , 1, 'kmt' , mask_global, .false.)
    call ice_close_nc(fid)

    write(il_out,*) 'After reading files'

    ! Send my details to the atm.
    tag = 0
    buf_int(1) = nx_global
    buf_int(2) = ny_global
    call MPI_send(buf_int, 2, MPI_INTEGER, 0, tag, il_commatm, ierror)

    allocate(buf_real(nx_global*ny_global))
    buf_real(:) = reshape(tlat_global(:, :), (/ size(tlat_global) /))
    call MPI_send(buf_real, nx_global*ny_global, MPI_DOUBLE, 0, tag, &
                  il_commatm, ierror)

    buf_real(:) = reshape(tlon_global(:, :), (/ size(tlon_global) /))
    call MPI_send(buf_real, nx_global*ny_global, MPI_DOUBLE, 0, tag, &
                  il_commatm, ierror)

    buf_real(:) = reshape(mask_global(:, :), (/ size(mask_global) /))
    call MPI_send(buf_real, nx_global*ny_global, MPI_DOUBLE, 0, tag, &
                  il_commatm, ierror)

    deallocate(buf_real)
    deallocate(tlat_global)
    deallocate(tlon_global)
    deallocate(mask_global)

    write(il_out,*) 'After send_grid_to_atm'
  endif

end subroutine send_grid_to_atm

subroutine from_atm(isteps)

  implicit none

  integer(kind=int_kind), intent(in) :: isteps

  integer(kind=int_kind) :: tag, request, info, i, j, n
  integer(kind=int_kind) :: buf(1)
  integer :: isc, iec, jsc, jec
  real(kind=dbl_kind) :: total_runoff

#if defined(DEBUG)
    write(il_out,*) '(from_atm) receiving coupling fields at rtime= ', isteps
#endif

  isc = 1+nghost
  iec = nx_block-nghost
  jsc = isc
  jec = ny_block-nghost

  call ice_timer_start(timer_from_atm)

  ! The return value 'info' is not checked, oasis does not return errors, it
  ! will abort if there is any problem.
  call ice_timer_start(timer_waiting_atm)
  call prism_get_proto(il_var_id_in(1), isteps, swflx0(isc:iec, jsc:jec, 1), info)
  call ice_timer_stop(timer_waiting_atm)
  call prism_get_proto(il_var_id_in(2), isteps, lwflx0(isc:iec, jsc:jec, 1), info)
  call prism_get_proto(il_var_id_in(3), isteps, rain0(isc:iec, jsc:jec, 1), info)
  call prism_get_proto(il_var_id_in(4), isteps, snow0(isc:iec, jsc:jec, 1), info)
  call prism_get_proto(il_var_id_in(5), isteps, press0(isc:iec, jsc:jec, 1), info)
  call prism_get_proto(il_var_id_in(6), isteps, runof0(isc:iec, jsc:jec, 1), info)
  call prism_get_proto(il_var_id_in(7), isteps, tair0(isc:iec, jsc:jec, 1), info)
  call prism_get_proto(il_var_id_in(8), isteps, qair0(isc:iec, jsc:jec, 1), info)
  call prism_get_proto(il_var_id_in(9), isteps, uwnd0(isc:iec, jsc:jec, 1), info)
  call prism_get_proto(il_var_id_in(10), isteps, vwnd0(isc:iec, jsc:jec, 1), info)

  ! need do t-grid to u-grid shift for vectors since all coupling occur on
  ! t-grid points: <==No! actually CICE requires the input wind on T grid! 
  ! (see comment in code ice_flux.F)
  !call t2ugrid(uwnd1)
  !call t2ugrid(vwnd1)
  ! ...and, as we use direct o-i communication and o-i share the same grid, 
  ! no need for any t2u and/or u2t shift before/after i-o coupling!

  if ( chk_a2i_fields ) then
    call check_a2i_fields('fields_a2i_in_ice.nc',isteps)
  endif

  ! Allow atm to progress. It is waiting on a receive.
  if (my_commatm_task == 0) then
    request = MPI_REQUEST_NULL
    tag = 0
    call MPI_Isend(buf, 1, MPI_INTEGER, 0, tag, il_commatm, request, ierror)
  endif

  call ice_timer_stop(timer_from_atm)

end subroutine from_atm

!=======================================================================
  subroutine from_ocn(isteps)
!-------------------------------------------!  

  integer(kind=int_kind), intent(in) :: isteps
 
  integer(kind=int_kind) :: jf, field_type

  call ice_timer_start(timer_from_ocn)

#if defined(DEBUG)
    write(il_out,*) '(from_ocn) receiving coupling fields at rtime: ', isteps
#endif

  do jf = n_a2i+1, jpfldin          !no 11-17 from ocn

      if (jf == n_a2i+1) then
        call ice_timer_start(timer_waiting_ocn)
      elseif (jf == n_a2i+2) then
        call ice_timer_stop(timer_waiting_ocn)
      endif

      !jf-th field in
#if defined(DEBUG)
      write(il_out,*) '*** receiving coupling fields No. ', jf, cl_read(jf)
#endif

        call prism_get_proto (il_var_id_in(jf), isteps, vwork2d(l_ilo:l_ihi, l_jlo:l_jhi), ierror)
      if ( ierror /= PRISM_Ok .and. ierror < PRISM_Recvd) then
        write(il_out,*) 'Err in _get_ sst at time with error: ', isteps, ierror
        call prism_abort_proto(il_comp_id, 'cice from_atm','stop 1')
      else
#if defined(DEBUG)
        write(il_out,*)'(from_ocn) rcvd at time with err: ',cl_read(jf),isteps,ierror
#endif
      endif

    ! Copy over non-ghost part of coupled field.
    select case (jf)
        case (n_a2i+1)
            ssto(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1) = vwork2d
        case (n_a2i+2)
            ssso(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1) = vwork2d
        case (n_a2i+3)
            ssuo(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1) = vwork2d
        case (n_a2i+4)
            ssvo(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1) = vwork2d
        case (n_a2i+5)
            sslx(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1) = vwork2d
        case (n_a2i+6)
            ssly(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1) = vwork2d
        case (n_a2i+7)
            pfmice(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1) = vwork2d
        case default
            stop "Error: invalid case in subroutine from_ocn()"
    end select

  enddo

  if (chk_o2i_fields) then
    call check_o2i_fields('fields_o2i_in_ice.nc',isteps)
  endif

  call ice_timer_stop(timer_from_ocn)  ! atm/ocn coupling

  end subroutine from_ocn

!=======================================================================
  subroutine into_ocn(isteps, scale)
!
! Note dummy 'scale', if /= 1 (then must be 1/coef_ic), is used here for the very 
! first-time-step-of-exp i2o fluxes scaling-up, because routine 'tavg_i2o_fluxes' 
! has scaled-down the current step i2o fluxes (calculated in get_i2o_fluxes) by 
! * coef_ic.  
! 
  integer(kind=int_kind), intent(in) :: isteps
  real, intent(in) :: scale             !only 1 or 1/coef_ic allowed! 
  integer(kind=int_kind) :: jf

#if defined(DEBUG)
  write(il_out,*) '(into_ocn) sending coupling fields at stime= ', isteps

  write(il_out,*)'chk iostrsu:', isteps, minval(iostrsu), maxval(iostrsu), sum(iostrsu)
  write(il_out,*)'chk iostrsv:', isteps, minval(iostrsv), maxval(iostrsv), sum(iostrsv)
  write(il_out,*)'chk iorain:', isteps, minval(iorain), maxval(iorain), sum(iorain)
  write(il_out,*)'chk iosnow:', isteps, minval(iosnow), maxval(iosnow), sum(iosnow)
  write(il_out,*)'chk iostflx:', isteps, minval(iostflx), maxval(iostflx), sum(iostflx)
  write(il_out,*)'chk iohtflx:', isteps, minval(iohtflx), maxval(iohtflx), sum(iohtflx)
  write(il_out,*)'chk ioswflx:', isteps, minval(ioswflx), maxval(ioswflx), sum(ioswflx)
  write(il_out,*)'chk ioqflux:', isteps, minval(ioqflux), maxval(ioqflux), sum(ioqflux)
  write(il_out,*)'chk ioshflx:', isteps, minval(ioshflx), maxval(ioshflx), sum(ioshflx)
  write(il_out,*)'chk iolwflx:', isteps, minval(iolwflx), maxval(iolwflx), sum(iolwflx)
  write(il_out,*)'chk iopress:', isteps, minval(iopress), maxval(iopress), sum(iopress)
  write(il_out,*)'chk ioaice:', isteps, minval(ioaice), maxval(ioaice), sum(ioaice)
  write(il_out,*)'chk iomelt:', isteps, minval(iomelt), maxval(iomelt), sum(iomelt)
  write(il_out,*)'chk ioform:', isteps, minval(ioform), maxval(ioform), sum(ioform)
#endif

  do jf = n_i2a+1, jpfldout       !no 2-14 are for the ocn

    if (jf == n_i2a+1 ) then
        vwork2d(:,:) = scale * iostrsu(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+2 ) then
        vwork2d(:,:) = scale * iostrsv(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+3 ) then
        vwork2d(:,:) = scale * iorain(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+4 ) then
        vwork2d(:,:) = scale * iosnow(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+5 ) then
        vwork2d(:,:) = scale * iostflx(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+6 ) then
         vwork2d(:,:) = scale * iohtflx(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+7 ) then
         vwork2d(:,:) = scale * ioswflx(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+8 ) then
        vwork2d(:,:) = scale * ioqflux(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+9 ) then
        vwork2d(:,:) = scale * ioshflx(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+10 ) then
        vwork2d(:,:) = scale * iolwflx(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+11) then 
       if ( use_core_nyf_runoff .or. use_core_iaf_runoff ) then 
         vwork2d(:,:) = core_runoff(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
       else 
         vwork2d(:,:) = scale * iorunof(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
       endif 
    elseif (jf == n_i2a+12) then
        vwork2d(:,:) = scale * iopress(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+13) then
        vwork2d(:,:) = scale * ioaice(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+14) then
        vwork2d(:,:) = scale * iomelt(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    elseif (jf == n_i2a+15) then
        vwork2d(:,:) = scale * ioform(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)
    endif

    !vwork2d(l_ilo:l_ihi, l_jlo:l_jhi) = vwork(1+nghost:nx_block-nghost,1+nghost:ny_block-nghost, 1)

#if defined(DEBUG)
      write(il_out,*) '*** sending coupling field No. ', jf, cl_writ(jf)
      write(il_out,*) 'chk: vwork2d', isteps, minval(vwork2d), maxval(vwork2d), sum(vwork2d)
#endif
      call prism_put_proto(il_var_id_out(jf), isteps, vwork2d(l_ilo:l_ihi, l_jlo:l_jhi), ierror)
      if ( ierror /= PRISM_Ok .and. ierror < PRISM_Sent) then
        write(il_out,*) '(into_ocn) Err in _put_ ', cl_writ(jf), isteps, ierror
        call prism_abort_proto(il_comp_id, 'cice into_ocn','STOP 1') 
      else
#if defined(DEBUG)
        write(il_out,*)'(into_ocn) sent: ', cl_writ(jf), isteps, ierror
#endif
      endif

  enddo     !jf = 6, jpfldout

  if (chk_i2o_fields) then
    call check_i2o_fields('fields_i2o_in_ice.nc',isteps, scale)
  endif

  end subroutine into_ocn

subroutine update_halos_from_ocn(time)

  integer(kind=int_kind), intent(in) :: time

  ! Fields from ocean.
  call ice_timer_start(timer_from_ocn_halos)
  call ice_HaloUpdate(ssto, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(ssso, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(ssuo, halo_info, field_loc_center, field_type_vector)
  call ice_HaloUpdate(ssvo, halo_info, field_loc_center, field_type_vector)
  call ice_HaloUpdate(sslx, halo_info, field_loc_center, field_type_vector)
  call ice_HaloUpdate(ssly, halo_info, field_loc_center, field_type_vector)
  call ice_HaloUpdate(pfmice, halo_info, field_loc_center, field_type_scalar)
  call ice_timer_stop(timer_from_ocn_halos)

#if defined(DEBUG)
  write(il_out,*)'chk ssto:', time, minval(ssto), maxval(ssto), sum(ssto)
  write(il_out,*)'chk ssso:', time, minval(ssso), maxval(ssso), sum(ssso)
  write(il_out,*)'chk ssuo:', time, minval(ssuo), maxval(ssuo), sum(ssuo)
  write(il_out,*)'chk ssvo:', time, minval(ssvo), maxval(ssvo), sum(ssvo)
  write(il_out,*)'chk sslx:', time, minval(sslx), maxval(sslx), sum(sslx)
  write(il_out,*)'chk ssly:', time, minval(ssly), maxval(ssly), sum(ssly)
  write(il_out,*)'chk pfmice:', time, minval(pfmice), maxval(pfmice), sum(pfmice)
#endif

end subroutine

subroutine update_halos_from_atm(time)

  integer(kind=int_kind), intent(in) :: time

  call ice_timer_start(timer_from_atm_halos)
  call ice_HaloUpdate(swflx0, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(lwflx0, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(rain0, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(snow0, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(press0, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(runof0, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(tair0, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(qair0, halo_info, field_loc_center, field_type_scalar)
  call ice_HaloUpdate(uwnd0, halo_info, field_loc_center, field_type_vector)
  call ice_HaloUpdate(vwnd0, halo_info, field_loc_center, field_type_vector)
  call ice_timer_stop(timer_from_atm_halos)

#if defined(DEBUG)
  write(il_out,*)'chk swflx0:', time, minval(swflx0), maxval(swflx0), sum(swflx0)
  write(il_out,*)'chk lwflx0:', time, minval(lwflx0), maxval(lwflx0), sum(lwflx0)
  write(il_out,*)'chk rain0:', time, minval(rain0), maxval(rain0), sum(rain0)
  write(il_out,*)'chk snow0:', time, minval(snow0), maxval(snow0), sum(snow0)
  write(il_out,*)'chk press0:', time, minval(press0), maxval(press0), sum(press0)
  write(il_out,*)'chk runof0:', time, minval(runof0), maxval(runof0), sum(runof0)
  write(il_out,*)'chk tair0:', time, minval(tair0), maxval(tair0), sum(tair0)
  write(il_out,*)'chk qair0:', time, minval(qair0), maxval(qair0), sum(qair0)
  write(il_out,*)'chk uwnd0:', time, minval(uwnd0), maxval(uwnd0), sum(uwnd0)
  write(il_out,*)'chk vwnd0:', time, minval(vwnd0), maxval(vwnd0), sum(vwnd0)
#endif

end subroutine update_halos_from_atm

!=======================================================================
  subroutine coupler_termination
!-------------------------------!
  !
  ! Detach from MPI buffer
  !
  call MPI_Buffer_Detach(rla_bufsend, il_bufsize, ierror) 
  write(il_out,*) 'il_bufsize =', il_bufsize
  deallocate (rla_bufsend)

  deallocate (tair0, swflx0, lwflx0, uwnd0, vwnd0, qair0, rain0, snow0, runof0, press0)
  deallocate (runof, press)
  deallocate (core_runoff)
  deallocate (ssto, ssso, ssuo, ssvo, sslx, ssly, pfmice)
  deallocate (iostrsu, iostrsv, iorain, iosnow, iostflx, iohtflx, ioswflx, &
              ioqflux, iolwflx, ioshflx, iorunof, iopress)
  deallocate (tiostrsu, tiostrsv, tiorain, tiosnow, tiostflx, tiohtflx, tioswflx, &
              tioqflux, tiolwflx, tioshflx, tiorunof, tiopress) 
  deallocate (iomelt, ioform, tiomelt, tioform)
  deallocate (gwork, vwork, sicemass)
  !  
  ! PSMILe termination 
  !   
  call prism_terminate_proto (ierror)
  if (ierror /= PRISM_Ok) then
    if (my_task == 0) then
      write (il_out,*) 'An error occured in prism_terminate = ', ierror
    endif
  else 
    if (my_task == 0) then
      write(il_out,*) '(main) calling prism_terminate_proto done!'
      write(il_out,*) '==================*** END ***================='
    endif
  endif
  ! 
  print *
  print *, '********** End of CICE **********'
  print *

  close(il_out)
  call MPI_Finalize (ierror)

  end subroutine coupler_termination

!=======================================================================
  subroutine decomp_def(id_part_id, id_length, id_imjm, &
   id_rank, id_nbcplproc, ld_comparal, ld_mparout)
!-------------------------------------------------------!
  !
  !use mod_prism_proto
  !use mod_prism_def_partition_proto

  implicit none

  integer(kind=int_kind), dimension(:), allocatable :: il_paral ! Decomposition for each proc
  integer(kind=int_kind) :: ig_nsegments  ! Number of segments of process decomposition 
  integer(kind=int_kind) :: ig_parsize    ! Size of array decomposition
  integer(kind=int_kind) :: id_nbcplproc  ! Number of processes involved in the coupling
  integer(kind=int_kind) :: id_part_id    ! Local partition ID
  integer(kind=int_kind) :: id_imjm       ! Total grid dimension, ib, ierror, my_task
  integer(kind=int_kind) :: id_length     ! Size of partial field for each process
  integer(kind=int_kind) :: id_rank       ! Rank of process
  integer(kind=int_kind) :: ld_mparout    ! Unit of log file
  logical :: ld_comparal
  integer(kind=int_kind) :: ib, ierror
  character(len=80), parameter :: cdec='BOX'
  !
  integer(kind=int_kind) :: ilo, ihi, jlo, jhi
  !
  !
  ! Refer to oasis/psmile/prism/modules/mod_prism_proto.F90 for integer(kind=int_kind) value
  ! of clim_xxxx parameters
  !
  if ( .not. ld_comparal .and. id_rank == 0) then
      ! Monoprocess model, or parallel model with only master process involved 
      ! in coupling: the entire field will be exchanged by the process. 
      ig_nsegments = 1
      ig_parsize = 3
      allocate(il_paral(ig_parsize))
      !
      il_paral ( clim_strategy ) = clim_serial
      il_paral ( clim_offset   ) = 0
      il_paral ( clim_length   ) = id_imjm
      id_length = id_imjm
      !
      call prism_def_partition_proto (id_part_id, il_paral, ierror)
      deallocate(il_paral)
      !
  else
      ! Parallel atm with all process involved in the coupling
      !
      if (cdec == 'APPLE') then
          ! Each process is responsible for a part of field defined by
          ! the number of grid points and the offset of the first point
          !
          write (ld_mparout,*) 'APPLE partitioning'
          ig_nsegments = 1
          ig_parsize = 3
          allocate(il_paral(ig_parsize))
          write(ld_mparout,*)'ig_parsize',ig_parsize
          !
          if (id_rank .LT. (id_nbcplproc-1)) then
              il_paral ( clim_strategy ) = clim_apple
              il_paral ( clim_length   ) = id_imjm/id_nbcplproc
              il_paral ( clim_offset   ) = id_rank*(id_imjm/id_nbcplproc)
          else
              il_paral ( clim_strategy ) = clim_apple
              il_paral ( clim_length   ) = id_imjm-(id_rank*(id_imjm/id_nbcplproc))
              il_paral ( clim_offset   ) = id_rank*(id_imjm/id_nbcplproc)
          endif
          id_length = il_paral(clim_length) 
          !
          call prism_def_partition_proto (id_part_id, il_paral, ierror)
          deallocate(il_paral)
          !
      else if (cdec == 'BOX') then
          !B: CICE uses a kind of Cartisian decomposition which actually may NOT
          !   be simply taken as "BOX" decomposition described here !!!
          !   (there is an issue associated with the 'halo' boundary for each
          !    segment and may NOT be treated as what we do below! 
          !    It needs further consideration to make this work correctly 
          !    for 'paralell coupling' if really needed in the future ...)
          !  
          ! Each process is responsible for a rectangular box 
          !
          write (ld_mparout,*) 'BOX partitioning'
          ig_parsize = 5
          allocate(il_paral(ig_parsize))
          write(ld_mparout,*)'ig_parsize',ig_parsize
          
          il_paral ( clim_strategy ) = clim_Box
          il_paral ( clim_offset   ) = (l_ilo-1)+(l_jlo-1)*nx_global
          il_paral ( clim_SizeX    ) = l_ihi-l_ilo+1
          il_paral ( clim_SizeY    ) = l_jhi-l_jlo+1
          il_paral ( clim_LdX      ) = nx_global

          id_length = il_paral(clim_sizeX) * il_paral(clim_sizeY)

          write(ld_mparout,*)'il_paral: ',il_paral
          
          call prism_def_partition_proto (id_part_id, il_paral, ierror)
          deallocate(il_paral)
          !
      else if (cdec == 'ORANGE') then
          !B: NOT FOR COMMON USE!
          ! Each process is responsible for arbitrarily distributed
          ! pieces of the field (here two segments by process)
          !
          write (ld_mparout,*) 'ORANGE partitioning'
          ig_nsegments = 2
          ig_parsize = 2 * ig_nsegments + 2
          write(ld_mparout,*)'ig_parsize',ig_parsize
          allocate(il_paral(ig_parsize))
          !
          il_paral ( clim_strategy   ) = clim_orange
          il_paral ( clim_segments   ) = 2
          il_paral ( clim_segments+1 ) = id_rank*768
          il_paral ( clim_segments+2 ) = 768
          il_paral ( clim_segments+3 ) = (id_rank+3)*768
          il_paral ( clim_segments+4 ) = 768
          id_length = 0
          do ib=1,2*il_paral(clim_segments) 
            if (mod(ib,2).eq.0) then
                id_length = id_length + il_paral(clim_segments+ib)
            endif
          enddo
          !
          call prism_def_partition_proto (id_part_id, il_paral, ierror)
          deallocate(il_paral)
          !
      else
          write (ld_mparout,*) 'incorrect decomposition '
      endif
  endif

  end subroutine decomp_def

  end module cpl_interface
