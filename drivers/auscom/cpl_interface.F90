!============================================================================
  module cpl_interface
!============================================================================
! coupling interface between CICE and the oasis3_25 coupler (via MPI2) using
! the PRISM System Model Interface (PSMILe).
!----------------------------------------------------------------------------

  !prism stuff
  use mpi
  use mod_prism
  use mod_oasis, only : oasis_def_partition, oasis_def_var
  use mod_oasis, only : oasis_put, oasis_get, oasis_enddef

  !cice stuff
  use ice_kinds_mod
  use ice_communicate, only : my_task, master_task, MPI_COMM_ICE
  use ice_blocks,      only : nx_block, ny_block, nghost
  use ice_grid,        only : grid_file, kmt_file
  use ice_read_write,  only : ice_read_global_nc, ice_open_nc, ice_close_nc
  use ice_domain_size
  use ice_exit, only: abort_ice
  use ice_gather_scatter
  use ice_constants
  use ice_boundary, only : ice_HaloUpdate
  use ice_domain, only: distribution_type, distrb_info
  use ice_domain, only: nblocks, blocks_ice
  use ice_domain, only: ew_boundary_type, ns_boundary_type, halo_info

  !cpl stuff
  use cpl_parameters
  use cpl_netcdf_setup
  use cpl_arrays_setup
  use cpl_forcing_handler

  use ice_timers, only: ice_timer_start, ice_timer_stop
  use ice_timers, only: timer_from_atm_halos, timer_from_ocn_halos
  use ice_timers, only: timer_from_atm, timer_waiting_atm, timer_waiting_ocn
  use ice_timers, only: timer_from_ocn, timer_runoff_remap

  use coupler_mod, only: coupler_type => coupler

  implicit none

  public :: prism_init, init_cpl, coupler_termination, get_time0_sstsss, &
            from_atm, into_ocn, from_ocn, il_commlocal
  public :: update_halos_from_ocn, update_halos_from_atm
  public :: write_boundary_checksums
  public :: coupler
  public :: segment

  private

  ! Nomenclature:
  ! 'partition': a specification of all grid points (in global space) that this
  ! rank is responsible for. Made up from a list of segments.
  ! 'segment': a single element of the partition. A segment is defined by a
  ! global offset and a size. In this case the size is fixed so we don't bother
  ! keeping track of it.

  ! Segment specifies a list of points using a global offset and size. It also
  ! tracks the CICE block that these points come from.
  type segment
    integer :: global_offset
    integer :: block_index
  endtype

  ! List of all segments that make up the partition. This is kept in sorted
  ! order according to ascending global_offset of segments.
  type(segment), dimension(:), allocatable :: part_def

  integer(kind=int_kind), dimension(jpfldout) :: il_var_id_out ! ID for fields sent
  integer(kind=int_kind), dimension(jpfldin)  :: il_var_id_in  ! ID for fields rcvd

  character(len=6), parameter :: cp_modnam='cicexx' ! Component model name
  integer, parameter :: ORANGE = 3

  integer(kind=int_kind) :: il_commlocal  ! Component internal communicator
  integer(kind=int_kind) :: ierror
  integer(kind=int_kind) :: il_comp_id    ! Component ID
  integer(kind=int_kind) :: il_nbtotproc   ! Total number of processes
  integer(kind=int_kind) :: il_nbcplproc   ! Number of processes involved in the coupling
  integer(kind=int_kind) :: l_ilo, l_ihi, l_jlo, l_jhi !local partition

  real(kind=dbl_kind), dimension(:,:), allocatable :: vwork2d

  type(coupler_type) :: coupler

  contains

!======================================================================
  subroutine prism_init(accessom2_config_dir)

    character(len=*), intent(in) :: accessom2_config_dir

  logical :: mpiflag

  character(len=12) :: chiceout
  character(len=6) :: chout

  ! NOTE: This function can probably be replaced by coupler%init_begin, but
  !       let's move slowly for now.

  !-----------------------------------
  ! 'define' the model global domain:
  !-----------------------------------
  nt_cells = nx_global * ny_global

  !-------------------
  ! Initialize PSMILe.
  !-------------------

  call coupler%init_begin('cicexx', config_dir=trim(accessom2_config_dir))

  il_commlocal = coupler%localcomm

  !
  ! Inquire if model is parallel or not and open the process log file
  !

  call MPI_Comm_Size(il_commlocal, il_nbtotproc, ierror)
  my_task = coupler%my_local_pe

  il_nbcplproc = il_nbtotproc   !multi-process coupling

  ! Open the process log file
#if defined(DEBUG)
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
#endif

end subroutine prism_init

!> Sort segments in ascending order using an exchange sort
subroutine sort_segments(seg_list)
    type(segment), dimension(:), intent(inout) :: seg_list

    integer :: i, j
    type(segment) :: tmp

    do i=1, size(seg_list) - 1
        do j=i+1, size(seg_list)
            if (seg_list(i).global_offset > seg_list(j).global_offset) then
                tmp = seg_list(i)
                seg_list(i) = seg_list(j)
                seg_list(j) = tmp
            endif
        enddo
    enddo

endsubroutine sort_segments

subroutine init_cpl(runtime_seconds, coupling_field_timesteps)

    integer, intent(in) :: runtime_seconds
    integer, dimension(:), intent(in) :: coupling_field_timesteps

    integer(kind=int_kind) :: jf

    integer(kind=int_kind), dimension(:), allocatable :: oasis_part_def
    integer :: part_id, part_idx

    integer(kind=int_kind), dimension(2) :: il_var_nodims ! see below
    integer(kind=int_kind), dimension(4) :: il_var_shape  ! see below

    integer(kind=int_kind) :: ilo,ihi,jlo,jhi,iblk,i,j,n,m

    integer(2), external :: cmp_function

    integer :: err
    type(block) :: this_block

    ! Send ice grid details to atmosphere. This is used to regrid runoff.
    call send_grid_to_atm()

    ! Define oasis partition and variables using orange partition. This is
    ! fairly general so other partition types should not be needed.

    ! Orange partitioning allows us to define the partition for a particular PE
    ! as a collection of contiguous segments. Each segment is described by its
    ! global offset and local extent. We define a segment as being a single row
    ! of grid points in a block, i.e. the size of each segment is block_size_x.
    allocate(part_def(block_size_y*nblocks))
    part_idx = 1
    do n=1, nblocks
        this_block = get_block(blocks_ice(n), n)
        ilo = this_block%ilo
        jlo = this_block%jlo
        jhi = this_block%jhi

        do j = jlo, jhi
            ! Oasis uses zero-indexing to define the global offset, hence the final - 1
            part_def(part_idx).global_offset = ((this_block%j_glob(j) - 1) * nx_global) + this_block%i_glob(ilo) - 1
            part_def(part_idx).block_index = n
            part_idx = part_idx + 1
        enddo
    enddo

    ! The segments may not be in increasing order of global index. MCT
    ! doesn't like this and will try to reorder things which causes confusion
    ! because OASIS does not reverse the ordering when delivering the field.
    ! We keep this sorted.
    call sort_segments(part_def)

    ! Set up the partition definition in OASIS format.
    allocate(oasis_part_def(2 + 2*block_size_y*nblocks))
    ! Partition type
    oasis_part_def(1) = ORANGE
    ! The total number of segments
    oasis_part_def(2) = block_size_y*nblocks
    oasis_part_def(3::2) = part_def(:).global_offset
    oasis_part_def(4::2) = block_size_x

    call oasis_def_partition(part_id, oasis_part_def, err, nx_global * ny_global)

    ! Define coupling fields
    il_var_nodims(1) = 1 ! rank of coupling field
    il_var_nodims(2) = 1 ! number of bundles in coupling field (always 1)
    il_var_shape(1) = 1 ! min index for the coupling field local dimension
    il_var_shape(2) = block_size_x*block_size_y*nblocks

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

    do jf=1, jpfldout
        call oasis_def_var(il_var_id_out(jf),cl_writ(jf), part_id, &
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

    do jf=1, jpfldin
      call oasis_def_var(il_var_id_in(jf), cl_read(jf), part_id, &
                         il_var_nodims, PRISM_In, il_var_shape, PRISM_Real, ierror)
    enddo
    !
    ! 7- PSMILe end of declaration phase
    !
    call oasis_enddef(ierror, runtime=runtime_seconds, &
                      coupling_field_timesteps=coupling_field_timesteps)

    !
    ! Allocate the 'coupling' fields (to be used) for EACH PROCESS:!
    !

    ! fields in: (local domain)
    allocate (tair0(nx_block, ny_block, max_blocks));  tair0(:,:,:) = 0
    allocate (swflx0(nx_block, ny_block, max_blocks)); swflx0(:,:,:) = 0
    allocate (lwflx0(nx_block, ny_block, max_blocks)); lwflx0(:,:,:) = 0
    allocate (uwnd0(nx_block, ny_block, max_blocks));  uwnd0(:,:,:) = 0
    allocate (vwnd0(nx_block, ny_block, max_blocks));  vwnd0(:,:,:) = 0
    allocate (qair0(nx_block, ny_block, max_blocks));  qair0(:,:,:) = 0
    allocate (rain0(nx_block, ny_block, max_blocks));  rain0(:,:,:) = 0
    allocate (snow0(nx_block, ny_block, max_blocks));  snow0(:,:,:) = 0
    allocate (runof0(nx_block, ny_block, max_blocks)); runof0(:,:,:) = 0
    allocate (press0(nx_block, ny_block, max_blocks)); press0(:,:,:) = 0

    allocate (runof(nx_block, ny_block, max_blocks)); runof(:,:,:) = 0
    allocate (press(nx_block, ny_block, max_blocks)); press(:,:,:) = 0

    allocate (core_runoff(nx_block, ny_block, max_blocks));  core_runoff(:,:,:) = 0.

    allocate (ssto(nx_block, ny_block, max_blocks));  ssto(:,:,:) = 0
    allocate (ssso(nx_block, ny_block, max_blocks));  ssso(:,:,:) = 0
    allocate (ssuo(nx_block, ny_block, max_blocks));  ssuo(:,:,:) = 0
    allocate (ssvo(nx_block, ny_block, max_blocks));  ssvo(:,:,:) = 0
    allocate (sslx(nx_block, ny_block, max_blocks));  sslx(:,:,:) = 0
    allocate (ssly(nx_block, ny_block, max_blocks));  ssly(:,:,:) = 0
    allocate (pfmice(nx_block, ny_block, max_blocks));  pfmice(:,:,:) = 0

    allocate (iostrsu(nx_block, ny_block, max_blocks)); iostrsu(:,:,:) = 0
    allocate (iostrsv(nx_block, ny_block, max_blocks)); iostrsv(:,:,:) = 0
    allocate (iorain(nx_block, ny_block, max_blocks));  iorain(:,:,:) = 0
    allocate (iosnow(nx_block, ny_block, max_blocks));  iosnow(:,:,:) = 0
    allocate (iostflx(nx_block, ny_block, max_blocks)); iostflx(:,:,:) = 0
    allocate (iohtflx(nx_block, ny_block, max_blocks)); iohtflx(:,:,:) = 0
    allocate (ioswflx(nx_block, ny_block, max_blocks)); ioswflx(:,:,:) = 0
    allocate (ioqflux(nx_block, ny_block, max_blocks)); ioqflux(:,:,:) = 0
    allocate (iolwflx(nx_block, ny_block, max_blocks)); iolwflx(:,:,:) = 0
    allocate (ioshflx(nx_block, ny_block, max_blocks)); ioshflx(:,:,:) = 0
    allocate (iorunof(nx_block, ny_block, max_blocks)); iorunof(:,:,:) = 0
    allocate (iopress(nx_block, ny_block, max_blocks)); iopress(:,:,:) = 0
    allocate (ioaice (nx_block, ny_block, max_blocks)); ioaice(:,:,:) = 0

    allocate (iomelt (nx_block, ny_block, max_blocks)); iomelt(:,:,:) = 0
    allocate (ioform (nx_block, ny_block, max_blocks)); ioform(:,:,:) = 0

    allocate (tiostrsu(nx_block, ny_block, max_blocks)); tiostrsu(:,:,:) = 0
    allocate (tiostrsv(nx_block, ny_block, max_blocks)); tiostrsv(:,:,:) = 0
    allocate (tiorain(nx_block, ny_block, max_blocks));  tiorain(:,:,:) = 0
    allocate (tiosnow(nx_block, ny_block, max_blocks));  tiosnow(:,:,:) = 0
    allocate (tiostflx(nx_block, ny_block, max_blocks)); tiostflx(:,:,:) = 0
    allocate (tiohtflx(nx_block, ny_block, max_blocks)); tiohtflx(:,:,:) = 0
    allocate (tioswflx(nx_block, ny_block, max_blocks)); tioswflx(:,:,:) = 0
    allocate (tioqflux(nx_block, ny_block, max_blocks)); tioqflux(:,:,:) = 0
    allocate (tiolwflx(nx_block, ny_block, max_blocks)); tiolwflx(:,:,:) = 0
    allocate (tioshflx(nx_block, ny_block, max_blocks)); tioshflx(:,:,:) = 0
    allocate (tiorunof(nx_block, ny_block, max_blocks)); tiorunof(:,:,:) = 0
    allocate (tiopress(nx_block, ny_block, max_blocks)); tiopress(:,:,:) = 0
    allocate (tioaice(nx_block, ny_block, max_blocks));  tioaice(:,:,:) = 0

    allocate (tiomelt(nx_block, ny_block, max_blocks));  tiomelt(:,:,:) = 0
    allocate (tioform(nx_block, ny_block, max_blocks));  tioform(:,:,:) = 0

    allocate (vwork(nx_block, ny_block, max_blocks)); vwork(:,:,:) = 0
    allocate (gwork(nx_global, ny_global)); gwork(:,:) = 0
    allocate (vwork2d(l_ilo:l_ihi, l_jlo:l_jhi)); vwork2d(:,:) = 0.

    allocate (sicemass(nx_block, ny_block, max_blocks)); sicemass(:,:,:) = 0.
    allocate (u_star0(nx_block, ny_block, max_blocks)); u_star0(:,:,:) = 0.
    allocate (rough_mom0(nx_block, ny_block, max_blocks)); rough_mom0(:,:,:) = 0.
    allocate (rough_heat0(nx_block, ny_block, max_blocks)); rough_heat0(:,:,:) = 0.
    allocate (rough_moist0(nx_block, ny_block, max_blocks)); rough_moist0(:,:,:) = 0.

endsubroutine init_cpl

subroutine send_grid_to_atm()

  integer(kind=int_kind) :: tag, buf_int(2)
  real(kind=dbl_kind), dimension(:), allocatable :: buf_real

  integer :: fid
  real(kind=dbl_kind), dimension(:, :), allocatable :: tlat_global, tlon_global
  real(kind=dbl_kind), dimension(:, :), allocatable :: mask_global

  if (my_task == master_task) then
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

    ! Send my details to the atm.
    tag = 4982
    buf_int(1) = nx_global
    buf_int(2) = ny_global
    call MPI_send(buf_int, 2, MPI_INTEGER, coupler%atm_root, tag, &
                  MPI_COMM_WORLD, ierror)

    allocate(buf_real(nx_global*ny_global))
    buf_real(:) = reshape(tlat_global(:, :), (/ size(tlat_global) /))
    call MPI_send(buf_real, nx_global*ny_global, MPI_DOUBLE, &
                  coupler%atm_root, tag, MPI_COMM_WORLD, ierror)

    buf_real(:) = reshape(tlon_global(:, :), (/ size(tlon_global) /))
    call MPI_send(buf_real, nx_global*ny_global, MPI_DOUBLE, &
                  coupler%atm_root, tag, MPI_COMM_WORLD, ierror)

    buf_real(:) = reshape(mask_global(:, :), (/ size(mask_global) /))
    call MPI_send(buf_real, nx_global*ny_global, MPI_DOUBLE, &
                  coupler%atm_root, tag, MPI_COMM_WORLD, ierror)

    deallocate(buf_real)
    deallocate(tlat_global)
    deallocate(tlon_global)
    deallocate(mask_global)

  endif

end subroutine send_grid_to_atm

!> Convert 1d coupling arrays into 3d CICE arrays. The 1d coupling array is in
! the format set by the partition definition.
subroutine unpack_coupling_array(input, output)
    real, dimension(:), intent(in) :: input
    real, dimension(:, :, :), intent(inout) :: output

    integer :: isc, iec
    integer :: iseg
    integer :: iblk, offset
    integer, dimension(nblocks) :: blk_seg_num

    isc = 1+nghost
    iec = nx_block-nghost

    blk_seg_num(:) = 1 + nghost
    offset = 0

    do iseg=1, size(part_def)
        iblk = part_def(iseg).block_index

        output(isc:iec, blk_seg_num(iblk), iblk) = input((offset + 1):(offset + block_size_x))

        offset = offset + block_size_x
        blk_seg_num(iblk) = blk_seg_num(iblk) + 1
    enddo

endsubroutine unpack_coupling_array

!> Convert 3d CICE arrays into a 1d coupling array. The 1d coupling array needs
! to be in a specific format to match the partition definition that has been set
! up with Oasis.
subroutine pack_coupling_array(input, output)
    real, dimension(:, :, :), intent(in) :: input
    real, dimension(:), intent(inout) :: output

    integer :: isc, iec
    integer :: iseg
    integer :: iblk, offset
    integer, dimension(nblocks) :: blk_seg_num

    isc = 1+nghost
    iec = nx_block-nghost

    blk_seg_num(:) = 1 + nghost
    offset = 0

    ! Load the coupling array one segemnt at a time.  The code relies on the
    ! part_def being sorted so simply incrementing the blk_seg_num gets data
    ! from the next segmennt.
    do iseg=1, size(part_def)
        iblk = part_def(iseg).block_index

        output((offset + 1):(offset + block_size_x)) = input(isc:iec, blk_seg_num(iblk), iblk)

        offset = offset + block_size_x
        blk_seg_num(iblk) = blk_seg_num(iblk) + 1
    enddo

endsubroutine pack_coupling_array

subroutine from_atm(isteps)
  integer(kind=int_kind), intent(in) :: isteps

  integer(kind=int_kind) :: tag, request, info
  integer(kind=int_kind) :: buf(1)
  real(kind=dbl_kind), dimension(block_size_x*block_size_y*nblocks) :: work

#if defined(DEBUG)
    write(il_out,*) '(from_atm) receiving coupling fields at rtime= ', isteps
#endif

  call ice_timer_start(timer_from_atm)

  call ice_timer_start(timer_waiting_atm)
  call oasis_get(il_var_id_in(1), isteps, work, info)
  call unpack_coupling_array(work, swflx0)
  call ice_timer_stop(timer_waiting_atm)

  call oasis_get(il_var_id_in(2), isteps, work, info)
  call unpack_coupling_array(work, lwflx0)

  call oasis_get(il_var_id_in(3), isteps, work, info)
  call unpack_coupling_array(work, rain0)

  call oasis_get(il_var_id_in(4), isteps, work, info)
  call unpack_coupling_array(work, snow0)

  call oasis_get(il_var_id_in(5), isteps, work, info)
  call unpack_coupling_array(work, press0)

  call oasis_get(il_var_id_in(6), isteps, work, info)
  call unpack_coupling_array(work, runof0)

  call oasis_get(il_var_id_in(7), isteps, work, info)
  call unpack_coupling_array(work, tair0)

  call oasis_get(il_var_id_in(8), isteps, work, info)
  call unpack_coupling_array(work, qair0)

  call oasis_get(il_var_id_in(9), isteps, work, info)
  call unpack_coupling_array(work, uwnd0)

  call oasis_get(il_var_id_in(10), isteps, work, info)
  call unpack_coupling_array(work, vwnd0)

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
  if (my_task == master_task) then
    request = MPI_REQUEST_NULL
    tag = 5793
    call MPI_Isend(buf, 1, MPI_INTEGER, coupler%atm_root, tag, &
                   MPI_COMM_WORLD, request, ierror)
  endif

  call ice_timer_stop(timer_from_atm)

end subroutine from_atm

subroutine from_ocn(isteps)
    integer(kind=int_kind), intent(in) :: isteps

    integer :: info
    real(kind=dbl_kind), dimension(block_size_x*block_size_y*nblocks) :: work

#if defined(DEBUG)
    write(il_out,*) '(from_ocn) receiving coupling fields at rtime: ', isteps
#endif

    call ice_timer_start(timer_from_ocn)
    call ice_timer_start(timer_waiting_ocn)

    call oasis_get(il_var_id_in(11), isteps, work, info)
    call unpack_coupling_array(work, ssto)
    call ice_timer_stop(timer_waiting_ocn)

    call oasis_get(il_var_id_in(12), isteps, work, info)
    call unpack_coupling_array(work, ssso)

    call oasis_get(il_var_id_in(13), isteps, work, info)
    call unpack_coupling_array(work, ssuo)

    call oasis_get(il_var_id_in(14), isteps, work, info)
    call unpack_coupling_array(work, ssvo)

    call oasis_get(il_var_id_in(15), isteps, work, info)
    call unpack_coupling_array(work, sslx)

    call oasis_get(il_var_id_in(16), isteps, work, info)
    call unpack_coupling_array(work, ssly)

    call oasis_get(il_var_id_in(17), isteps, work, info)
    call unpack_coupling_array(work, pfmice)

    if (chk_o2i_fields) then
      call check_o2i_fields('fields_o2i_in_ice.nc',isteps)
    endif

    call ice_timer_stop(timer_from_ocn)  ! ice/ocn coupling

end subroutine from_ocn

!
! Note dummy 'scale', if /= 1 (then must be 1/coef_ic), is used here for the very 
! first-time-step-of-exp i2o fluxes scaling-up, because routine 'tavg_i2o_fluxes' 
! has scaled-down the current step i2o fluxes (calculated in get_i2o_fluxes) by 
! * coef_ic.  
!
subroutine into_ocn(isteps, scale)
 
    integer(kind=int_kind), intent(in) :: isteps
    real, intent(in) :: scale             !only 1 or 1/coef_ic allowed! 

    real(kind=dbl_kind), dimension(block_size_x*block_size_y*nblocks) :: work

    call pack_coupling_array(iostrsu*scale, work)
    call oasis_put(il_var_id_out(2), isteps, work, ierror)

    call pack_coupling_array(iostrsv*scale, work)
    call oasis_put(il_var_id_out(3), isteps, work, ierror)

    call pack_coupling_array(iorain*scale, work)
    call oasis_put(il_var_id_out(4), isteps, work, ierror)

    call pack_coupling_array(iosnow*scale, work)
    call oasis_put(il_var_id_out(5), isteps, work, ierror)

    call pack_coupling_array(iostflx*scale, work)
    call oasis_put(il_var_id_out(6), isteps, work, ierror)

    call pack_coupling_array(iohtflx*scale, work)
    call oasis_put(il_var_id_out(7), isteps, work, ierror)

    call pack_coupling_array(ioswflx*scale, work)
    call oasis_put(il_var_id_out(8), isteps, work, ierror)

    call pack_coupling_array(ioqflux*scale, work)
    call oasis_put(il_var_id_out(9), isteps, work, ierror)

    call pack_coupling_array(ioshflx*scale, work)
    call oasis_put(il_var_id_out(10), isteps, work, ierror)

    call pack_coupling_array(iolwflx*scale, work)
    call oasis_put(il_var_id_out(11), isteps, work, ierror)

    call pack_coupling_array(iorunof*scale, work)
    call oasis_put(il_var_id_out(12), isteps, work, ierror)

    call pack_coupling_array(iopress*scale, work)
    call oasis_put(il_var_id_out(13), isteps, work, ierror)

    call pack_coupling_array(ioaice*scale, work)
    call oasis_put(il_var_id_out(14), isteps, work, ierror)

    call pack_coupling_array(iomelt*scale, work)
    call oasis_put(il_var_id_out(15), isteps, work, ierror)

    call pack_coupling_array(ioform*scale, work)
    call oasis_put(il_var_id_out(16), isteps, work, ierror)

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

end subroutine update_halos_from_atm

!=======================================================================
  subroutine coupler_termination
!-------------------------------!

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
      print *, 'CICE: an error occured in prism_terminate = ', ierror
    endif
  else 
    if (my_task == 0) then
      print *,  '==================*** CICE END ***================='
    endif
  endif

#if defined(DEBUG)
  close(il_out)
#endif

    call MPI_Finalize (ierror)

  end subroutine coupler_termination

subroutine write_boundary_checksums(time)
   integer, intent(in) :: time

   integer isc, iec, jsc, jec

   if (my_task == master_task) then
     isc = 1+nghost
     iec = nx_block-nghost
     jsc = 1+nghost
     jec = ny_block-nghost

     print*, '[ice chksum] time:', time
     print*,   '[ice chksum] iostrsu:', sum(iostrsu(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] iostrsv:', sum(iostrsv(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] iorain:', sum(iorain(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] iosnow:', sum(iosnow(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] iostflx:', sum(iostflx(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] iohtflx:', sum(iohtflx(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ioswflx:', sum(ioswflx(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ioqflux:', sum(ioqflux(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ioshflx:', sum(ioshflx(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] iolwflx:', sum(iolwflx(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] iopress:', sum(iopress(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ioaice:', sum(ioaice(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] iomelt:', sum(iomelt(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ioform:', sum(ioform(isc:iec, jsc:jec, 1))

     print*,   '[ice chksum] ssto:', sum(ssto(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ssso:', sum(ssso(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ssuo:', sum(ssuo(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ssvo:', sum(ssvo(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] sslx:', sum(sslx(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] ssly:', sum(ssly(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] pfmice:', sum(pfmice(isc:iec, jsc:jec, 1))

     print*,   '[ice chksum] swflx0:', sum(swflx0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] lwflx0:', sum(lwflx0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] rain0:', sum(rain0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] snow0:', sum(snow0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] press0:',  sum(press0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] runof0:',  sum(runof0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] tair0:',  sum(tair0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] qair0:',  sum(qair0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] uwnd0:',  sum(uwnd0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] vwnd0:',  sum(vwnd0(isc:iec, jsc:jec, 1))

     print*,   '[ice chksum] u_star0:', sum(u_star0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] rough_mom0:', sum(rough_mom0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] rough_heat0:', sum(rough_heat0(isc:iec, jsc:jec, 1))
     print*,   '[ice chksum] rough_moist0:', sum(rough_moist0(isc:iec, jsc:jec, 1))
  endif

end subroutine

end module cpl_interface


