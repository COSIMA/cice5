!============================================================================
  module cpl_interface
!============================================================================
! coupling interface between CICE and the oasis3_25 coupler (via MPI2) using 
! the PRISM System Model Interface (PSMILe).
!----------------------------------------------------------------------------

  !prism stuff
  use mod_prism

  !cice stuff
  use ice_kinds_mod
  use ice_communicate !, only : my_task, master_task
  use ice_broadcast
  use ice_blocks       !, only : nx_block, ny_block, nghost
  use ice_domain_size  !, only : max_blocks, nx_global, ny_global, ncat
  use ice_distribution, only : distrb, nprocsX, nprocsY
  use ice_gather_scatter
  use ice_constants
  use ice_domain       !, only : distrb_info
  use ice_grid,        only : u2tgrid_vector
  use ice_grid,        only : ANGLE, ANGLET 
  use ice_exit,        only : abort_ice

  !cpl stuff
  use cpl_parameters
  use cpl_netcdf_setup
  use cpl_arrays_setup
  use cpl_forcing_handler

  implicit none

  public :: prism_init, init_cpl, coupler_termination, get_time0_sstsss, &
            from_atm, into_ocn, from_ocn, into_atm, save_restart_i2a

  private

  logical :: mpiflag
  integer(kind=int_kind)  :: ierror, ibou
  character(len=9) :: chiceout
  character(len=3) :: chout
  logical :: ll_comparal    ! paralell or mono-cpl coupling
  integer(kind=int_kind) :: il_comp_id     ! Component ID
  integer(kind=int_kind) :: il_nbtotproc   ! Total number of processes
  integer(kind=int_kind) :: il_nbcplproc   ! No of processes involved in coupling
  integer(kind=int_kind) :: il_part_id     ! Local partition ID
  integer(kind=int_kind) :: il_length      ! Size of partial field for each process
  integer(kind=int_kind), dimension(2) :: il_var_nodims
  integer(kind=int_kind), dimension(4) :: il_var_shape

  integer(kind=int_kind) :: l_ilo, l_ihi, l_jlo, l_jhi !local partition
  integer(kind=int_kind) :: gh_ilo, gh_ihi, gh_jlo, gh_jhi !local ghost outline 
  integer :: sendsubarray, recvsubarray , resizedrecvsubarray
  integer, dimension(:), allocatable :: counts, disps 

  integer(kind=int_kind) :: il_flag        ! Flag for grid writing
  integer(kind=int_kind) :: il_status, il_fileid, il_varid
  integer(kind=int_kind)      :: io_size, ii, il_bufsize, il_real, il_bufsizebyt
  integer(kind=int_kind)      :: integer_byte_size, integer_io_size
  real(kind=dbl_kind), dimension(:,:), allocatable :: rla_array
  real(kind=dbl_kind), dimension(:), allocatable :: rla_bufsend
  real(kind=dbl_kind), dimension(:,:), allocatable :: vwork2d
    !local domain work array, 4 coupling data passing 
  contains

!======================================================================
  subroutine prism_init
!-----------------------!

  include 'mpif.h'

  !-----------------------------------
  ! 'define' the model global domain: 
  !-----------------------------------
  il_im = nx_global
  il_jm = ny_global
  il_imjm = il_im * il_jm

  !allocate rla_array to be used below
  allocate (rla_array(il_im,il_jm) )

  !print *, 'CICE: (prism_init) dbl_kind, ip_realwp_p= ',dbl_kind, ip_realwp_p 

  !-------------------
  ! Initialize PSMILe.
  !-------------------

  ! Initialise MPI
  mpiflag = .FALSE.
  call MPI_Initialized (mpiflag, ierror)
  print *,'CICE: (prism_init) BF MPI_INIT, mpiflag = ',mpiflag

  if ( .not. mpiflag ) then
    call MPI_INIT(ierror)
  endif

  call MPI_Initialized (mpiflag, ierror)
  print *, 'CICE: (prism_init) AF MPI_INIT, mpiflag = ',mpiflag

  print *
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
  il_bufsize = il_imjm + MPI_BSEND_OVERHEAD/il_real + 1
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
 ! print *, '* CICE: Entering init_cpl.....'

  print *, '* CICE (prism_init) calling MPI_Comm_Size ...'
  call MPI_Comm_Size(il_commlocal, il_nbtotproc, ierror)
  print *, '* CICE (prism_init) calling MPI_Comm_Rank ...'
  call MPI_Comm_Rank(il_commlocal, my_task, ierror)

  print *, '* CICE (prism_init) il_commlocal, il_nbtotproc, my_task = '
  print *, '* CICE (prism_init) ', il_commlocal, il_nbtotproc, my_task
  !
  il_nbcplproc = il_nbtotproc   !multi-process coupling (real parallel cpl)!
  !il_nbcplproc = 1               !mono process coupling

  if (il_nbtotproc /= 1 .and. il_nbcplproc == il_nbtotproc ) then
    ll_comparal = .TRUE.          ! multi-cpl coupling!
  else
    ll_comparal = .FALSE.          !mono-cpl coupling!
  endif

  print *, '* CICE: prism_init called OK!'  

  end subroutine prism_init

!=======================================================================
  subroutine init_cpl

  use mpi
  use ice_communicate
!--------------------!
  integer(kind=int_kind) :: jf, jfs
  integer(kind=int_kind), dimension(2) :: il_var_nodims ! see below
  integer(kind=int_kind), dimension(4) :: il_var_shape  ! see below

  integer(kind=int_kind) :: ilo,ihi,jlo,jhi,iblk,i,j, n
  type (block) ::  this_block           ! block information for current block

  integer, dimension(2) :: starts,sizes,subsizes
  integer(kind=mpi_address_kind) :: start, extent
!  integer, dimension(:), allocatable :: counts, disps 
  real(kind=dbl_kind) :: realvalue
  integer (int_kind) :: nprocs
  integer (int_kind),dimension(:), allocatable :: vilo, vjlo 

  nprocs = get_num_procs()
  allocate(vilo(nprocs))
  allocate(vjlo(nprocs))
!initialise partition to inexisting region
 l_ilo=nx_global
 l_ihi=0
 l_jlo=ny_global
 l_jhi=0
 gh_ilo=nx_global
 gh_ihi=0
 gh_jlo=ny_global
 gh_jhi=0
  ! Open the process log file
!20100406  if (my_task == 0 .or. ll_comparal) then
    il_out = 85 + my_task
    write(chout,'(I3.3)')il_out
    chiceout='iceout'//chout
    open(il_out,file=chiceout,form='formatted')
  
    write(il_out,*) 'Number of processes:', il_nbtotproc
    write(il_out,*) 'Local process number:', my_task
    write(il_out,*) 'Local communicator is : ',il_commlocal
    write(il_out,*) 'Grid layout: nx_global,ny_global= ',nx_global,ny_global
    write(il_out,*) 'Grid decomposition: nx_block,ny_block,max_blocks= ',&
                     nx_block,ny_block,max_blocks
!20100406  endif

!  write(il_out,*) 'Number of blocks :', nblocks
!  do iblk = 1, nblocks
!
!    this_block = get_block(blocks_ice(iblk),iblk)
!    ilo = this_block%ilo
!    ihi = this_block%ihi
!    jlo = this_block%jlo
!    jhi = this_block%jhi
!!         do j=this_block%jlo,this_block%jhi
!!         do i=this_block%ilo,this_block%ihi
!!           ARRAY_G(this_block%i_glob(i), &
!!                   this_block%j_glob(j)) = &
!!                  ARRAY(i,j,src_dist%blockLocalID(n))
!
!    write(il_out,*) '   block:', iblk, "ilo, jlo, ihi, jhi=", ilo, jlo, ihi, jhi
!
!  end do
!find out partition of this processor, which is done by init_domain_blocks
  write(il_out,*) 'nblocks_x, nblocks_y, Number of tot blocks :', nblocks_x, nblocks_y, nblocks_tot
!!!!!!!!!!!!
!     do iblk=1,nblocks_tot
!        if (distrb_info%blockLocation(iblk) == my_task+1) then
!          this_block = get_block(iblk,iblk)
!          ilo = this_block%ilo
!          ihi = this_block%ihi
!          jlo = this_block%jlo
!          jhi = this_block%jhi
!          !print local to global mapping
!          write(il_out,*) 'block, local ilo ihi jlo jhi:', distrb_info%blockLocalID(iblk), ilo,ihi,jlo,jhi
!          write(il_out,*) 'block global:', this_block%i_glob(ilo),this_block%i_glob(ihi),  &
!              this_block%j_glob(jlo),this_block%j_glob(jhi)
!        endif
!     end do

!!debug
  if (my_task == 0) then
  write(il_out,*) "all block info:"
  do iblk=1,nblocks_tot
      this_block = get_block(iblk,iblk)
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      write(il_out,*) '   this block: cpu, iblock, jblock=', distrb_info%blockLocation(iblk)-1, this_block%iblock, this_block%jblock
      write(il_out,*) '   block:', iblk, "gilo, gjlo, gihi, gjhi=", this_block%i_glob(ilo), this_block%j_glob(jlo), this_block%i_glob(ihi), this_block%j_glob(jhi)
  end do
  end if

  do iblk=1,nblocks_tot

    if (distrb_info%blockLocation(iblk) == my_task+1) then

      this_block = get_block(iblk,iblk)
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi
      write(il_out,*) '   this block: iblock, jblock=', this_block%iblock, this_block%jblock
!      write(il_out,*) '   block:', iblk, "ilo, jlo, ihi, jhi=", ilo, jlo, ihi, jhi
      write(il_out,*) '   block:', iblk, "gilo, gjlo, gihi, gjhi=", this_block%i_glob(ilo), this_block%j_glob(jlo), this_block%i_glob(ihi), this_block%j_glob(jhi)
      if (this_block%i_glob(ilo) < l_ilo) then
        l_ilo = this_block%i_glob(ilo)
        gh_ilo = this_block%i_glob(ilo-nghost)
      endif
      if (this_block%j_glob(jlo) < l_jlo) then
        l_jlo = this_block%j_glob(jlo)
        gh_jlo = this_block%j_glob(jlo-nghost)
      endif
      if (this_block%i_glob(ihi) > l_ihi) then
        l_ihi = this_block%i_glob(ihi)
        gh_ihi = this_block%i_glob(ihi+nghost)
      endif
      if (this_block%j_glob(jhi) > l_jhi) then
        l_jhi = this_block%j_glob(jhi)
        gh_jhi = this_block%j_glob(jhi+nghost)
      endif
!      l_ilo = min(l_ilo, this_block%i_glob(ilo))
!      l_ihi = max(l_ihi, this_block%i_glob(ihi))
!      l_jlo = min(l_jlo, this_block%j_glob(jlo))
!      l_jhi = max(l_jhi, this_block%j_glob(jhi))
!    else if (distrb_info%blockLocation(n) == 0) then 
!       write(il_out,*) '  land block:', n
       
    endif
  end do
  write(il_out,*) '  local partion, ilo, ihi, jlo, jhi=', l_ilo, l_ihi, l_jlo, l_jhi
  write(il_out,*) '  partition x,y sizes:', l_ihi-l_ilo+1, l_jhi-l_jlo+1
!print ghost info
  write(il_out,*) '  ghost global:',gh_ilo, gh_ihi, gh_jlo, gh_jhi 

!calculate partition using nprocsX and nprocsX
  l_ilo=mod(my_task,nprocsX)*nx_global/nprocsX+1
  l_ihi=l_ilo + nx_global/nprocsX -1
  l_jlo=int(my_task/nprocsX) * ny_global/nprocsY+1
  l_jhi=l_jlo+ny_global/nprocsY - 1
 
  write(il_out,*) '  2local partion, ilo, ihi, jlo, jhi=', l_ilo, l_ihi, l_jlo, l_jhi
  write(il_out,*) '  2partition x,y sizes:', l_ihi-l_ilo+1, l_jhi-l_jlo+1
 
  call mpi_gather(l_ilo, 1, mpi_integer, vilo, 1, mpi_integer, 0, MPI_COMM_ICE, ierror)
  call broadcast_array(vilo, 0)  
  call mpi_gather(l_jlo, 1, mpi_integer, vjlo, 1, mpi_integer, 0, MPI_COMM_ICE, ierror)
  call broadcast_array(vjlo, 0)  
 
!create subarray of this rank
    sizes(1)=l_ihi-l_ilo+1; sizes(2)=l_jhi-l_jlo+1
    subsizes(1)=l_ihi-l_ilo+1; subsizes(2)=l_jhi-l_jlo+1
    starts(1)=0; starts(2)=0
    call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, &
                                MPI_REAL8, sendsubarray, ierror)
    call mpi_type_commit(sendsubarray,ierror)
    if (my_task == 0) then ! create recv buffer in main cpu
     sizes(1)=nx_global; sizes(2)=ny_global
     subsizes(1)=l_ihi-l_ilo+1; subsizes(2)=l_jhi-l_jlo+1
     starts(1)=0; starts(2)=0
     call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, &
                                   MPI_REAL8, recvsubarray, ierror)
     call mpi_type_commit(recvsubarray, ierror)
     extent = sizeof(realvalue)
     start = 0
     call mpi_type_create_resized(recvsubarray, start, extent, resizedrecvsubarray, ierror)
     call mpi_type_commit(resizedrecvsubarray,ierror)
    end if
    allocate(counts(nprocs),disps(nprocs))
    forall (n=1:nprocs) counts(n) = 1
    do n=1,  nprocs
      disps(n) = ((vjlo(n)-1)*nx_global + (vilo(n)-1)) 
      !disps(n) = ((vilo(n)-1)*ny_global + (vjlo(n)-1)) 
    end do
    
    write(il_out,*) ' vilo ', vilo
    write(il_out,*) ' vjlo ', vjlo
    write(il_out,*) ' counts ', counts
    write(il_out,*) ' disps ', disps
 
!  if ( ll_comparal ) then 
!    il_im = l_ihi-l_ilo+1 !nx_global
!    il_jm = l_jhi-l_jlo+1 !ny_global
!    il_imjm = il_im * il_jm
!  endif
    if (ll_comparal) then
      xdim=l_ihi-l_ilo+1
      ydim=l_jhi-l_jlo+1
    else
      xdim=il_im
      ydim=il_jm
    endif
  

!-----------------------------------------------------------------------
  if (my_task == 0 .or. ll_comparal) then
    !
    ! The following steps need to be done:
    ! -> by the process if cice is monoprocess;
    ! -> only by the master process, if cice is parallel and only 
    !    master process is involved in the coupling;
    ! -> by all processes, if cice is parallel and all processes 
    ! are involved in the coupling.
    
    call decomp_def (il_part_id, il_length, il_imjm, &
         my_task, il_nbcplproc, ll_comparal, il_out)

    write(il_out,*)'(init_cpl) called decomp_def, my_task, ierror = ',my_task, ierror

    !
    ! PSMILe coupling fields declaration
    !

    il_var_nodims(1) = 2 ! rank of coupling field
    il_var_nodims(2) = 1 ! number of bundles in coupling field (always 1)
    !il_var_shape(1)= 1     ! min index for the coupling field local dim
    !il_var_shape(2)= xdim !il_im ! max index for the coupling field local dim  
    !il_var_shape(3)= 1   
    !il_var_shape(4)= ydim !il_jm
    if (ll_comparal) then 
      il_var_shape(1)= 1 !l_ilo ! min index for the coupling field local dimension
      il_var_shape(2)= l_ihi-l_ilo+1 ! max index for the coupling field local dim
      il_var_shape(3)= 1 !l_jlo ! min index for the coupling field local dim
      il_var_shape(4)= l_jhi-l_jlo+1 ! max index for the coupling field local dim
    else
      il_var_shape(1)= 1 ! min index for the coupling field local dimension
      il_var_shape(2)= il_im ! max index for the coupling field local dim
      il_var_shape(3)= 1 ! min index for the coupling field local dim
      il_var_shape(4)= il_jm ! max index for the coupling field local dim
    endif

    ! ?Does this help?
    !il_var_shape(1)= 2       ! min index for the coupling field local dim
    !il_var_shape(2)= il_im+1 ! max index for the coupling field local dim
    !il_var_shape(3)= 2  
    !il_var_shape(4)= il_jm+1

  endif   !my_task==0 .or. ll_comparal

  !*** 								***!
  !***B: we now define cl_writ/cl_read on all ranks! (20090403) ***!
  !***  							***!   
  
  !
  ! Define name (as in namcouple) and declare each field sent by ice 
  !

    !
    ! ice ==> atm
    !
    nsend_i2a = 0

    nsend_i2a = nsend_i2a + 1
    cl_writ(nsend_i2a)='isst_ia'
    do jf = 1, ncat
      nsend_i2a = nsend_i2a + 1
      write(cl_writ(nsend_i2a), '(a6,i2.2)')'icecon',jf
    enddo
    do jf = 1, ncat
      nsend_i2a = nsend_i2a + 1
      write(cl_writ(nsend_i2a), '(a6,i2.2)')'snwthk',jf
    enddo
    do jf = 1, ncat
      nsend_i2a = nsend_i2a + 1
      write(cl_writ(nsend_i2a), '(a6,i2.2)')'icethk',jf
    enddo
    nsend_i2a = nsend_i2a + 1
    cl_writ(nsend_i2a)='uvel_ia'
    nsend_i2a = nsend_i2a + 1
    cl_writ(nsend_i2a)='vvel_ia'
    nsend_i2a = nsend_i2a + 1
    cl_writ(nsend_i2a)='co2_i2'
    nsend_i2a = nsend_i2a + 1
    cl_writ(nsend_i2a)='co2fx_i2'
 
    if (my_task == 0) then
      write(il_out,*) 'init_cpl: Number of fields sent to atm: ',nsend_i2a
    endif
    !
    ! ice ==> ocn
    !
    nsend_i2o = nsend_i2a

    nsend_i2o = nsend_i2o + 1 
    cl_writ(nsend_i2o)='strsu_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='strsv_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='rain_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='snow_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='stflx_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='htflx_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='swflx_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='qflux_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='shflx_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='lwflx_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='runof_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='press_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='aice_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='melt_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='form_io'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='co2_i1'
    nsend_i2o = nsend_i2o + 1
    cl_writ(nsend_i2o)='wnd_i1'

    if (my_task == 0 .or. ll_comparal) then

      write(il_out,*) 'init_cpl: Number of fields sent to ocn: ',nsend_i2o - nsend_i2a

      if (nsend_i2o /= jpfldout) then
        write(il_out,*)
        write(il_out,*)'!!! Fatal Error: (init_cpl) nsend = ',nsend_i2o
        write(il_out,*)'!!!           It should be  nsend = ',jpfldout
        call abort_ice('CICE: Number of outgoing coupling fields incorrect!') 
      endif

      write(il_out,*) 'init_cpl: Total number of fields sent from ice: ',jpfldout

      !jpfldout == nsend_i2o!
      !---------------------!

      do jf=1, jpfldout
        call prism_def_var_proto (il_var_id_out(jf),cl_writ(jf), il_part_id, &
           il_var_nodims, PRISM_Out, il_var_shape, PRISM_Real, ierror)
      enddo 

    endif

    !
    ! Define name (as in namcouple) and declare each field received by ice
    !

    !
    ! atm ==> ice
    !
    nrecv_a2i = 0

    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'thflx_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'pswflx_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'runoff_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'wme_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'rain_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'snow_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'evap_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'lhflx_i'
    do jf = 1, ncat
      nrecv_a2i = nrecv_a2i + 1
      write(cl_read(nrecv_a2i), '(a4,i2.2,a2)')'tmlt',jf,'_i'
    enddo
    do jf = 1, ncat
      nrecv_a2i = nrecv_a2i + 1
      write(cl_read(nrecv_a2i), '(a4,i2.2,a2)')'bmlt',jf,'_i'
    enddo
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'taux_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'tauy_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'swflx_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'lwflx_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'shflx_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'press_i'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'co2_ai'
    nrecv_a2i = nrecv_a2i + 1
    cl_read(nrecv_a2i) = 'wnd_ai'

    if (my_task==0 .or. ll_comparal) then
      write(il_out,*) 'init_cpl: Number of fields rcvd from atm: ',nrecv_a2i
    endif

    !
    ! ocn ==> ice
    !
    nrecv_o2i = nrecv_a2i

    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'sst_i'
    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'sss_i'
    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'ssu_i'
    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'ssv_i'
    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'sslx_i'
    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'ssly_i'
    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'pfmice_i'
    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'co2_oi'
    nrecv_o2i = nrecv_o2i + 1
    cl_read(nrecv_o2i) = 'co2fx_oi'

    if (my_task==0 .or. ll_comparal) then

      write(il_out,*) 'init_cpl: Number of fields rcvd from ocn: ',nrecv_o2i-nrecv_a2i

      if (nrecv_o2i /= jpfldin) then
        write(il_out,*)
        write(il_out,*)'!!! Fatal Error: (init_cpl) nrecv = ',nrecv_o2i
        write(il_out,*)'!!!           It should be  nrecv = ',jpfldin
        call abort_ice('CICE: Number of incoming coupling fields incorrect!')
      endif
      !jpfldin == nrecv_o2i!
      !--------------------!
    
      write(il_out,*) 'init_cpl: Total number of fields rcvd by ice: ',jpfldin

      do jf=1, jpfldin
        call prism_def_var_proto (il_var_id_in(jf), cl_read(jf), il_part_id, &
           il_var_nodims, PRISM_In, il_var_shape, PRISM_Real, ierror)
      enddo 

      !
      ! PSMILe end of declaration phase 
      !
      call prism_enddef_proto (ierror)

    endif     !my_task==0

  !
  ! Allocate the 'coupling' fields (to be used) for EACH PROCESS:! 
  !

  ! fields in: (local domain)
  !
  ! from atm:
  allocate (um_thflx(nx_block,ny_block,max_blocks));  um_thflx(:,:,:) = 0
  allocate (um_pswflx(nx_block,ny_block,max_blocks)); um_pswflx(:,:,:) = 0
  allocate (um_runoff(nx_block,ny_block,max_blocks)); um_runoff(:,:,:) = 0
  allocate (um_wme(nx_block,ny_block,max_blocks)); um_wme(:,:,:) = 0
  allocate (um_snow(nx_block,ny_block,max_blocks)); um_snow(:,:,:) = 0
  allocate (um_rain(nx_block,ny_block,max_blocks)); um_rain(:,:,:) = 0
  allocate (um_evap(nx_block,ny_block,max_blocks)); um_evap(:,:,:) = 0
  allocate (um_lhflx(nx_block,ny_block,max_blocks)); um_lhflx(:,:,:) = 0
  allocate (um_taux(nx_block,ny_block,max_blocks)); um_taux(:,:,:) = 0
  allocate (um_tauy(nx_block,ny_block,max_blocks)); um_tauy(:,:,:) = 0
  allocate (um_swflx(nx_block,ny_block,max_blocks)); um_swflx(:,:,:) = 0
  allocate (um_lwflx(nx_block,ny_block,max_blocks)); um_lwflx(:,:,:) = 0
  allocate (um_shflx(nx_block,ny_block,max_blocks)); um_shflx(:,:,:) = 0
  allocate (um_press(nx_block,ny_block,max_blocks)); um_press(:,:,:) = 0
  allocate (um_tmlt(nx_block,ny_block,ncat,max_blocks)); um_tmlt(:,:,:,:) = 0
  allocate (um_bmlt(nx_block,ny_block,ncat,max_blocks)); um_bmlt(:,:,:,:) = 0
  allocate (um_co2(nx_block,ny_block,max_blocks)); um_co2(:,:,:) = 0
  allocate (um_wnd(nx_block,ny_block,max_blocks)); um_wnd(:,:,:) = 0

  !
  allocate ( core_runoff(nx_block,ny_block,max_blocks));  core_runoff(:,:,:) = 0.
  !

  ! from ocn:
  allocate (ocn_sst(nx_block,ny_block,max_blocks)); ocn_sst(:,:,:) = 0
  allocate (ocn_sss(nx_block,ny_block,max_blocks)); ocn_sss(:,:,:) = 0
  allocate (ocn_ssu(nx_block,ny_block,max_blocks)); ocn_ssu(:,:,:) = 0
  allocate (ocn_ssv(nx_block,ny_block,max_blocks)); ocn_ssv(:,:,:) = 0
  allocate (ocn_sslx(nx_block,ny_block,max_blocks)); ocn_sslx(:,:,:) = 0
  allocate (ocn_ssly(nx_block,ny_block,max_blocks)); ocn_ssly(:,:,:) = 0
  allocate (ocn_pfmice(nx_block,ny_block,max_blocks)); ocn_pfmice(:,:,:) = 0
  allocate (ocn_co2(nx_block,ny_block,max_blocks)); ocn_co2(:,:,:) = 0
  allocate (ocn_co2fx(nx_block,ny_block,max_blocks)); ocn_co2fx(:,:,:) = 0

  ! fields out: (local domain)
  !
  ! to atm:
  allocate (ia_sst(nx_block,ny_block,max_blocks));  ia_sst(:,:,:) = 0
  allocate (ia_uvel(nx_block,ny_block,max_blocks)); ia_uvel(:,:,:) = 0
  allocate (ia_vvel(nx_block,ny_block,max_blocks)); ia_vvel(:,:,:) = 0
  allocate (ia_aicen(nx_block,ny_block,ncat,max_blocks)); ia_aicen(:,:,:,:) = 0
  allocate (ia_snown(nx_block,ny_block,ncat,max_blocks)); ia_snown(:,:,:,:) = 0
  allocate (ia_thikn(nx_block,ny_block,ncat,max_blocks)); ia_thikn(:,:,:,:) = 0
  allocate (ia_co2(nx_block,ny_block,max_blocks)); ia_co2(:,:,:) = 0
  allocate (ia_co2fx(nx_block,ny_block,max_blocks)); ia_co2fx(:,:,:) = 0
  !
  ! to ocn:
  allocate (io_strsu(nx_block,ny_block,max_blocks)); io_strsu(:,:,:) = 0
  allocate (io_strsv(nx_block,ny_block,max_blocks)); io_strsv(:,:,:) = 0
  allocate (io_rain (nx_block,ny_block,max_blocks)); io_rain (:,:,:) = 0
  allocate (io_snow (nx_block,ny_block,max_blocks)); io_snow (:,:,:) = 0
  allocate (io_stflx(nx_block,ny_block,max_blocks)); io_stflx(:,:,:) = 0
  allocate (io_htflx(nx_block,ny_block,max_blocks)); io_htflx(:,:,:) = 0
  allocate (io_swflx(nx_block,ny_block,max_blocks)); io_swflx(:,:,:) = 0
  allocate (io_qflux(nx_block,ny_block,max_blocks)); io_qflux(:,:,:) = 0
  allocate (io_lwflx(nx_block,ny_block,max_blocks)); io_lwflx(:,:,:) = 0
  allocate (io_shflx(nx_block,ny_block,max_blocks)); io_shflx(:,:,:) = 0
  allocate (io_runof(nx_block,ny_block,max_blocks)); io_runof(:,:,:) = 0
  allocate (io_press(nx_block,ny_block,max_blocks)); io_press(:,:,:) = 0
  allocate (io_aice(nx_block,ny_block,max_blocks));  io_aice(:,:,:) = 0
  allocate (io_melt(nx_block,ny_block,max_blocks));  io_melt(:,:,:) = 0
  allocate (io_form(nx_block,ny_block,max_blocks));  io_form(:,:,:) = 0
  allocate (io_co2(nx_block,ny_block,max_blocks));  io_co2(:,:,:) = 0
  allocate (io_wnd(nx_block,ny_block,max_blocks));  io_wnd(:,:,:) = 0

  ! temporary arrays:
  ! IO cpl int time-average
  allocate (maice(nx_block,ny_block,max_blocks)); maice(:,:,:) = 0
  allocate (mstrocnxT(nx_block,ny_block,max_blocks)); mstrocnxT(:,:,:) = 0
  allocate (mstrocnyT(nx_block,ny_block,max_blocks)); mstrocnyT(:,:,:) = 0
  allocate (mfresh(nx_block,ny_block,max_blocks)); mfresh(:,:,:) = 0
  allocate (mfsalt(nx_block,ny_block,max_blocks)); mfsalt(:,:,:) = 0
  allocate (mfhocn(nx_block,ny_block,max_blocks)); mfhocn(:,:,:) = 0
  allocate (mfswthru(nx_block,ny_block,max_blocks)); mfswthru(:,:,:) = 0
  allocate (msicemass(nx_block,ny_block,max_blocks)); msicemass(:,:,:) = 0
  ! IA cpl int time-average (3D)
  allocate (maiu(nx_block,ny_block,max_blocks));  maiu(:,:,:) = 0
  allocate (muvel(nx_block,ny_block,max_blocks)); muvel(:,:,:) = 0
  allocate (mvvel(nx_block,ny_block,max_blocks)); mvvel(:,:,:) = 0
  allocate (msst(nx_block,ny_block,max_blocks));  msst(:,:,:) = 0    
  allocate (mssu(nx_block,ny_block,max_blocks));  mssu(:,:,:) = 0
  allocate (mssv(nx_block,ny_block,max_blocks));  mssv(:,:,:) = 0
  allocate (mco2(nx_block,ny_block,max_blocks));  mco2(:,:,:) = 0
  allocate (mco2fx(nx_block,ny_block,max_blocks));  mco2fx(:,:,:) = 0
! IA cpl int time-average (4D)
  allocate (maicen(nx_block,ny_block,ncat,max_blocks)); maicen(:,:,:,:) = 0
  allocate (msnown(nx_block,ny_block,ncat,max_blocks)); msnown(:,:,:,:) = 0
  allocate (mthikn(nx_block,ny_block,ncat,max_blocks)); mthikn(:,:,:,:) = 0

  allocate (vwork(nx_block,ny_block,max_blocks)); vwork(:,:,:) = 0
  allocate (gwork(nx_global,ny_global)); gwork(:,:) = 0
  allocate (sicemass(nx_block,ny_block,max_blocks)); sicemass(:,:,:) = 0.
  allocate (vwork2d(l_ilo:l_ihi, l_jlo:l_jhi)); vwork2d(:,:) = 0. !l_ihi-l_ilo+1, l_jhi-l_jlo+1
  end subroutine init_cpl

!=======================================================================
  subroutine from_atm(isteps)
!----------------------------!

  implicit none

  integer(kind=int_kind), intent(in) :: isteps

  real(kind=dbl_kind) :: tmpu, tmpv
  integer(kind=int_kind) :: ilo,ihi,jlo,jhi,iblk,i,j
  type (block) :: this_block           ! block information for current block

  integer(kind=int_kind) :: jf
  integer(kind=int_kind) :: ncid,currstep,ll,ilout
  
  data currstep/0/
  save currstep

  currstep=currstep+1

  if (my_task == 0) then  
    write(il_out,*)
    write(il_out,*) '(from_atm) receiving coupling fields at rtime= ', isteps
    if (chk_a2i_fields) then
      if ( .not. file_exist('fields_a2i_in_ice.nc') ) then
        call create_ncfile('fields_a2i_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
      endif
      write(il_out,*) 'opening file fields_a2i_in_ice.nc at nstep = ', isteps
      call ncheck( nf_open('fields_a2i_in_ice.nc',nf_write,ncid) )
      call write_nc_1Dtime(real(isteps),currstep,'time',ncid)
    endif
    write(il_out,*)
    write(il_out,*) '(from_atm) Total number of fields to be rcvd: ', nrecv_a2i
  endif
  
  write(il_out,*) "prism_get from_atm at sec: ", isteps
  do jf = 1, nrecv_a2i

    if (my_task==0 .or. ll_comparal ) then

      !jf-th field in
      write(il_out,*)
      write(il_out,*) '*** receiving coupling field No. ', jf, cl_read(jf)
      !call flush(il_out)

      if (ll_comparal) then 
        call prism_get_proto (il_var_id_in(jf), isteps, vwork2d(l_ilo:l_ihi, l_jlo:l_jhi), ierror) !vwork(2:,2:,my_task+1), 
        call mpi_gatherv(vwork2d(l_ilo:l_ihi, l_jlo:l_jhi),1,sendsubarray,gwork,counts,disps,resizedrecvsubarray, &
                     0,MPI_COMM_ICE,ierror)
        call broadcast_array(gwork, 0)
!         gwork(l_ilo:l_ihi, l_jlo:l_jhi) = vwork2d(l_ilo:l_ihi, l_jlo:l_jhi)
      else
        call prism_get_proto (il_var_id_in(jf), isteps, gwork, ierror)
      endif 

      if ( ierror /= PRISM_Ok .and. ierror < PRISM_Recvd) then
        write(il_out,*) 'Err in _get_ sst at time with error: ', isteps, ierror
        call prism_abort_proto(il_comp_id, 'cice from_atm','stop 1') 
      else 
        write(il_out,*)
        write(il_out,*)'(from_atm) rcvd at time with err: ',cl_read(jf),isteps,ierror
     
        if (ll_comparal .and. chk_a2i_fields) then
           call mpi_gatherv(vwork2d(l_ilo:l_ihi, l_jlo:l_jhi),1,sendsubarray,gwork, &
                     counts,disps, resizedrecvsubarray, 0,MPI_COMM_ICE,ierror)
           call MPI_Barrier(MPI_COMM_ICE, ierror)
        endif
        if (my_task==0 .and. chk_a2i_fields) then
            call write_nc2D(ncid, trim(cl_read(jf)), gwork, 1, il_im,il_jm, &
                          currstep,ilout=il_out)
        endif
      endif

    endif

    if (.not. ll_comparal ) then
      call scatter_global(vwork,gwork,master_task,distrb_info, &
                        field_loc_center, field_type_scalar)
    else
      call unpack_global_dbl(vwork,gwork,master_task,distrb_info, &
                        field_loc_center, field_type_scalar)
    endif ! not ll_comparal

    !***Note following "select case" works only if cl_read(:) is defined at ALL ranks***!
    !-----------------------------------------------------------------------------------!
    select case (trim(cl_read(jf)))
    case ('thflx_i');  um_thflx(:,:,:)  = vwork(:,:,:)
    case ('pswflx_i'); um_pswflx(:,:,:) = vwork(:,:,:)
    case ('runoff_i'); um_runoff(:,:,:) = vwork(:,:,:)
    case ('wme_i');   um_wme(:,:,:) = vwork(:,:,:)
!    case ('rain_i');  um_rain(:,:,:) = vwork(:,:,:)
!    case ('snow_i');  um_snow(:,:,:) = vwork(:,:,:)
!---20100825 -- just be cauious: -------------------------
    case ('rain_i');  um_rain(:,:,:) = max(0.0,vwork(:,:,:))
    case ('snow_i');  um_snow(:,:,:) = max(0.0,vwork(:,:,:))
!---------------------------------------------------------   
    case ('evap_i');  um_evap(:,:,:) = vwork(:,:,:)
    case ('lhflx_i'); um_lhflx(:,:,:) = vwork(:,:,:)
    case ('tmlt01_i'); um_tmlt(:,:,1,:) = vwork(:,:,:) 
    case ('tmlt02_i'); um_tmlt(:,:,2,:) = vwork(:,:,:)
    case ('tmlt03_i'); um_tmlt(:,:,3,:) = vwork(:,:,:)
    case ('tmlt04_i'); um_tmlt(:,:,4,:) = vwork(:,:,:)
    case ('tmlt05_i'); um_tmlt(:,:,5,:) = vwork(:,:,:)
    case ('bmlt01_i'); um_bmlt(:,:,1,:) = vwork(:,:,:)
    case ('bmlt02_i'); um_bmlt(:,:,2,:) = vwork(:,:,:)
    case ('bmlt03_i'); um_bmlt(:,:,3,:) = vwork(:,:,:)
    case ('bmlt04_i'); um_bmlt(:,:,4,:) = vwork(:,:,:)
    case ('bmlt05_i'); um_bmlt(:,:,5,:) = vwork(:,:,:)
    case ('taux_i');  um_taux(:,:,:)  = vwork(:,:,:)
    case ('tauy_i');  um_tauy(:,:,:)  = vwork(:,:,:)
    case ('swflx_i'); um_swflx(:,:,:) = vwork(:,:,:)
    case ('lwflx_i'); um_lwflx(:,:,:) = vwork(:,:,:)
    case ('shflx_i'); um_shflx(:,:,:) = vwork(:,:,:)
    case ('press_i'); um_press(:,:,:) = vwork(:,:,:)
    case ('co2_ai'); um_co2(:,:,:) = vwork(:,:,:)
    case ('wnd_ai'); um_wnd(:,:,:) = vwork(:,:,:)
    end select 

    if (my_task == 0 .or. ll_comparal) then
      write(il_out,*) 
      write(il_out,*)'(from_atm) done: ', jf, trim(cl_read(jf))
    endif

  enddo

  IF (rotate_winds) THEN   !rotate_winds=.t. means oasis does not do the vector rotation.

  do iblk = 1, nblocks

    this_block = get_block(blocks_ice(iblk),iblk)
    ilo = this_block%ilo
    ihi = this_block%ihi
    jlo = this_block%jlo
    jhi = this_block%jhi

    do j = jlo, jhi
    do i = ilo, ihi
      tmpu = um_taux(i,j,iblk)   ! on geographical coord. (T cell)
      tmpv = um_tauy(i,j,iblk)    
      um_taux(i,j,iblk) = tmpu*cos(ANGLET(i,j,iblk)) &   ! converted onto model curvelear
                   + tmpv*sin(ANGLET(i,j,iblk))          ! coord. (T cell)
      um_tauy(i,j,iblk) = tmpv*cos(ANGLET(i,j,iblk)) &   ! 
                   - tmpu*sin(ANGLET(i,j,iblk))
    enddo
    enddo

  enddo

  ENDIF  !rotate_winds

  ! need do t-grid to u-grid shift for vectors since all coupling occur on
  ! t-grid points: <==No! actually CICE requires the input wind on T grid! 
  ! (see comment in code ice_flux.F)
  !call t2ugrid(uwnd1)
  !call t2ugrid(vwnd1)

  !-------------------------------
  !if ( chk_a2i_fields ) then
  !  call check_a2i_fields(isteps)
  !endif
  !-------------------------------

  if ( chk_a2i_fields .and. my_task == 0 ) then
    call ncheck(nf_close(ncid))
  endif

  end subroutine from_atm

!=======================================================================
  subroutine from_ocn(isteps)
!----------------------------!  

  integer(kind=int_kind), intent(in) :: isteps
 
  integer(kind=int_kind) :: jf
  integer(kind=int_kind) :: ilo,ihi,jlo,jhi,iblk,i,j
  type (block) :: this_block           ! block information for current block

  integer(kind=int_kind) :: ncid,currstep,ll,ilout
  data currstep/0/
  save currstep
  
  currstep=currstep+1

  if (my_task == 0) then
    write(il_out,*) '(from_ocn) receiving coupling fields at rtime: ', isteps
    if (chk_o2i_fields) then
      if ( .not. file_exist('fields_o2i_in_ice.nc') ) then
        call create_ncfile('fields_o2i_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
      endif
      write(il_out,*) 'opening file fields_o2i_in_ice.nc at nstep = ', isteps
      call ncheck( nf_open('fields_o2i_in_ice.nc',nf_write,ncid) )
      call write_nc_1Dtime(real(isteps),currstep,'time',ncid)
    endif
  endif

  write(il_out,*) "prism_get from_ocn at sec: ", isteps
  do jf = nrecv_a2i + 1, jpfldin 
  
    if (my_task==0 .or. ll_comparal) then

      !jf-th field in
      write(il_out,*)
      write(il_out,*) '*** receiving coupling fields No. ', jf, cl_read(jf)
      if(ll_comparal) then
        call prism_get_proto (il_var_id_in(jf), isteps, vwork2d(l_ilo:l_ihi, l_jlo:l_jhi), ierror)
        call mpi_gatherv(vwork2d(l_ilo:l_ihi, l_jlo:l_jhi),1,sendsubarray,gwork,counts,disps,resizedrecvsubarray, &
                     0,MPI_COMM_ICE,ierror)
        call broadcast_array(gwork, 0)
!         gwork(l_ilo:l_ihi, l_jlo:l_jhi) = vwork2d(l_ilo:l_ihi, l_jlo:l_jhi)
      else
        call prism_get_proto (il_var_id_in(jf), isteps, gwork, ierror)
      endif 

      if ( ierror /= PRISM_Ok .and. ierror < PRISM_Recvd) then
        write(il_out,*) 'Err in _get_ sst at time with error: ', isteps, ierror
        call prism_abort_proto(il_comp_id, 'cice from_ocn','stop 1')
      else
        write(il_out,*)
        write(il_out,*)'(from_ocn) rcvd at time with err: ',cl_read(jf),isteps,ierror
          if(ll_comparal .and. chk_o2i_fields) then
            call mpi_gatherv(vwork2d(l_ilo:l_ihi, l_jlo:l_jhi),1,sendsubarray,gwork,  &
                       counts,disps,resizedrecvsubarray, 0,MPI_COMM_ICE,ierror)
           call MPI_Barrier(MPI_COMM_ICE, ierror)
          endif
        if (my_task==0 .and. chk_o2i_fields) then
          call write_nc2D(ncid, trim(cl_read(jf)), gwork, 1, il_im,il_jm, &
                          currstep,ilout=il_out)
        endif
       endif

    endif
    
    if (.not. ll_comparal) then
      call scatter_global(vwork, gwork, master_task, distrb_info, &
                        field_loc_center, field_type_scalar)
    else
      call unpack_global_dbl(vwork, gwork, master_task, distrb_info, &
                        field_loc_center, field_type_scalar)
    endif 

    !Q: 'field_type_scalar' all right for 'vector' (ssu/ssv, sslx/ssly))?! 
    select case (trim(cl_read(jf)))
    case ('sst_i'); ocn_sst = vwork
    case ('sss_i'); ocn_sss = vwork 
    case ('ssu_i'); ocn_ssu = vwork
    case ('ssv_i'); ocn_ssv = vwork
    case ('sslx_i'); ocn_sslx = vwork
    case ('ssly_i'); ocn_ssly = vwork
    case ('pfmice_i'); ocn_pfmice = vwork
    case ('co2_oi'); ocn_co2 = vwork
    case ('co2fx_oi'); ocn_co2fx = vwork
    end select

  enddo

  !-------------------------------
  !if (chk_o2i_fields) then
  !  call check_o2i_fields(isteps)
  !endif
  !-------------------------------

  if ( chk_o2i_fields .and. my_task == 0 ) then
    call ncheck(nf_close(ncid))
  endif

  end subroutine from_ocn

!=======================================================================
  subroutine into_ocn(isteps)
!-----------------------------------!

  implicit none  
 
  integer(kind=int_kind), intent(in) :: isteps
  integer(kind=int_kind) :: jf
  integer(kind=int_kind) :: ilo,ihi,jlo,jhi,iblk,i,j
  type (block) :: this_block           ! block information for current block

  integer(kind=int_kind) :: ncid,currstep,ll,ilout
  data currstep/0/
  save currstep

  currstep=currstep+1

  if (my_task == 0) then  
    write(il_out,*)
    write(il_out,*) '(into_ocn) sending coupling fields at stime= ', isteps
    if (chk_i2o_fields) then
      if ( .not. file_exist('fields_i2o_in_ice.nc') ) then
        call create_ncfile('fields_i2o_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
      endif
      write(il_out,*) 'opening file fields_i2o_in_ice.nc at nstep = ', isteps
      call ncheck( nf_open('fields_i2o_in_ice.nc',nf_write,ncid) )
      call write_nc_1Dtime(real(isteps),currstep,'time',ncid)
    endif
  endif

  write(il_out,*) "prism_put into_ocn at sec: ", isteps
  do jf = nsend_i2a + 1, jpfldout

!CH: make sure the 'LIMITS' are to be released!

    select case(trim(cl_writ(jf)))
!20100531 for MYM's test (iostress_factor) .............
    case('strsu_io'); vwork = io_strsu * iostress_factor
    case('strsv_io'); vwork = io_strsv * iostress_factor
!.......................................................
    case('rain_io'); vwork = io_rain
    case('snow_io'); vwork = io_snow
    !case('stflx_io'); vwork = io_stflx
    case('stflx_io')
        if (limit_stflx) then
          vwork = max(-5.e-6, min(io_stflx, 5.e-6))
        else
          vwork = io_stflx
        endif
    !case('htflx_io'); vwork = io_htflx
    !case('htflx_io'); vwork = max(io_htflx, -450.0)
    !Jan2010:
    case('htflx_io'); vwork = min(io_htflx,0.0)
    case('swflx_io'); vwork = io_swflx
    case('qflux_io'); vwork = io_qflux
    case('shflx_io'); vwork = io_shflx
    case('lwflx_io'); vwork = io_lwflx
    case('runof_io')
       if (use_core_runoff) then 
         vwork = core_runoff
       else 
         vwork = io_runof
       endif 
    case('press_io'); vwork = io_press
    case('aice_io');  vwork = io_aice
    case('melt_io');  vwork = io_melt
    case('form_io');  vwork = io_form
    case('co2_i1'); vwork = io_co2
    case('wnd_i1'); vwork = io_wnd
    end select
    
    if(.not. ll_comparal) then 
      call gather_global(gwork, vwork, master_task, distrb_info)
    else
      call pack_global_dbl(gwork, vwork, master_task, distrb_info)
      vwork2d(l_ilo:l_ihi, l_jlo:l_jhi) = gwork(l_ilo:l_ihi, l_jlo:l_jhi)
!      do iblk=1,nblocks_tot
!
!        if (distrb_info%blockLocation(iblk) == my_task+1) then
!
!          this_block = get_block(iblk,iblk)
!          ilo = this_block%ilo
!          ihi = this_block%ihi
!          jlo = this_block%jlo
!          jhi = this_block%jhi
!
!          vwork2d(this_block%i_glob(ilo):this_block%i_glob(ihi),   &
!             this_block%j_glob(jlo):this_block%j_glob(jhi)) =      &
!          vwork(ilo:ihi,jlo:jhi,distrb_info%blockLocalID(iblk)) 
!        endif
!      end do
      
    endif
    if (my_task == 0 .or. ll_comparal) then   

      write(il_out,*)
      write(il_out,*) '*** sending coupling field No. ', jf, cl_writ(jf)
      if(ll_comparal) then
        call prism_put_proto(il_var_id_out(jf), isteps, vwork2d(l_ilo:l_ihi, l_jlo:l_jhi), ierror)
      else
        call prism_put_proto(il_var_id_out(jf), isteps, gwork, ierror)
      endif
      if ( ierror /= PRISM_Ok .and. ierror < PRISM_Sent) then
        write(il_out,*)
        write(il_out,*) '(into_ocn) Err in _put_ ', cl_writ(jf), isteps, ierror
        call prism_abort_proto(il_comp_id, 'cice into_ocn','STOP 1') 
      else
        write(il_out,*)
        write(il_out,*)'(into_ocn) sent: ', cl_writ(jf), isteps, ierror
        if(chk_i2o_fields .and. ll_comparal) then
          call mpi_gatherv(vwork2d(l_ilo:l_ihi, l_jlo:l_jhi),1,sendsubarray,gwork, &
                     counts,disps,resizedrecvsubarray, 0,MPI_COMM_ICE,ierror)
           call MPI_Barrier(MPI_COMM_ICE, ierror)
        endif
        if (my_task==0 .and. chk_i2o_fields) then
            call write_nc2D(ncid, trim(cl_writ(jf)), gwork, 1, il_im,il_jm, &
                          currstep,ilout=il_out)
        endif
      endif

    endif !my_task == 0

  enddo   !jf 

  !--------------------------------------
  !if (chk_i2o_fields) then
  !  call check_i2o_fields(isteps)
  !endif
  !--------------------------------------

  if ( chk_i2o_fields .and. my_task == 0 ) then
    call ncheck(nf_close(ncid))
  endif

  end subroutine into_ocn

!=======================================================================
  subroutine into_atm(isteps)
!----------------------------!    

  integer(kind=int_kind), intent(in) :: isteps
  integer(kind=int_kind) :: jf

  real(kind=dbl_kind) :: tmpu, tmpv
  integer(kind=int_kind) :: ilo,ihi,jlo,jhi,iblk,i,j
  type (block) :: this_block           ! block information for current block

  integer(kind=int_kind) :: ncid,currstep,ll,ilout
  data currstep/0/
  save currstep

  currstep=currstep+1

!debug hxy599
!if (isteps==runtime-3600) then
! chk_i2a_fields=.true. !save the last step
! currstep = 1
! write (il_out,*) 'hxy599 debug: save i2a fields at time ', isteps
!end if

  if (my_task == 0) then
    write(il_out,*)
    write(il_out,*) '(into_atm) sending coupling fields at stime= ', isteps
    if (chk_i2a_fields) then
      if ( .not. file_exist('fields_i2a_in_ice.nc') ) then
        call create_ncfile('fields_i2a_in_ice.nc',ncid,il_im,il_jm,ll=1,ilout=il_out)
      else
        write(il_out,*) 'opening file fields_i2a_in_ice.nc at nstep = ', isteps
        call ncheck( nf_open('fields_i2a_in_ice.nc',nf_write,ncid) )
      end if
      call write_nc_1Dtime(real(isteps),currstep,'time',ncid)
    endif
  endif

  IF (rotate_winds) THEN

  do iblk = 1, nblocks

    this_block = get_block(blocks_ice(iblk),iblk)
    ilo = this_block%ilo
    ihi = this_block%ihi
    jlo = this_block%jlo
    jhi = this_block%jhi

    do j = jlo, jhi
    do i = ilo, ihi
      !note note uvel/vvel are on the U-cell here.
      tmpu = ia_uvel(i,j,iblk); tmpv = ia_vvel(i,j,iblk)   ! ice/ocn velocity, m/s
      ia_uvel(i,j,iblk) = tmpu*cos(ANGLE(i,j,iblk)) &      ! remapped on to geographical
                   - tmpv*sin(ANGLE(i,j,iblk))             ! grid. 
      ia_vvel(i,j,iblk) = tmpv*cos(ANGLE(i,j,iblk)) &      ! they also need be shifted 
                   + tmpu*sin(ANGLE(i,j,iblk))             ! on to T-cell (below).
    enddo
    enddo

  enddo 

  ENDIF  !rotate_winds

  !shift ia_uvel/ia_vvel onto T points before passing into coupler
  call u2tgrid_vector(ia_uvel)
  call u2tgrid_vector(ia_vvel) 

  write(il_out,*) "prism_put into_atm at sec: ", isteps
  do jf = 1, nsend_i2a

    select case (trim(cl_writ(jf)))
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
    case('icethk04'); vwork(:,:,:) = ia_thikn(:,:,4,:)
    case('icethk05'); vwork(:,:,:) = ia_thikn(:,:,5,:)
    !20100305: test effect of ssuv on the tropical cooling biases (as per Harry Henden)
    case('uvel_ia');  vwork = ia_uvel * ocn_ssuv_factor     !note ice u/v are also 
    case('vvel_ia');  vwork = ia_vvel * ocn_ssuv_factor     !     included here.
    case('co2_i2');  vwork = ia_co2
    case('co2fx_i2');  vwork = ia_co2fx
    end select
    
    if (.not. ll_comparal) then
      call gather_global(gwork, vwork, master_task, distrb_info)
    else
      call pack_global_dbl(gwork, vwork, master_task, distrb_info)
      vwork2d(l_ilo:l_ihi, l_jlo:l_jhi) = gwork(l_ilo:l_ihi, l_jlo:l_jhi)
!      do iblk=1,nblocks_tot
!
!        if (distrb_info%blockLocation(iblk) == my_task+1) then
!
!          this_block = get_block(iblk,iblk)
!          ilo = this_block%ilo
!          ihi = this_block%ihi
!          jlo = this_block%jlo
!          jhi = this_block%jhi
!
!          vwork2d(this_block%i_glob(ilo):this_block%i_glob(ihi),   &
!             this_block%j_glob(jlo):this_block%j_glob(jhi)) =      &
!          vwork(ilo:ihi,jlo:jhi,distrb_info%blockLocalID(iblk))
!        endif
!      end do

    end if
    if (my_task == 0 .or. ll_comparal) then
  
      write(il_out,*)
      write(il_out,*) '*** sending coupling field No. ', jf, cl_writ(jf)

      !call prism_put_inquire_proto(il_var_id_out(jf),isteps,ierror)
  
      write(il_out,*)
      write(il_out,*) '(into_atm) what to do with this var==> Err= ',ierror
      if(ll_comparal) then 
        call prism_put_proto(il_var_id_out(jf), isteps, vwork2d(l_ilo:l_ihi, l_jlo:l_jhi), ierror)
      else
        call prism_put_proto(il_var_id_out(jf), isteps, gwork, ierror)
      endif 
      if ( ierror /= PRISM_Ok .and. ierror < PRISM_Sent) then
        write(il_out,*)
        write(il_out,*) '(into_atm) Err in _put_ ', cl_writ(jf), isteps, ierror
        call prism_abort_proto(il_comp_id, 'cice into_atm','STOP 1')
      else
        write(il_out,*)
        write(il_out,*)'(into_atm) sent: ', cl_writ(jf), isteps, ierror
        if(chk_i2a_fields .and. ll_comparal) then
          call mpi_gatherv(vwork2d(l_ilo:l_ihi, l_jlo:l_jhi),1,sendsubarray,gwork, &
                     counts,disps,resizedrecvsubarray, 0,MPI_COMM_ICE,ierror)
           call MPI_Barrier(MPI_COMM_ICE, ierror)
        endif
        if (my_task==0 .and. chk_i2a_fields) then
            call write_nc2D(ncid, trim(cl_writ(jf)), gwork, 1, il_im,il_jm, &
                          currstep,ilout=il_out)
        endif
      endif

    endif !my_task == 0

  enddo     !jf = 1, jpfldout

  !-------------------------------
  !if (chk_i2a_fields) then
  !  call check_i2a_fields(isteps)
  !endif 
  !-------------------------------

  if ( chk_i2a_fields .and. my_task == 0 ) then
    call ncheck(nf_close(ncid))
  endif

  end subroutine into_atm

!=======================================================================
  subroutine coupler_termination
!-------------------------------!
  !
  ! Detach from MPI buffer
  !
  call MPI_Buffer_Detach(rla_bufsend, il_bufsize, ierror) 
  deallocate (rla_bufsend)
  !deallocate all the coupling associated arrays... (no bother...) 
  !  
  ! 9- PSMILe termination 
  !   

  call MPI_Barrier(MPI_COMM_ICE, ierror)
  call prism_terminate_proto (ierror)
  if (ierror /= PRISM_Ok) then
    if (my_task == 0) then
      write (il_out,*) 'An error occured in prism_terminate = ', ierror
    endif
  else 
    if (my_task == 0) then
      write(il_out,*)
      write(il_out,*) '==================*** END ***================='
      close(il_out)
    endif
  endif
  ! 
  print *
  print *, '********** End of CICE **********'
  print *

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
          
          !ilo = 1 + nghost
          !ihi = nx_block - nghost
          !jlo = 1 + nghost
          !jhi = ny_block - nghost
          
          il_paral ( clim_strategy ) = clim_Box
          il_paral ( clim_offset   ) = nx_global * (l_jlo-1) + (l_ilo-1)
          !il_paral ( clim_offset   ) = (l_ilo-1)
          il_paral ( clim_SizeX    ) = l_ihi-l_ilo+1
          il_paral ( clim_SizeY    ) = l_jhi-l_jlo+1
          il_paral ( clim_LdX      ) = nx_global
          
          write(ld_mparout,*)'il_paral=',il_paral
 
          id_length = il_paral(clim_sizeX) * il_paral(clim_sizeY)
          
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

!============================================================================

!============================================================================
 subroutine unpack_global_dbl(ARRAY, ARRAY_G, src_task, dst_dist, &
                               field_loc, field_type)

! !DESCRIPTION:
!  This subroutine scatters a global-sized array to a distributed array.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific interface for double precision arrays
!  corresponding to the generic interface scatter_global.

! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     src_task       ! task from which array should be scattered

   type (distrb), intent(in) :: &
     dst_dist       ! distribution of resulting blocks

   real (dbl_kind), dimension(:,:), intent(in) :: &
     ARRAY_G        ! array containing global field on src_task

   integer (int_kind), intent(in) :: &
      field_type,               &! id for type of field (scalar, vector, angle)
      field_loc                  ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

! !OUTPUT PARAMETERS:

   real (dbl_kind), dimension(:,:,:), intent(inout) :: &
     ARRAY          ! array containing distributed field
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n,bid,          &! dummy loop indices
     nrecvs,             &! actual number of messages received
     isrc, jsrc,         &! source addresses
     dst_block,          &! location of block in dst array
     xoffset, yoffset,   &! offsets for tripole boundary conditions
     yoffset2,           &!
     isign,              &! sign factor for tripole boundary conditions
     ierr                 ! MPI error flag

   type (block) :: &
     this_block  ! block info for current block

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     rcv_request     ! request array for receives

   integer (int_kind), dimension(:,:), allocatable :: &
     rcv_status      ! status array for receives

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer      ! buffer for sending blocks

!-----------------------------------------------------------------------
!
!  initialize return array to zero and set up tripole quantities
!
!-----------------------------------------------------------------------
   ARRAY = c0

   this_block = get_block(1,1) ! for the tripoleTflag - all blocks have it
   if (this_block%tripoleTFlag) then
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 2
        yoffset = 0
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 1
        yoffset = -1
     case (field_loc_Eface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_Nface)    ! cell face location
        xoffset = 2
        yoffset = -1
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   else
     select case (field_loc)
     case (field_loc_center)   ! cell center location
        xoffset = 1
        yoffset = 1
     case (field_loc_NEcorner) ! cell corner (velocity) location
        xoffset = 0
        yoffset = 0
     case (field_loc_Eface)    ! cell face location
        xoffset = 0
        yoffset = 1
     case (field_loc_Nface)    ! cell face location
        xoffset = 1
        yoffset = 0
     case (field_loc_noupdate) ! ghost cells never used - use cell center
        xoffset = 1
        yoffset = 1
     end select
   endif
   select case (field_type)
   case (field_type_scalar)
      isign =  1
   case (field_type_vector)
      isign = -1
   case (field_type_angle)
      isign = -1
   case (field_type_noupdate) ! ghost cells never used - use cell center
      isign =  1
   case default
      call abort_ice('Unknown field type in scatter')
   end select


     !*** copy any local blocks

     do n=1,nblocks_tot
       if (dst_dist%blockLocation(n) == my_task+1) then
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         !*** if this is an interior block, then there is no
         !*** padding or update checking required

         if (this_block%iblock > 1         .and. &
             this_block%iblock < nblocks_x .and. &
             this_block%jblock > 1         .and. &
             this_block%jblock < nblocks_y) then

!            do j=1,ny_block
!            do i=1,nx_block
!               ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
!                                              this_block%j_glob(j))
!            end do
!            end do
            ARRAY(1:nx_block,1:ny_block,dst_block) =                              &
                        ARRAY_G(this_block%i_glob(1):this_block%i_glob(nx_block), &
                                this_block%j_glob(1):this_block%j_glob(ny_block))

         !*** if this is an edge block but not a northern edge
         !*** we only need to check for closed boundaries and
         !*** padding (global index = 0)

         else if (this_block%jblock /= nblocks_y) then

            do j=1,ny_block
               if (this_block%j_glob(j) /= 0) then
                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do
               endif
            end do

         !*** if this is a northern edge block, we need to check
         !*** for and properly deal with tripole boundaries

         else

            do j=1,ny_block
               if (this_block%j_glob(j) > 0) then ! normal boundary

                  do i=1,nx_block
                     if (this_block%i_glob(i) /= 0) then
                        ARRAY(i,j,dst_block) = ARRAY_G(this_block%i_glob(i),&
                                                       this_block%j_glob(j))
                     endif
                  end do

               else if (this_block%j_glob(j) < 0) then  ! tripole
                  ! for yoffset=0 or 1, yoffset2=0,0
                  ! for yoffset=-1, yoffset2=0,1, for u-rows on T-fold grid
                  do yoffset2=0,max(yoffset,0)-yoffset
                    jsrc = ny_global + yoffset + yoffset2 + &
                         (this_block%j_glob(j) + ny_global)
                    do i=1,nx_block
                      if (this_block%i_glob(i) /= 0) then
                         isrc = nx_global + xoffset - this_block%i_glob(i)
                         if (isrc < 1) isrc = isrc + nx_global
                         if (isrc > nx_global) isrc = isrc - nx_global
                         ARRAY(i,j-yoffset2,dst_block) &
                           = isign * ARRAY_G(isrc,jsrc)
                      endif
                    end do
                  end do

               endif
            end do

         endif
       endif
     end do

   !-----------------------------------------------------------------
   ! Ensure unused ghost cell values are 0
   !-----------------------------------------------------------------

   if (field_loc == field_loc_noupdate) then
      do n=1,nblocks_tot
         dst_block = dst_dist%blockLocalID(n)
         this_block = get_block(n,n)

         if (dst_block > 0) then

         ! north edge
         do j = this_block%jhi+1,ny_block
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! east edge
         do j = 1, ny_block
         do i = this_block%ihi+1,nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! south edge
         do j = 1, this_block%jlo-1
         do i = 1, nx_block
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo
         ! west edge
         do j = 1, ny_block
         do i = 1, this_block%ilo-1
            ARRAY (i,j,dst_block) = c0
         enddo
         enddo

         endif
      enddo
   endif

 end subroutine unpack_global_dbl
!============================================================================

!============================================================================
 subroutine pack_global_dbl(ARRAY_G, ARRAY, dst_task, src_dist)

! !DESCRIPTION:
!  This subroutine gathers a distributed array to a global-sized
!  array on the processor dst_task.
!
! !REVISION HISTORY:
!  same as module
!
! !REMARKS:
!  This is the specific inteface for double precision arrays
!  corresponding to the generic interface gather_global.  It is shown
!  to provide information on the generic interface (the generic
!  interface is identical, but chooses a specific inteface based
!  on the data type of the input argument).


! !USES:

   include 'mpif.h'

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     dst_task   ! task to which array should be gathered

   type (distrb), intent(in) :: &
     src_dist   ! distribution of blocks in the source array

   real (dbl_kind), dimension(:,:,:), intent(in) :: &
     ARRAY      ! array containing horizontal slab of distributed field

! !OUTPUT PARAMETERS:


   real (dbl_kind), dimension(:,:), intent(inout) :: &
     ARRAY_G    ! array containing global horizontal field on dst_task

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,n          ,&! dummy loop counters
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (dbl_kind), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block  ! block info for current block

     do n=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%blockLocation(n) == my_task+1) then

         this_block = get_block(n,n)

!         do j=this_block%jlo,this_block%jhi
!         do i=this_block%ilo,this_block%ihi
!           ARRAY_G(this_block%i_glob(i), &
!                   this_block%j_glob(j)) = &
!                  ARRAY(i,j,src_dist%blockLocalID(n))
!         end do
!         end do
          ARRAY_G(this_block%i_glob(this_block%ilo):this_block%i_glob(this_block%ihi), &
                 this_block%j_glob(this_block%jlo):this_block%j_glob(this_block%jhi)) = &
                 ARRAY(this_block%ilo:this_block%ihi,this_block%jlo:this_block%jhi,src_dist%blockLocalID(n))

       !*** fill land blocks with special values

       else if (src_dist%blockLocation(n) == 0) then

         this_block = get_block(n,n)

!         do j=this_block%jlo,this_block%jhi
!         do i=this_block%ilo,this_block%ihi
!           ARRAY_G(this_block%i_glob(i), &
!                   this_block%j_glob(j)) = spval_dbl
!         end do
!         end do
         ARRAY_G(this_block%i_glob(this_block%ilo):this_block%i_glob(this_block%ihi), &
                 this_block%j_glob(this_block%jlo):this_block%j_glob(this_block%jhi)) = spval_dbl
       endif

     end do

  end subroutine pack_global_dbl
!============================================================================

!==============================================================================
subroutine save_restart_i2a(fname, nstep)
! output the last i2a forcing data in cice by the end of the run,
! to be read in at the beginning of next run by cice and sent to atm

implicit none

character*(*), intent(in) :: fname
integer(kind=int_kind), intent(in) :: nstep
integer(kind=int_kind) :: ncid
integer(kind=int_kind) :: jf, jfs, ll, ilout

if (my_task == 0) then
  call open_ncfile(fname, ncid, il_im, il_jm, ll=1, ilout=il_out)
endif

do jf = 1, nsend_i2a
    select case (trim(cl_writ(jf)))
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
    case('icethk04'); vwork(:,:,:) = ia_thikn(:,:,4,:)
    case('icethk05'); vwork(:,:,:) = ia_thikn(:,:,5,:)
    !20100305: test effect of ssuv on the tropical cooling biases (as per Harry Henden)
    case('uvel_ia');  vwork = ia_uvel !* ocn_ssuv_factor     !note ice u/v are also
    case('vvel_ia');  vwork = ia_vvel !* ocn_ssuv_factor     !     included here.
    end select

!    if (.not. ll_comparal) then
      call gather_global(gwork, vwork, master_task, distrb_info)
!    else
!      call pack_global_dbl(gwork, vwork, master_task, distrb_info)
!    end if
  if (my_task == 0) then
    call modify_nc2D(ncid, cl_writ(jf), gwork, 2, il_im, il_jm, 1, ilout=il_out)
  endif

enddo

if (my_task == 0) call ncheck( nf_close(ncid) )

end subroutine save_restart_i2a
!==========================================================

  end module cpl_interface
