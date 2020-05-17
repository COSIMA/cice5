!============================================================================
!  Writes netcdf files
!    Created by Mariana Vertenstein, June 2009

  module ice_pio

  use ice_kinds_mod
  use ice_blocks
  use ice_broadcast
  use ice_communicate, only : master_task, my_task, MPI_COMM_ICE, get_num_procs
  use ice_domain, only : nblocks, blocks_ice
  use ice_domain_size
  use ice_fileunits
  use ice_exit
  use pio
  use pio, only: pio_set_buffer_size_limit
  use pio_types, only: pio_iotype_netcdf4p, PIO_rearr_box

  implicit none

  private
  public :: ice_pio_subsystem
  save

  interface ice_pio_initdecomp
     module procedure ice_pio_initdecomp_2d
     module procedure ice_pio_initdecomp_3d
     module procedure ice_pio_initdecomp_4d
     module procedure ice_pio_initdecomp_3d_inner
  end interface

  public ice_pio_init
  public ice_pio_initfile
  public ice_pio_initdecomp
  logical :: pio_initialized = .false.
  integer :: pio_iotype
  type(iosystem_desc_t) :: ice_pio_subsystem

!===============================================================================

  contains

!===============================================================================

   subroutine ice_pio_init(io_stride)
        integer, intent(in), optional :: io_stride
        integer :: num_iotasks, stride, ierr

        character(*),parameter :: subName = '(ice_pio_init) '

        if (pio_initialized) then
            return
        endif

        if (present(io_stride)) then
            stride = io_stride
        else
            stride = 1
        endif

        pio_iotype = pio_iotype_netcdf4p

        num_iotasks = get_num_procs() / stride

        call pio_init(my_task, MPI_COMM_ICE, num_iotasks, 0, stride, PIO_rearr_box, ice_pio_subsystem)

        call pio_set_buffer_size_limit(256*1024*1024)

        pio_initialized = .true.
   end subroutine ice_pio_init

   subroutine ice_pio_initfile(mode, filename, File, clobber, cdf64)


   implicit none
   character(len=*)     , intent(in),    optional :: mode
   character(len=*)     , intent(in),    optional :: filename
   type(file_desc_t)    , intent(inout), optional :: File
   logical              , intent(in),    optional :: clobber
   logical              , intent(in),    optional :: cdf64

   ! local variables

   integer (int_kind) :: &
      nml_error          ! namelist read error flag

   logical :: exists
   logical :: lclobber
   integer :: status
   character(*),parameter :: subName = '(ice_pio_wopen) '

   if (present(mode) .and. present(filename) .and. present(File)) then

      if (trim(mode) == 'write') then
         lclobber = .false.
         if (present(clobber)) lclobber=clobber

         if (File%fh<0) then
            ! filename not open
            inquire(file=trim(filename),exist=exists)
            if (exists) then
               if (lclobber) then
                  status = pio_createfile(ice_pio_subsystem, File, pio_iotype, trim(filename), PIO_clobber)
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,' create file ',trim(filename)
                  end if
               else
                  status = pio_openfile(ice_pio_subsystem, File, pio_iotype, trim(filename), pio_write)
                  if (my_task == master_task) then
                     write(nu_diag,*) subname,' open file ',trim(filename)
                  end if
               endif
            else
               status = pio_createfile(ice_pio_subsystem, File, pio_iotype, trim(filename), pio_noclobber)
               if (my_task == master_task) then
                  write(nu_diag,*) subname,' create file ',trim(filename)
               end if
            endif
         else
            ! filename is already open, just return
         endif
      end if

      if (trim(mode) == 'read') then
         inquire(file=trim(filename),exist=exists)
         if (exists) then
            status = pio_openfile(ice_pio_subsystem, File, pio_iotype, trim(filename), pio_nowrite)
         else
            if(my_task==master_task) then
               write(nu_diag,*) 'ice_pio_ropen ERROR: file invalid ',trim(filename)
            end if
            call abort_ice('aborting in ice-pio_ropen with invalid file')
         endif
      end if

   end if

   end subroutine ice_pio_initfile

!================================================================================

   subroutine ice_pio_initdecomp_2d(iodesc, use_double)

      type(io_desc_t), intent(out) :: iodesc
      logical, intent(in), optional :: use_double

      logical :: luse_double
      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k

      type(block) :: this_block

      integer(kind=int_kind), pointer :: dof2d(:)

      allocate(dof2d(nx_block*ny_block*nblocks))

      luse_double = .false.
      if (present(use_double)) then
          luse_double = .true.
      endif

      n=0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         do j=1,ny_block
         do i=1,nx_block
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof2d(n) = 0
            else if (i < ilo .or. i > ihi) then
               dof2d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof2d(n) = (lat-1)*nx_global + lon
            endif
         enddo !i
         enddo !j
      end do

      if (luse_double) then
          call pio_initdecomp(ice_pio_subsystem, pio_double, (/nx_global,ny_global/), &
               dof2d, iodesc)
      else
          call pio_initdecomp(ice_pio_subsystem, pio_real, (/nx_global,ny_global/), &
               dof2d, iodesc)
      endif

      deallocate(dof2d)
 
   end subroutine ice_pio_initdecomp_2d

!================================================================================

   subroutine ice_pio_initdecomp_3d (ndim3, iodesc, remap, use_double)

      integer(kind=int_kind), intent(in) :: ndim3
      type(io_desc_t), intent(out) :: iodesc
      logical, optional :: remap
      logical, intent(in), optional :: use_double
      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k 

      type(block) :: this_block 
      logical :: lremap, luse_double
      integer(kind=int_kind), pointer :: dof3d(:)

      allocate(dof3d(nx_block*ny_block*nblocks*ndim3))

      lremap = .false.
      if (present(remap)) then
          lremap = remap
      endif

      luse_double = .false.
      if (present(use_double)) then
          luse_double = use_double
      endif

      if(lremap) then
         ! Reorder the ndim3 and nblocks loops to avoid a temporary array in restart read/write
         n=0
         do iblk = 1, nblocks
            this_block = get_block(blocks_ice(iblk),iblk)         
            ilo = this_block%ilo
            ihi = this_block%ihi
            jlo = this_block%jlo
            jhi = this_block%jhi
            do k=1,ndim3         
               do j=1,ny_block
                  do i=1,nx_block
                     n = n+1
                     if (j < jlo .or. j>jhi) then
                        dof3d(n)=0
                     else if (i < ilo .or. i > ihi) then
                        dof3d(n) = 0
                     else
                        lon = this_block%i_glob(i)
                        lat = this_block%j_glob(j)
                        dof3d(n) = ((lat-1)*nx_global + lon) + (k-1)*nx_global*ny_global 
                     endif
                  enddo !i
               enddo !j
            enddo !ndim3
         enddo ! iblk
   else
         n=0
         do k=1,ndim3         
            do iblk = 1, nblocks
               this_block = get_block(blocks_ice(iblk),iblk)         
               ilo = this_block%ilo
               ihi = this_block%ihi
               jlo = this_block%jlo
               jhi = this_block%jhi
               do j=1,ny_block
                  do i=1,nx_block
                     n = n+1
                     if (j < jlo .or. j>jhi) then
                        dof3d(n)=0
                     else if (i < ilo .or. i > ihi) then
                        dof3d(n) = 0
                     else
                        lon = this_block%i_glob(i)
                        lat = this_block%j_glob(j)
                        dof3d(n) = ((lat-1)*nx_global + lon) + (k-1)*nx_global*ny_global 
                     endif
                  enddo !i
               enddo !j
            enddo ! iblk
         enddo !ndim3
      endif

      if (luse_double) then
          call pio_initdecomp(ice_pio_subsystem, pio_double, (/nx_global,ny_global,ndim3/), &
               dof3d, iodesc)
      else
          call pio_initdecomp(ice_pio_subsystem, pio_real, (/nx_global,ny_global,ndim3/), &
               dof3d, iodesc)
      endif

      deallocate(dof3d)

   end subroutine ice_pio_initdecomp_3d

!================================================================================

   subroutine ice_pio_initdecomp_3d_inner(ndim3, inner_dim, iodesc)

      integer(kind=int_kind), intent(in) :: ndim3
      logical, intent(in) :: inner_dim
      type(io_desc_t), intent(out) :: iodesc

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k 

      type(block) :: this_block 

      integer(kind=int_kind), pointer :: dof3d(:)

      allocate(dof3d(nx_block*ny_block*nblocks*ndim3))

      n=0
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         
         do j=1,ny_block
         do i=1,nx_block
         do k=1,ndim3
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof3d(n) = 0
            else if (i < ilo .or. i > ihi) then
               dof3d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof3d(n) = k + ((lon-1) + (lat-1)*nx_global)*ndim3
            endif
         end do !ndim3
         enddo  !i
         enddo  !j
      end do    !iblk

      call pio_initdecomp(ice_pio_subsystem, pio_real, (/ndim3,nx_global,ny_global/), &
           dof3d, iodesc)

      deallocate(dof3d)

   end subroutine ice_pio_initdecomp_3d_inner

   subroutine ice_pio_initdecomp_4d (ndim3, ndim4, iodesc)

      integer(kind=int_kind), intent(in) :: ndim3, ndim4
      type(io_desc_t), intent(out) :: iodesc

      integer (kind=int_kind) :: &
          iblk,ilo,ihi,jlo,jhi,lon,lat,i,j,n,k,l 

      type(block) :: this_block 

      integer(kind=int_kind), pointer :: dof4d(:)

      allocate(dof4d(nx_block*ny_block*nblocks*ndim3*ndim4))

      n=0
      do l=1,ndim4
      do k=1,ndim3
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         
         do j=1,ny_block
         do i=1,nx_block
            n = n+1
            if (j < jlo .or. j>jhi) then
               dof4d(n)=0
            else if (i < ilo .or. i > ihi) then
               dof4d(n) = 0
            else
               lon = this_block%i_glob(i)
               lat = this_block%j_glob(j)
               dof4d(n) = ((lat-1)*nx_global + lon) &
                        + (k-1)*nx_global*ny_global & 
                        + (l-1)*nx_global*ny_global*ndim3 
            endif
         enddo !i
         enddo !j
      enddo ! iblk
      enddo !ndim3
      enddo !ndim4

      call pio_initdecomp(ice_pio_subsystem, pio_real, &
          (/nx_global,ny_global,ndim3,ndim4/), dof4d, iodesc)

      deallocate(dof4d)

   end subroutine ice_pio_initdecomp_4d
   
!================================================================================

  end module ice_pio

!================================================================================
