!  SVN:$Id: ice_history_write.F90 567 2013-01-07 02:57:36Z eclare $
!=======================================================================
!
! Writes history in netCDF format
!
! authors Tony Craig and Bruce Briegleb, NCAR
!         Elizabeth C. Hunke and William H. Lipscomb, LANL
!         C. M. Bitz, UW
!
! 2004 WHL: Block structure added
! 2006 ECH: Accepted some CCSM code into mainstream CICE
!           Added ice_present, aicen, vicen; removed aice1...10, vice1...1.
!           Added histfreq_n and histfreq='h' options, removed histfreq='w'
!           Converted to free source form (F90)
!           Added option for binary output instead of netCDF
! 2009 D Bailey and ECH: Generalized for multiple frequency output
! 2010 Alison McLaren and ECH: Added 3D capability
! 2013 ECH split from ice_history.F90

module ice_history_write

    use netcdf
    use mpi, only: MPI_INFO_NULL, MPI_COMM_WORLD

    use ice_kinds_mod
    use ice_constants, only: c0, c360, secday, spval, rad_to_deg
    use ice_blocks, only: nx_block, ny_block, block, get_block
    use ice_exit, only: abort_ice
    use ice_domain, only: distrb_info, nblocks, blocks_ice
    use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
    use ice_broadcast, only: broadcast_scalar
    use ice_gather_scatter, only: gather_global
    use ice_domain_size, only: nx_global, ny_global, max_nstrm, max_blocks
    use ice_grid, only: TLON, TLAT, ULON, ULAT, hm, bm, tarea, uarea, &
        dxu, dxt, dyu, dyt, HTN, HTE, ANGLE, ANGLET, &
        lont_bounds, latt_bounds, lonu_bounds, latu_bounds
    use ice_history_shared
    use ice_itd, only: hin_max
    use ice_calendar, only: write_ic, histfreq


    implicit none
    private
    public :: ice_write_hist
    save

    type coord_attributes         ! netcdf coordinate attributes
      character (len=11)   :: short_name
      character (len=45)   :: long_name
      character (len=20)   :: units
    end type coord_attributes

    type req_attributes         ! req'd netcdf attributes
      type (coord_attributes) :: req
      character (len=20)   :: coordinates
    end type req_attributes

    ! 4 coordinate variables: TLON, TLAT, ULON, ULAT
    INTEGER (kind=int_kind), PARAMETER :: ncoord = 4

    ! 4 vertices in each grid cell
    INTEGER (kind=int_kind), PARAMETER :: nverts = 4

    ! 4 variables describe T, U grid boundaries:
    ! lont_bounds, latt_bounds, lonu_bounds, latu_bounds
    INTEGER (kind=int_kind), PARAMETER :: nvar_verts = 4

    contains


subroutine check(status, msg)
    integer, intent (in) :: status
    character(len=*), intent (in) :: msg

    if(status /= nf90_noerr) then
        call abort_ice('ice: NetCDF error '//trim(nf90_strerror(status)//' '//trim(msg)))
    end if
end subroutine check


!=======================================================================
!
! write average ice quantities or snapshots
!
! author:   Elizabeth C. Hunke, LANL

subroutine ice_write_hist (ns)

#ifdef ncdf
    use ice_calendar, only: time, sec, idate, idate0, &
        dayyr, days_per_year, use_leap_years
    use ice_fileunits, only: nu_diag
    use ice_restart_shared, only: runid
#endif

    integer (kind=int_kind), intent(in) :: ns

    ! local variables

#ifdef ncdf
    real (kind=dbl_kind),  dimension(:,:),   allocatable :: work_g1
    real (kind=real_kind), dimension(:,:),   allocatable :: work_gr
    real (kind=real_kind), dimension(:,:,:), allocatable :: work_gr3
    real (kind=dbl_kind),  dimension(nx_block,ny_block,max_blocks) :: &
       work1

    integer (kind=int_kind) :: i,k,ic,n,nn, &
       ncid,status,imtid,jmtid,kmtidi,kmtids,kmtidb, cmtid,timid,varid, &
       nvertexid,ivertex
    integer (kind=int_kind), dimension(3) :: dimid
    integer (kind=int_kind), dimension(4) :: dimidz
    integer (kind=int_kind), dimension(5) :: dimidcz
    integer (kind=int_kind), dimension(3) :: dimid_nverts
    integer (kind=int_kind), dimension(4) :: dimidex
    real (kind=real_kind) :: ltime
    character (char_len) :: title
    character (char_len_long) :: ncfile(max_nstrm)

    integer (kind=int_kind) :: shuffle, deflate, deflate_level

    integer (kind=int_kind) :: ind,boundid

    character (char_len) :: start_time,current_date,current_time
    character (len=8) :: cdate

    TYPE(req_attributes), dimension(nvar) :: var
    TYPE(coord_attributes), dimension(ncoord) :: coord_var
    TYPE(coord_attributes), dimension(nvar_verts) :: var_nverts
    TYPE(coord_attributes), dimension(nvarz) :: var_nz
    CHARACTER (char_len), dimension(ncoord) :: coord_bounds

    logical :: do_parallel_io

    do_parallel_io = .true.

    ! We leave shuffle at 0, this is only useful for integer data.
    shuffle = 0

    ! If history_deflate_level < 0 then don't do deflation,
    ! otherwise it sets the deflate level
    if (history_deflate_level < 0) then
        deflate = 0
        deflate_level = 0
    else
        deflate = 1
        deflate_level = history_deflate_level
    endif

    if (my_task == master_task .or. do_parallel_io) then

        ltime=time/int(secday)

        call construct_filename(ncfile(ns),'nc',ns)

        ! add local directory path name to ncfile
        if (write_ic) then
            ncfile(ns) = trim(incond_dir)//ncfile(ns)
        else
            ncfile(ns) = trim(history_dir)//ncfile(ns)
        endif

        ! create file
        if (do_parallel_io) then
            call check(nf90_create(ncfile(ns), ior(NF90_NETCDF4, NF90_MPIIO), ncid, &
                                   comm=MPI_COMM_ICE, info=MPI_INFO_NULL), &
                        'create history ncfile '//ncfile(ns))
            ! FIXME: put in a check that each PE has the same number of blocks
        else
            call check(nf90_create(ncfile(ns), ior(NF90_CLASSIC_MODEL, NF90_HDF5), ncid), &
                        'create history ncfile '//ncfile(ns))
        endif

        !-----------------------------------------------------------------
        ! define dimensions
        !-----------------------------------------------------------------

        if (hist_avg) then
            call check(nf90_def_dim(ncid,'d2',2,boundid), 'def dim d2')
        endif

        call check(nf90_def_dim(ncid, 'ni', nx_global, imtid), &
                    'def dim ni')
        call check(nf90_def_dim(ncid, 'nj', ny_global, jmtid), &
                    'def dim nj')
        call check(nf90_def_dim(ncid, 'nc', ncat_hist, cmtid), &
                    'def dim nc')
        call check(nf90_def_dim(ncid, 'nkice', nzilyr, kmtidi), &
                    'def dim nkice')
        call check(nf90_def_dim(ncid, 'nksnow', nzslyr, kmtids), &
                    'def dim nksnow')
        call check(nf90_def_dim(ncid, 'nkbio', nzblyr, kmtidb), &
                    'def dim nkbio')
        call check(nf90_def_dim(ncid, 'time', 1, timid), &
                    'def dim time')
        call check(nf90_def_dim(ncid, 'nvertices', nverts, nvertexid), &
                    'def dim nverts')

        !-----------------------------------------------------------------
        ! define coordinate variables
        !-----------------------------------------------------------------

        call check(nf90_def_var(ncid,'time',nf90_float,timid,varid), &
                      'def var time')
        call check(nf90_put_att(ncid,varid,'long_name','model time'), &
                    'put_att long_name')

        write(cdate,'(i8.8)') idate0
        write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
              cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
        call check(nf90_put_att(ncid,varid,'units',title), &
                    'put_att time units')

        if (days_per_year == 360) then
            status = nf90_put_att(ncid,varid,'calendar','360_day')
            if (status /= nf90_noerr) call abort_ice( &
                         'ice Error: time calendar')
        elseif (days_per_year == 365 .and. .not.use_leap_years ) then
            status = nf90_put_att(ncid,varid,'calendar','NoLeap')
            if (status /= nf90_noerr) call abort_ice( &
                         'ice Error: time calendar')
        elseif (use_leap_years) then
            status = nf90_put_att(ncid,varid,'calendar','Gregorian')
            if (status /= nf90_noerr) call abort_ice( &
                         'ice Error: time calendar')
        else
            call abort_ice( 'ice Error: invalid calendar settings')
        endif

        if (hist_avg) then
            status = nf90_put_att(ncid,varid,'bounds','time_bounds')
            if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: time bounds')
        endif

        !-----------------------------------------------------------------
        ! Define attributes for time bounds if hist_avg is true
        !-----------------------------------------------------------------

        if (hist_avg) then
            dimid(1) = boundid
            dimid(2) = timid
            status = nf90_def_var(ncid, 'time_bounds', &
                        nf90_float,dimid(1:2),varid)
            if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error defining var time_bounds')

            status = nf90_put_att(ncid,varid,'long_name', &
                                'boundaries for time-averaging interval')
            if (status /= nf90_noerr) call abort_ice( &
                        'ice Error: time_bounds long_name')
            write(cdate,'(i8.8)') idate0
            write(title,'(a,a,a,a,a,a,a,a)') 'days since ', &
                    cdate(1:4),'-',cdate(5:6),'-',cdate(7:8),' 00:00:00'
            status = nf90_put_att(ncid,varid,'units',title)
            if (status /= nf90_noerr) call abort_ice( &
                        'ice Error: time_bounds units')
        endif

        !-----------------------------------------------------------------
        ! define information for required time-invariant variables
        !-----------------------------------------------------------------

        ind = 0
        ind = ind + 1
        coord_var(ind) = coord_attributes('TLON', &
                         'T grid center longitude', 'degrees_east')
        coord_bounds(ind) = 'lont_bounds'
        ind = ind + 1
        coord_var(ind) = coord_attributes('TLAT', &
                         'T grid center latitude',  'degrees_north')
        coord_bounds(ind) = 'latt_bounds'
        ind = ind + 1
        coord_var(ind) = coord_attributes('ULON', &
                         'U grid center longitude', 'degrees_east')
        coord_bounds(ind) = 'lonu_bounds'
        ind = ind + 1
        coord_var(ind) = coord_attributes('ULAT', &
                         'U grid center latitude',  'degrees_north')
        coord_bounds(ind) = 'latu_bounds'

        var_nz(1) = coord_attributes('NCAT', 'category maximum thickness', 'm')
        var_nz(2) = coord_attributes('VGRDi', 'vertical ice levels', '1')
        var_nz(3) = coord_attributes('VGRDs', 'vertical snow levels', '1')
        var_nz(4) = coord_attributes('VGRDb', 'vertical ice-bio levels', '1')

        !-----------------------------------------------------------------
        ! define information for optional time-invariant variables
        !-----------------------------------------------------------------

        var(n_tarea)%req = coord_attributes('tarea', &
                    'area of T grid cells', 'm^2')
        var(n_tarea)%coordinates = 'TLON TLAT'
        var(n_uarea)%req = coord_attributes('uarea', &
                    'area of U grid cells', 'm^2')
        var(n_uarea)%coordinates = 'ULON ULAT'
        var(n_dxt)%req = coord_attributes('dxt', &
                    'T cell width through middle', 'm')
        var(n_dxt)%coordinates = 'TLON TLAT'
        var(n_dyt)%req = coord_attributes('dyt', &
                    'T cell height through middle', 'm')
        var(n_dyt)%coordinates = 'TLON TLAT'
        var(n_dxu)%req = coord_attributes('dxu', &
                    'U cell width through middle', 'm')
        var(n_dxu)%coordinates = 'ULON ULAT'
        var(n_dyu)%req = coord_attributes('dyu', &
                    'U cell height through middle', 'm')
        var(n_dyu)%coordinates = 'ULON ULAT'
        var(n_HTN)%req = coord_attributes('HTN', &
                    'T cell width on North side','m')
        var(n_HTN)%coordinates = 'TLON TLAT'
        var(n_HTE)%req = coord_attributes('HTE', &
                    'T cell width on East side', 'm')
        var(n_HTE)%coordinates = 'TLON TLAT'
        var(n_ANGLE)%req = coord_attributes('ANGLE', &
                    'angle grid makes with latitude line on U grid', &
                    'radians')
        var(n_ANGLE)%coordinates = 'ULON ULAT'
        var(n_ANGLET)%req = coord_attributes('ANGLET', &
                    'angle grid makes with latitude line on T grid', &
                    'radians')
        var(n_ANGLET)%coordinates = 'TLON TLAT'

        ! These fields are required for CF compliance
        ! dimensions (nx,ny,nverts)
        var_nverts(n_lont_bnds) = coord_attributes('lont_bounds', &
                    'longitude boundaries of T cells', 'degrees_east')
        var_nverts(n_latt_bnds) = coord_attributes('latt_bounds', &
                    'latitude boundaries of T cells', 'degrees_north')
        var_nverts(n_lonu_bnds) = coord_attributes('lonu_bounds', &
                    'longitude boundaries of U cells', 'degrees_east')
        var_nverts(n_latu_bnds) = coord_attributes('latu_bounds', &
                    'latitude boundaries of U cells', 'degrees_north')

        !-----------------------------------------------------------------
        ! define attributes for time-invariant variables
        !-----------------------------------------------------------------

        dimid(1) = imtid
        dimid(2) = jmtid
        dimid(3) = timid

        do i = 1, ncoord
          status = nf90_def_var(ncid, coord_var(i)%short_name, nf90_float, &
                                dimid(1:2), varid)
          if (status /= nf90_noerr) call abort_ice( &
               'Error defining short_name for '//coord_var(i)%short_name)

          status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                        deflate_level)
          if (status /= nf90_noerr) call abort_ice( &
               'Error deflate short_name for '//coord_var(i)%short_name)

          status = nf90_put_att(ncid,varid,'long_name',coord_var(i)%long_name)
          if (status /= nf90_noerr) call abort_ice( &
               'Error defining long_name for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid, varid, 'units', coord_var(i)%units)
          if (status /= nf90_noerr) call abort_ice( &
                  'Error defining units for '//coord_var(i)%short_name)
          status = nf90_put_att(ncid,varid,'missing_value',spval)
          if (status /= nf90_noerr) call abort_ice( &
             'Error defining missing_value for '//coord_var(i)%short_name)

          call check(nf90_put_att(ncid, varid, '_FillValue', spval), &
                        'put att _FillValue for '//coord_var(i)%short_name)

          if (coord_var(i)%short_name == 'ULAT') then
             status = nf90_put_att(ncid,varid,'comment', &
                  'Latitude of NE corner of T grid cell')
             if (status /= nf90_noerr) call abort_ice( &
                  'Error defining comment for '//coord_var(i)%short_name)
          endif
          if (f_bounds) then
              status = nf90_put_att(ncid, varid, 'bounds', coord_bounds(i))
              if (status /= nf90_noerr) call abort_ice( &
                  'Error defining bounds for '//coord_var(i)%short_name)
          endif

        enddo


        ! Extra dimensions (NCAT, NZILYR, NZSLYR, NZBLYR)
        dimidex(1)=cmtid
        dimidex(2)=kmtidi
        dimidex(3)=kmtids
        dimidex(4)=kmtidb

        do i = 1, nvarz
            if (igrdz(i)) then
                status = nf90_def_var(ncid, var_nz(i)%short_name, &
                                   nf90_float, dimidex(i), varid)
            if (status /= nf90_noerr) call abort_ice( &
                'Error defining short_name for '//var_nz(i)%short_name)

            status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                            deflate_level)
            if (status /= nf90_noerr) call abort_ice( &
                'Error defining short_name for '//var_nz(i)%short_name)

           status = nf90_put_att(ncid,varid,'long_name',var_nz(i)%long_name)
           if (status /= nf90_noerr) call abort_ice( &
                'Error defining long_name for '//var_nz(i)%short_name)
           status = nf90_put_att(ncid, varid, 'units', var_nz(i)%units)
           if (Status /= nf90_noerr) call abort_ice( &
                'Error defining units for '//var_nz(i)%short_name)
            endif
        enddo

        ! Attributes for tmask, blkmask defined separately, since they have no units
        if (igrd(n_tmask)) then
            status = nf90_def_var(ncid, 'tmask', nf90_float, dimid(1:2), varid)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error defining var tmask')

            status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                          deflate_level)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error deflating var tmask')

            status = nf90_put_att(ncid,varid, 'long_name', 'ocean grid mask')
            if (status /= nf90_noerr) call abort_ice('ice Error: tmask long_name')
            status = nf90_put_att(ncid, varid, 'coordinates', 'TLON TLAT')
            if (status /= nf90_noerr) call abort_ice('ice Error: tmask units')
            status = nf90_put_att(ncid,varid,'comment', '0 = land, 1 = ocean')
            if (status /= nf90_noerr) call abort_ice('ice Error: tmask comment')
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice('Error defining missing_value for tmask')
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice('Error defining _FillValue for tmask')
        endif

        if (igrd(n_blkmask)) then
            status = nf90_def_var(ncid, 'blkmask', nf90_float, dimid(1:2), varid)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error defining var blkmask')

            status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                          deflate_level)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error deflating var blkmask')

            status = nf90_put_att(ncid,varid, 'long_name', 'ice grid block mask')
            if (status /= nf90_noerr) call abort_ice('ice Error: blkmask long_name')
            status = nf90_put_att(ncid, varid, 'coordinates', 'TLON TLAT')
            if (status /= nf90_noerr) call abort_ice('ice Error: blkmask units')
            status = nf90_put_att(ncid,varid,'comment', 'mytask + iblk/100')
            if (status /= nf90_noerr) call abort_ice('ice Error: blkmask comment')
            status = nf90_put_att(ncid,varid,'missing_value',spval)
            if (status /= nf90_noerr) call abort_ice('Error defining missing_value for blkmask')
            status = nf90_put_att(ncid,varid,'_FillValue',spval)
            if (status /= nf90_noerr) call abort_ice('Error defining _FillValue for blkmask')
        endif

        do i = 3, nvar      ! note n_tmask=1, n_blkmask=2
            if (igrd(i)) then
                status = nf90_def_var(ncid, var(i)%req%short_name, &
                                      nf90_float, dimid(1:2), varid)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error defining variable '//var(i)%req%short_name)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error deflating variable '//var(i)%req%short_name)

                status = nf90_put_att(ncid,varid, 'long_name', var(i)%req%long_name)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error defining long_name for '//var(i)%req%short_name)
                status = nf90_put_att(ncid, varid, 'units', var(i)%req%units)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error defining units for '//var(i)%req%short_name)
                status = nf90_put_att(ncid, varid, 'coordinates', var(i)%coordinates)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error defining coordinates for '//var(i)%req%short_name)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining missing_value for '//var(i)%req%short_name)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining _FillValue for '//var(i)%req%short_name)
            endif
        enddo

        ! Fields with dimensions (nverts,nx,ny)
        dimid_nverts(1) = nvertexid
        dimid_nverts(2) = imtid
        dimid_nverts(3) = jmtid
        do i = 1, nvar_verts
            if (f_bounds) then
                status = nf90_def_var(ncid, var_nverts(i)%short_name, &
                                      nf90_float,dimid_nverts, varid)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error defining variable '//var_nverts(i)%short_name)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error deflating variable '//var_nverts(i)%short_name)

                status = nf90_put_att(ncid,varid, 'long_name', &
                            var_nverts(i)%long_name)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error defining long_name for '//var_nverts(i)%short_name)
                status = nf90_put_att(ncid, varid, 'units', var_nverts(i)%units)
                if (status /= nf90_noerr) call abort_ice( &
                     'Error defining units for '//var_nverts(i)%short_name)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining missing_value for '//var_nverts(i)%short_name)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining _FillValue for '//var_nverts(i)%short_name)
            endif
        enddo

        do n=1,num_avail_hist_fields_2D
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
                status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                             nf90_float, dimid, varid)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining variable '//avail_hist_fields(n)%vname)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error deflating variable '//avail_hist_fields(n)%vname)

                status = nf90_put_att(ncid,varid,'units', &
                            avail_hist_fields(n)%vunit)
                if (status /= nf90_noerr) call abort_ice( &
                    'Error defining units for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid, 'long_name', &
                            avail_hist_fields(n)%vdesc)
                if (status /= nf90_noerr) call abort_ice( &
                    'Error defining long_name for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'coordinates', &
                            avail_hist_fields(n)%vcoord)
                if (status /= nf90_noerr) call abort_ice( &
                    'Error defining coordinates for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'cell_measures', &
                            avail_hist_fields(n)%vcellmeas)
                if (status /= nf90_noerr) call abort_ice( &
                    'Error defining cell measures for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                    'Error defining missing_value for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                    'Error defining _FillValue for '//avail_hist_fields(n)%vname)

                !---------------------------------------------------------------
                ! Add cell_methods attribute to variables if averaged
                !---------------------------------------------------------------
                if (hist_avg) then
                    if (TRIM(avail_hist_fields(n)%vname)/='sig1' .or. &
                        TRIM(avail_hist_fields(n)%vname)/='sig2') then

                        status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                        if (status /= nf90_noerr) call abort_ice( &
                            'Error defining cell methods for '//avail_hist_fields(n)%vname)
                    endif
                endif

                if (histfreq(ns) == '1' .or. .not. hist_avg         &
                    .or. n==n_divu(ns)      .or. n==n_shear(ns)     &  ! snapshots
                    .or. n==n_sig1(ns)      .or. n==n_sig2(ns)      &
                    .or. n==n_trsig(ns)                             &
                    .or. n==n_mlt_onset(ns) .or. n==n_frz_onset(ns) &
                    .or. n==n_hisnap(ns)    .or. n==n_aisnap(ns)) then
                    status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
                else
                    status = nf90_put_att(ncid,varid,'time_rep','averaged')
                endif
            endif
        enddo  ! num_avail_hist_fields_2D

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = cmtid
        dimidz(4) = timid

        do n = n2D + 1, n3Dccum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
                status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                             nf90_float, dimidz, varid)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining variable '//avail_hist_fields(n)%vname)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error deflating variable '//avail_hist_fields(n)%vname)

                status = nf90_put_att(ncid,varid,'units', &
                            avail_hist_fields(n)%vunit)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining units for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid, 'long_name', &
                            avail_hist_fields(n)%vdesc)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining long_name for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'coordinates', &
                            avail_hist_fields(n)%vcoord)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining coordinates for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'cell_measures', &
                            avail_hist_fields(n)%vcellmeas)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining cell measures for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining missing_value for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining _FillValue for '//avail_hist_fields(n)%vname)

                !---------------------------------------------------------------
                ! Add cell_methods attribute to variables if averaged
                !---------------------------------------------------------------
                if (hist_avg) then
                    status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                    if (status /= nf90_noerr) call abort_ice( &
                     'Error defining cell methods for '//avail_hist_fields(n)%vname)
                endif

                if (histfreq(ns) == '1' .or. .not. hist_avg) then
                    status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
                else
                    status = nf90_put_att(ncid,varid,'time_rep','averaged')
                endif
            endif
        enddo  ! num_avail_hist_fields_3Dc

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = kmtidi
        dimidz(4) = timid

        do n = n3Dccum + 1, n3Dzcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
                status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                             nf90_float, dimidz, varid)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining variable '//avail_hist_fields(n)%vname)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error deflating variable '//avail_hist_fields(n)%vname)

                status = nf90_put_att(ncid,varid,'units', &
                            avail_hist_fields(n)%vunit)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining units for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid, 'long_name', &
                            avail_hist_fields(n)%vdesc)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining long_name for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'coordinates', &
                            avail_hist_fields(n)%vcoord)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining coordinates for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'cell_measures', &
                            avail_hist_fields(n)%vcellmeas)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining cell measures for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining missing_value for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining _FillValue for '//avail_hist_fields(n)%vname)

            endif
        enddo  ! num_avail_hist_fields_3Dz

        dimidz(1) = imtid
        dimidz(2) = jmtid
        dimidz(3) = kmtidb
        dimidz(4) = timid

        do n = n3Dzcum + 1, n3Dbcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
                status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                             nf90_float, dimidz, varid)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining variable '//avail_hist_fields(n)%vname)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error deflating variable '//avail_hist_fields(n)%vname)

                status = nf90_put_att(ncid,varid,'units', &
                            avail_hist_fields(n)%vunit)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining units for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid, 'long_name', &
                            avail_hist_fields(n)%vdesc)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining long_name for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'coordinates', &
                            avail_hist_fields(n)%vcoord)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining coordinates for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'cell_measures', &
                            avail_hist_fields(n)%vcellmeas)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining cell measures for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining missing_value for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining _FillValue for '//avail_hist_fields(n)%vname)
            endif
        enddo  ! num_avail_hist_fields_3Db

        dimidcz(1) = imtid
        dimidcz(2) = jmtid
        dimidcz(3) = kmtidi
        dimidcz(4) = cmtid
        dimidcz(5) = timid

        do n = n3Dbcum + 1, n4Dicum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
                status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                                 !nf90_float, dimidcz, varid)
                                 nf90_float, dimidcz(1:4), varid) ! ferret
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining variable '//avail_hist_fields(n)%vname)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error deflating variable '//avail_hist_fields(n)%vname)

                status = nf90_put_att(ncid,varid,'units', &
                            avail_hist_fields(n)%vunit)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining units for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid, 'long_name', &
                            avail_hist_fields(n)%vdesc)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining long_name for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'coordinates', &
                            avail_hist_fields(n)%vcoord)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining coordinates for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'cell_measures', &
                            avail_hist_fields(n)%vcellmeas)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining cell measures for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining missing_value for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining _FillValue for '//avail_hist_fields(n)%vname)

                !---------------------------------------------------------------
                ! Add cell_methods attribute to variables if averaged
                !---------------------------------------------------------------
                if (hist_avg) then
                    status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                    if (status /= nf90_noerr) call abort_ice( &
                        'Error defining cell methods for '//avail_hist_fields(n)%vname)
                endif

                if (histfreq(ns) == '1' .or. .not. hist_avg) then
                    status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
                else
                    status = nf90_put_att(ncid,varid,'time_rep','averaged')
                endif
            endif
        enddo  ! num_avail_hist_fields_4Di

        dimidcz(1) = imtid
        dimidcz(2) = jmtid
        dimidcz(3) = kmtids
        dimidcz(4) = cmtid
        dimidcz(5) = timid

        do n = n4Dicum + 1, n4Dscum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
                status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                                 !nf90_float, dimidcz, varid)
                                 nf90_float, dimidcz(1:4), varid) ! ferret
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining variable '//avail_hist_fields(n)%vname)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error deflating variable '//avail_hist_fields(n)%vname)

                status = nf90_put_att(ncid,varid,'units', &
                            avail_hist_fields(n)%vunit)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining units for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid, 'long_name', &
                            avail_hist_fields(n)%vdesc)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining long_name for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'coordinates', &
                            avail_hist_fields(n)%vcoord)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining coordinates for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'cell_measures', &
                            avail_hist_fields(n)%vcellmeas)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining cell measures for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining missing_value for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining _FillValue for '//avail_hist_fields(n)%vname)

                !---------------------------------------------------------------
                ! Add cell_methods attribute to variables if averaged
                !---------------------------------------------------------------
                if (hist_avg) then
                    status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                    if (status /= nf90_noerr) call abort_ice( &
                     'Error defining cell methods for '//avail_hist_fields(n)%vname)
                endif

                if (histfreq(ns) == '1' .or. .not. hist_avg) then
                   status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
                else
                   status = nf90_put_att(ncid,varid,'time_rep','averaged')
                endif
            endif
        enddo  ! num_avail_hist_fields_4Ds

        dimidcz(1) = imtid
        dimidcz(2) = jmtid
        dimidcz(3) = kmtidb
        dimidcz(4) = cmtid
        dimidcz(5) = timid

        do n = n4Dscum + 1, n4Dbcum
            if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
                status  = nf90_def_var(ncid, avail_hist_fields(n)%vname, &
                                 !nf90_float, dimidcz, varid)
                                 nf90_float, dimidcz(1:4), varid) ! ferret
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining variable '//avail_hist_fields(n)%vname)

                status = nf90_def_var_deflate(ncid, varid, shuffle, deflate, &
                                              deflate_level)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error deflating variable '//avail_hist_fields(n)%vname)

                status = nf90_put_att(ncid,varid,'units', &
                            avail_hist_fields(n)%vunit)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining units for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid, 'long_name', &
                            avail_hist_fields(n)%vdesc)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining long_name for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'coordinates', &
                            avail_hist_fields(n)%vcoord)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining coordinates for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'cell_measures', &
                            avail_hist_fields(n)%vcellmeas)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining cell measures for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'missing_value',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining missing_value for '//avail_hist_fields(n)%vname)
                status = nf90_put_att(ncid,varid,'_FillValue',spval)
                if (status /= nf90_noerr) call abort_ice( &
                   'Error defining _FillValue for '//avail_hist_fields(n)%vname)

                !---------------------------------------------------------------
                ! Add cell_methods attribute to variables if averaged
                !---------------------------------------------------------------
                if (hist_avg) then
                    status = nf90_put_att(ncid,varid,'cell_methods','time: mean')
                    if (status /= nf90_noerr) call abort_ice( &
                     'Error defining cell methods for '//avail_hist_fields(n)%vname)
                endif

                if (histfreq(ns) == '1' .or. .not. hist_avg) then
                   status = nf90_put_att(ncid,varid,'time_rep','instantaneous')
                else
                   status = nf90_put_att(ncid,varid,'time_rep','averaged')
                endif
            endif
        enddo  ! num_avail_hist_fields_4Db

        !-----------------------------------------------------------------
        ! global attributes
        !-----------------------------------------------------------------
        ! ... the user should change these to something useful ...
        !-----------------------------------------------------------------
#ifdef CCSMCOUPLED
        status = nf90_put_att(ncid,nf90_global,'title',runid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error in global attribute title')
#else
        title  = 'sea ice model output for CICE'
        status = nf90_put_att(ncid,nf90_global,'title',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error in global attribute title')
#endif
        title = 'Diagnostic and Prognostic Variables'
        status = nf90_put_att(ncid,nf90_global,'contents',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute contents')

        title  = 'Los Alamos Sea Ice Model (CICE) Version 5'
        status = nf90_put_att(ncid,nf90_global,'source',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute source')

#ifdef AusCOM
        write(title,'(a,i3,a)') 'This Year Has ',int(dayyr),' days'
#else
        if (use_leap_years) then
            write(title,'(a,i3,a)') 'This year has ',int(dayyr),' days'
        else
            write(title,'(a,i3,a)') 'All years have exactly ',int(dayyr),' days'
        endif
#endif
        status = nf90_put_att(ncid,nf90_global,'comment',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute comment')

        write(title,'(a,i8.8)') 'File written on model date ',idate
        status = nf90_put_att(ncid,nf90_global,'comment2',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute date1')

        write(title,'(a,i6)') 'seconds elapsed into model date: ',sec
        status = nf90_put_att(ncid,nf90_global,'comment3',title)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute date2')

        title = 'CF-1.0'
        status =  &
             nf90_put_att(ncid,nf90_global,'conventions',title)
        if (status /= nf90_noerr) call abort_ice( &
             'Error in global attribute conventions')

        call date_and_time(date=current_date, time=current_time)
        write(start_time,1000) current_date(1:4), current_date(5:6), &
                               current_date(7:8), current_time(1:2), &
                               current_time(3:4), current_time(5:8)
1000    format('This dataset was created on ', &
                a,'-',a,'-',a,' at ',a,':',a,':',a)

        status = nf90_put_att(ncid,nf90_global,'history',start_time)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute history')

        status = nf90_put_att(ncid,nf90_global,'io_flavor','io_netcdf')
        if (status /= nf90_noerr) call abort_ice( &
                      'ice Error: global attribute io_flavor')

        !-----------------------------------------------------------------
        ! end define mode
        !-----------------------------------------------------------------

        call check(nf90_enddef(ncid), 'enddef')

        !-----------------------------------------------------------------
        ! write time variable
        !-----------------------------------------------------------------

        status = nf90_inq_varid(ncid,'time',varid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error getting time varid')
        status = nf90_put_var(ncid,varid,ltime)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error writing time variable')

        !-----------------------------------------------------------------
        ! write time_bounds info
        !-----------------------------------------------------------------

        if (hist_avg) then
            status = nf90_inq_varid(ncid,'time_bounds',varid)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error getting time_bounds id')
            status = nf90_put_var(ncid,varid,time_beg(ns),start=(/1/))
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error writing time_beg')
            status = nf90_put_var(ncid,varid,time_end(ns),start=(/2/))
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error writing time_end')
        endif
    endif                     ! master_task or do_parallel_io

    !-----------------------------------------------------------------
    ! write coordinate variables
    !-----------------------------------------------------------------

    if (do_parallel_io) then
        call write_coordinate_variables_parallel(ncid, coord_var, var_nz)
    else
        call write_coordinate_variables(ncid, coord_var, var_nz)
    endif

    !-----------------------------------------------------------------
    ! write grid masks, area and rotation angle
    !-----------------------------------------------------------------

    if (do_parallel_io) then
        call write_grid_variables_parallel(ncid, var, var_nverts)
    else
        call write_grid_variables(ncid, var, var_nverts)
    endif


    !-----------------------------------------------------------------
    ! write 2d variable data
    !-----------------------------------------------------------------

    if (do_parallel_io) then
        call write_2d_variables_parallel(ns, ncid)
    else
        call write_2d_variables(ns, ncid)
    endif

    if (do_parallel_io) then
        call write_3d_and_4d_variables_parallel(ns, ncid)
    else
        call write_3d_and_4d_variables(ns, ncid)
    endif

    !-----------------------------------------------------------------
    ! close output dataset
    !-----------------------------------------------------------------

    if (my_task == master_task .or. do_parallel_io) then
        status = nf90_close(ncid)
        if (status /= nf90_noerr) call abort_ice( &
                      'ice: Error closing netCDF history file')
        write(nu_diag,*) ' '
        write(nu_diag,*) 'Finished writing ',trim(ncfile(ns))
    endif
#endif

end subroutine ice_write_hist

subroutine write_coordinate_variables(ncid, coord_var, var_nz)

    integer, intent(in) :: ncid
    type(coord_attributes), dimension(ncoord), intent(in) :: coord_var
    type(coord_attributes), dimension(nvarz) :: var_nz

    real (kind=dbl_kind),  dimension(:,:),   allocatable :: work_g1
    real (kind=real_kind), dimension(:,:),   allocatable :: work_gr
    real (kind=dbl_kind),  dimension(nx_block,ny_block,max_blocks) :: work1

    integer :: i, k, status
    integer :: varid
    character (len=len(coord_var(1)%short_name)) :: coord_var_name

    if (my_task==master_task) then
        allocate(work_g1(nx_global,ny_global))
        allocate(work_gr(nx_global,ny_global))
    else
        allocate(work_g1(1,1))
        allocate(work_gr(1,1))   ! to save memory
    endif

    work_g1(:,:) = c0

    do i = 1,ncoord
        coord_var_name = coord_var(i)%short_name

        call broadcast_scalar(coord_var_name, master_task)
        SELECT CASE (coord_var_name)
        CASE ('TLON')
            ! Convert T grid longitude from -180 -> 180 to 0 to 360
            work1 = TLON*rad_to_deg + c360
            where (work1 > c360) work1 = work1 - c360
            where (work1 < c0 )  work1 = work1 + c360
            call gather_global(work_g1,work1,master_task,distrb_info)
        CASE ('TLAT')
            work1 = TLAT*rad_to_deg
            call gather_global(work_g1,work1,master_task,distrb_info)
        CASE ('ULON')
            work1 = ULON*rad_to_deg
            call gather_global(work_g1,work1,master_task,distrb_info)
        CASE ('ULAT')
            work1 = ULAT*rad_to_deg
            call gather_global(work_g1,work1,master_task,distrb_info)
        END SELECT

        if (my_task == master_task) then
            work_gr = work_g1
            status = nf90_inq_varid(ncid, coord_var_name, varid)
            if (status /= nf90_noerr) call abort_ice( &
                 'ice: Error getting varid for '//coord_var_name)
            status = nf90_put_var(ncid,varid,work_gr)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error writing'//coord_var_name)
        endif
    enddo

    ! Extra dimensions (NCAT, VGRD*)

    do i = 1, nvarz
        if (igrdz(i)) then
            call broadcast_scalar(var_nz(i)%short_name,master_task)
            if (my_task == master_task) then
                status = nf90_inq_varid(ncid, var_nz(i)%short_name, varid)
                if (status /= nf90_noerr) call abort_ice( &
                     'ice: Error getting varid for '//var_nz(i)%short_name)
                SELECT CASE (var_nz(i)%short_name)
                CASE ('NCAT')
                    status = nf90_put_var(ncid,varid,hin_max(1:ncat_hist))
                CASE ('VGRDi') ! index - needed for Met Office analysis code
                    status = nf90_put_var(ncid,varid,(/(k, k=1,nzilyr)/))
                CASE ('VGRDs') ! index - needed for Met Office analysis code
                    status = nf90_put_var(ncid,varid,(/(k, k=1,nzslyr)/))
                CASE ('VGRDb')
                    status = nf90_put_var(ncid,varid,(/(k, k=1,nzblyr)/))
                END SELECT
                if (status /= nf90_noerr) call abort_ice( &
                              'ice: Error writing'//var_nz(i)%short_name)
            endif
        endif
    enddo

    deallocate(work_g1)
    deallocate(work_gr)

end subroutine write_coordinate_variables



subroutine write_grid_variables(ncid, var, var_nverts)

    integer, intent(in) :: ncid
    type(req_attributes), dimension(nvar), intent(in) :: var
    type(coord_attributes), dimension(nvar_verts), intent(in) :: var_nverts

    real (kind=dbl_kind),  dimension(:,:), allocatable :: work_g1
    real (kind=real_kind), dimension(:,:), allocatable :: work_gr
    real (kind=real_kind), dimension(:,:, :), allocatable :: work_gr3
    real (kind=dbl_kind),  dimension(nx_block,ny_block,max_blocks) :: work1

    integer :: ivertex, i, status
    integer :: varid
    character (len=len(var(1)%req%short_name)) :: var_name
    character (len=len(var_nverts(1)%short_name)) :: var_nverts_name

    if (my_task == master_task) then
        allocate(work_g1(nx_global,ny_global))
        allocate(work_gr(nx_global,ny_global))
        allocate(work_gr3(nverts,nx_global,ny_global))
    else
        allocate(work_g1(1,1))
        allocate(work_gr(1,1))   ! to save memory
        allocate(work_gr3(1,1,1))
    endif

    work_g1(:,:) = c0
    work_gr(:,:) = c0
    work_gr3(:,:,:) = c0

    if (igrd(n_tmask)) then
        call gather_global(work_g1, hm, master_task, distrb_info)
        if (my_task == master_task) then
            work_gr = work_g1
            status = nf90_inq_varid(ncid, 'tmask', varid)
            if (status /= nf90_noerr) call abort_ice( &
                                'ice: Error getting varid for tmask')
            status = nf90_put_var(ncid,varid,work_gr)
            if (status /= nf90_noerr) call abort_ice( &
                        'ice: Error writing variable tmask')
        endif
    endif

    if (igrd(n_blkmask)) then
        call gather_global(work_g1, bm, master_task, distrb_info)
        if (my_task == master_task) then
            work_gr=work_g1
            status = nf90_inq_varid(ncid, 'blkmask', varid)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error getting varid for blkmask')
            status = nf90_put_var(ncid,varid,work_gr)
            if (status /= nf90_noerr) call abort_ice( &
                          'ice: Error writing variable blkmask')
        endif
    endif

    do i = 3, nvar      ! note n_tmask=1, n_blkmask=2
        if (igrd(i)) then
            var_name = var(i)%req%short_name

            call broadcast_scalar(var_name,master_task)
            SELECT CASE (var_name)
            CASE ('tarea')
                call gather_global(work_g1, tarea, master_task, distrb_info)
            CASE ('uarea')
                call gather_global(work_g1, uarea, master_task, distrb_info)
            CASE ('dxu')
                call gather_global(work_g1,   dxu, master_task, distrb_info)
            CASE ('dyu')
              call gather_global(work_g1,   dyu, master_task, distrb_info)
            CASE ('dxt')
              call gather_global(work_g1,   dxt, master_task, distrb_info)
            CASE ('dyt')
              call gather_global(work_g1,   dyt, master_task, distrb_info)
            CASE ('HTN')
              call gather_global(work_g1,   HTN, master_task, distrb_info)
            CASE ('HTE')
              call gather_global(work_g1,   HTE, master_task, distrb_info)
            CASE ('ANGLE')
              call gather_global(work_g1, ANGLE, master_task, distrb_info)
            CASE ('ANGLET')
              call gather_global(work_g1, ANGLET,master_task, distrb_info)
            END SELECT

            if (my_task == master_task) then
              work_gr=work_g1
              status = nf90_inq_varid(ncid, var_name, varid)
              if (status /= nf90_noerr) call abort_ice( &
                            'ice: Error getting varid for '//var_name)
              status = nf90_put_var(ncid,varid,work_gr)
              if (status /= nf90_noerr) call abort_ice( &
                            'ice: Error writing variable '//var_name)
            endif
        endif
    enddo

    !----------------------------------------------------------------
    ! Write coordinates of grid box vertices
    !----------------------------------------------------------------

    if (f_bounds) then
        work_gr3(:,:,:) = c0
        work1   (:,:,:) = c0

        do i = 1, nvar_verts
            var_nverts_name = var_nverts(i)%short_name
            call broadcast_scalar(var_nverts_name,master_task)
            SELECT CASE (var_nverts_name)
            CASE ('lont_bounds')
                do ivertex = 1, nverts
                    work1(:,:,:) = lont_bounds(ivertex,:,:,:)
                    call gather_global(work_g1, work1, master_task, distrb_info)
                    if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
                enddo
            CASE ('latt_bounds')
                do ivertex = 1, nverts
                    work1(:,:,:) = latt_bounds(ivertex,:,:,:)
                    call gather_global(work_g1, work1, master_task, distrb_info)
                    if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
                enddo
            CASE ('lonu_bounds')
                do ivertex = 1, nverts
                    work1(:,:,:) = lonu_bounds(ivertex,:,:,:)
                    call gather_global(work_g1, work1, master_task, distrb_info)
                    if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
                enddo
            CASE ('latu_bounds')
                do ivertex = 1, nverts
                    work1(:,:,:) = latu_bounds(ivertex,:,:,:)
                    call gather_global(work_g1, work1, master_task, distrb_info)
                    if (my_task == master_task) work_gr3(ivertex,:,:) = work_g1(:,:)
                enddo
            END SELECT

            if (my_task == master_task) then
                status = nf90_inq_varid(ncid, var_nverts_name, varid)
                if (status /= nf90_noerr) call abort_ice( &
                   'ice: Error getting varid for '//var_nverts_name)
                status = nf90_put_var(ncid,varid,work_gr3)
                if (status /= nf90_noerr) call abort_ice( &
                   'ice: Error writing variable '//var_nverts_name)
            endif
        enddo
    endif

    deallocate(work_g1)
    deallocate(work_gr)
    deallocate(work_gr3)

end subroutine write_grid_variables


subroutine write_2d_variables(ns, ncid)

    integer, intent(in) :: ns
    integer, intent(in) :: ncid

    real (kind=dbl_kind),  dimension(:,:), allocatable :: work_g1
    real (kind=real_kind), dimension(:,:), allocatable :: work_gr

    integer :: n, status
    integer :: varid

    if (my_task == master_task) then
       allocate(work_g1(nx_global,ny_global))
       allocate(work_gr(nx_global,ny_global))
    else
       allocate(work_g1(1,1))
       allocate(work_gr(1,1))     ! to save memory
    endif

    work_gr(:,:) = c0
    work_g1(:,:) = c0

    do n=1, num_avail_hist_fields_2D
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call gather_global(work_g1, a2D(:,:,n,:), &
                               master_task, distrb_info)
            if (my_task == master_task) then
                work_gr(:,:) = work_g1(:,:)
                call check(nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid), &
                           'inq varid '//avail_hist_fields(n)%vname)
                call check(nf90_put_var(ncid,varid,work_gr(:,:), &
                                        count=(/nx_global,ny_global/)), &
                            'put var '//avail_hist_fields(n)%vname)
            endif
        endif
    enddo ! num_avail_hist_fields_2D

    deallocate(work_g1)
    deallocate(work_gr)

end subroutine write_2d_variables


subroutine write_3d_and_4d_variables(ns, ncid)

    integer, intent(in) :: ns
    integer, intent(in) :: ncid

    real (kind=dbl_kind),  dimension(:,:), allocatable :: work_g1
    real (kind=real_kind), dimension(:,:), allocatable :: work_gr

    integer :: varid
    integer :: status, n, nn, k, ic
    
    if (my_task == master_task) then
       allocate(work_g1(nx_global,ny_global))
       allocate(work_gr(nx_global,ny_global))
    else
       allocate(work_g1(1,1))
       allocate(work_gr(1,1))     ! to save memory
    endif

    work_gr(:,:) = c0
    work_g1(:,:) = c0

    do n = n2D + 1, n3Dccum
        nn = n - n2D
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
                call check(nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid), &
                           'inq varid '//avail_hist_fields(n)%vname)
            endif
            do k = 1, ncat_hist
                call gather_global(work_g1, a3Dc(:,:,k,nn,:), &
                                   master_task, distrb_info)
                work_gr(:,:) = work_g1(:,:)

                if (my_task == master_task) then
                    call check(nf90_put_var(ncid,varid,work_gr(:,:), &
                                            start=(/        1,        1, k/), &
                                            count=(/nx_global,ny_global, 1/)), &
                               'put var '//avail_hist_fields(n)%vname)
                endif
            enddo ! k
        endif
    enddo ! num_avail_hist_fields_3Dc

    work_gr(:,:) = c0
    work_g1(:,:) = c0

    do n = n3Dccum+1, n3Dzcum
        nn = n - n3Dccum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
                call check(nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid), &
                           'inq varid '//avail_hist_fields(n)%vname)
            endif
            do k = 1, nzilyr
                call gather_global(work_g1, a3Dz(:,:,k,nn,:), &
                                  master_task, distrb_info)
                work_gr(:,:) = work_g1(:,:)

                if (my_task == master_task) then
                    call check(nf90_put_var(ncid,varid,work_gr(:,:), &
                                            start=(/        1,        1,k/), &
                                            count=(/nx_global,ny_global,1/)), &
                               'put var '//avail_hist_fields(n)%vname)
                endif
            enddo ! k
        endif
    enddo ! num_avail_hist_fields_3Dz

    work_gr(:,:) = c0
    work_g1(:,:) = c0

    do n = n3Dzcum+1, n3Dbcum
        nn = n - n3Dzcum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
                call check(nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid), &
                           'inq varid '//avail_hist_fields(n)%vname)
            endif
            do k = 1, nzblyr
                call gather_global(work_g1, a3Db(:,:,k,nn,:), &
                                  master_task, distrb_info)
                work_gr(:,:) = work_g1(:,:)

                if (my_task == master_task) then
                    call check(nf90_put_var(ncid,varid,work_gr(:,:),    &
                                            start=(/        1,        1,k/), &
                                            count=(/nx_global,ny_global,1/)), &
                               'put var '//avail_hist_fields(n)%vname)
                endif
             enddo ! k
        endif
    enddo ! num_avail_hist_fields_3Db

    work_gr(:,:) = c0
    work_g1(:,:) = c0

    do n = n3Dbcum+1, n4Dicum
        nn = n - n3Dbcum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
                call check(nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid), &
                           'inq varid '//avail_hist_fields(n)%vname)
            endif
            do ic = 1, ncat_hist
                do k = 1, nzilyr
                    call gather_global(work_g1, a4Di(:,:,k,ic,nn,:), &
                                    master_task, distrb_info)
                    work_gr(:,:) = work_g1(:,:)
                    if (my_task == master_task) then
                        call check(nf90_put_var(ncid,varid,work_gr(:,:), &
                                                start=(/        1,        1,k,ic/), &
                                                count=(/nx_global,ny_global,1, 1/)), &
                                   'put var '//avail_hist_fields(n)%vname)
                    endif
                enddo ! k
            enddo ! ic
        endif
    enddo ! num_avail_hist_fields_4Di

    work_gr(:,:) = c0
    work_g1(:,:) = c0

    do n = n4Dicum+1, n4Dscum
        nn = n - n4Dicum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
                call check(nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid), &
                           'inq var '//avail_hist_fields(n)%vname)
            endif
            do ic = 1, ncat_hist
                do k = 1, nzslyr
                    call gather_global(work_g1, a4Ds(:,:,k,ic,nn,:), &
                                    master_task, distrb_info)
                    work_gr(:,:) = work_g1(:,:)
                    if (my_task == master_task) then
                        call check(nf90_put_var(ncid,varid,work_gr(:,:), &
                                                start=(/        1,        1,k,ic/), &
                                                count=(/nx_global,ny_global,1, 1/)), &
                                  'put var '//avail_hist_fields(n)%vname)
                    endif
                enddo ! k
            enddo ! ic
        endif
    enddo ! num_avail_hist_fields_4Ds

    work_gr(:,:) = c0
    work_g1(:,:) = c0

    do n = n4Dscum+1, n4Dbcum
        nn = n - n4Dscum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            if (my_task == master_task) then
                call check(nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid), &
                           'inq varid '//avail_hist_fields(n)%vname)
            endif
            do ic = 1, ncat_hist
                do k = 1, nzblyr
                    call gather_global(work_g1, a4Db(:,:,k,ic,nn,:), &
                                    master_task, distrb_info)
                    work_gr(:,:) = work_g1(:,:)
                    if (my_task == master_task) then
                        call check(nf90_put_var(ncid,varid,work_gr(:,:), &
                                                start=(/        1,        1,k,ic/), &
                                                count=(/nx_global,ny_global,1, 1/)), &
                                   'put var '//avail_hist_fields(n)%vname)
                    endif
                enddo ! k
            enddo ! ic
        endif
    enddo ! num_avail_hist_fields_4Db

    deallocate(work_gr)
    deallocate(work_g1)

end subroutine write_3d_and_4d_variables


subroutine write_coordinate_variables_parallel(ncid, coord_var, var_nz)

    integer, intent(in) :: ncid
    type(coord_attributes), dimension(ncoord), intent(in) :: coord_var
    type(coord_attributes), dimension(nvarz) :: var_nz

    integer :: varid
    integer :: iblk, i, k
    real(kind=dbl_kind), dimension(nx_block,ny_block, max_blocks) :: work1

    do i = 1,ncoord
        SELECT CASE (coord_var(i)%short_name)
        CASE ('TLON')
            ! Convert T grid longitude from -180 -> 180 to 0 to 360
            work1 = TLON*rad_to_deg + c360
            where (work1 > c360) work1 = work1 - c360
            where (work1 < c0 )  work1 = work1 + c360
        CASE ('TLAT')
            work1 = TLAT*rad_to_deg
        CASE ('ULON')
            work1 = ULON*rad_to_deg
        CASE ('ULAT')
            work1 = ULAT*rad_to_deg
        END SELECT

        call check(nf90_inq_varid(ncid, coord_var(i)%short_name, varid), &
                    'inq varid '//coord_var(i)%short_name)
        call put_2d_with_blocks(ncid, varid, coord_var(i)%short_name, work1)
    enddo

    ! Extra dimensions (NCAT, VGRD*)
    do i = 1, nvarz
        if (igrdz(i)) then
            call check(nf90_inq_varid(ncid, var_nz(i)%short_name, varid), &
                        'inq_varid '//var_nz(i)%short_name)
            SELECT CASE (var_nz(i)%short_name)
            CASE ('NCAT')
                call check(nf90_put_var(ncid, varid, hin_max(1:ncat_hist)), &
                            'put var NCAT')
            CASE ('VGRDi') ! index - needed for Met Office analysis code
                call check(nf90_put_var(ncid, varid, (/(k, k=1, nzilyr)/)), &
                            'put var VGRDi')
            CASE ('VGRDs') ! index - needed for Met Office analysis code
                call check(nf90_put_var(ncid, varid, (/(k, k=1, nzslyr)/)), &
                            'put var VGRDs')
            CASE ('VGRDb')
                call check(nf90_put_var(ncid, varid, (/(k, k=1, nzblyr)/)), &
                            'put var VGRDb')
            END SELECT
        endif
    enddo

end subroutine write_coordinate_variables_parallel


subroutine write_grid_variables_parallel(ncid, var, var_nverts)

    integer, intent(in) :: ncid
    type(req_attributes), dimension(nvar), intent(in) :: var
    type(coord_attributes), dimension(nvar_verts), intent(in) :: var_nverts

    real (kind=dbl_kind),  dimension(nx_block, ny_block, max_blocks) :: work1
    real (kind=dbl_kind),  dimension(nverts, nx_block, ny_block, max_blocks) :: work2

    integer :: iblk
    integer :: ilo, jlo, ihi, jhi, gilo, gjlo, gihi, gjhi
    integer, dimension(3) :: start, count
    type(block) :: the_block

    integer :: i
    integer :: varid

    if (igrd(n_tmask)) then
        call check(nf90_inq_varid(ncid, 'tmask', varid), 'inq var for tmask')
        call put_2d_with_blocks(ncid, varid, 'tmask', hm)
    endif

    if (igrd(n_blkmask)) then
        call check(nf90_inq_varid(ncid, 'blkmask', varid), 'inq var for blkmask')
        call put_2d_with_blocks(ncid, varid, 'blkmask', bm)
    endif

    do i = 3, nvar      ! note n_tmask=1, n_blkmask=2
        if (igrd(i)) then
            SELECT CASE (var(i)%req%short_name)
            CASE ('tarea')
                work1 = tarea
            CASE ('uarea')
                work1 = uarea
            CASE ('dxu')
                work1 = dxu
            CASE ('dyu')
                work1 = dyu
            CASE ('dxt')
                work1 = dxt
            CASE ('dyt')
                work1 = dyt
            CASE ('HTN')
                work1 = HTN
            CASE ('HTE')
                work1 = HTE
            CASE ('ANGLE')
                work1 = ANGLE
            CASE ('ANGLET')
                work1 = ANGLET
            END SELECT

            call check(nf90_inq_varid(ncid, var(i)%req%short_name, varid), &
                        'inq var '//var(i)%req%short_name)
            call put_2d_with_blocks(ncid, varid, var(i)%req%short_name, work1)
        endif
    enddo

    !----------------------------------------------------------------
    ! Write coordinates of grid box vertices
    !----------------------------------------------------------------

    if (f_bounds) then
        do i = 1, nvar_verts
            SELECT CASE (var_nverts(i)%short_name)
            CASE ('lont_bounds')
                work2(:, :, :, :) = lont_bounds(:, :, :, :)
            CASE ('latt_bounds')
                work2(:, :, :, :) = latt_bounds(:, :, :, :)
            CASE ('lonu_bounds')
                work2(:, :, :, :) = lonu_bounds(:, :, :, :)
            CASE ('latu_bounds')
                work2(:, :, :, :) = lonu_bounds(:, :, :, :)
            END SELECT

            call check(nf90_inq_varid(ncid, var_nverts(i)%short_name, varid), &
                       'inq varid '//var_nverts(i)%short_name)

            do iblk=1, nblocks
                the_block = get_block(blocks_ice(iblk), iblk)
                ilo = the_block%ilo
                jlo = the_block%jlo
                ihi = the_block%ihi
                jhi = the_block%jhi

                gilo = the_block%i_glob(ilo)
                gjlo = the_block%j_glob(jlo)
                gihi = the_block%i_glob(ihi)
                gjhi = the_block%j_glob(jhi)

                start = (/ 1, gilo, gjlo /)
                count = (/ nverts, gihi - gilo + 1, gjhi - gjlo + 1 /)
                call check(nf90_put_var(ncid, varid, &
                                        work2(1:nverts, ilo:ihi, jlo:jhi, iblk), &
                                        start=start, count=count), &
                            'grid vars _put '//trim(var_nverts(i)%short_name))
            enddo
        enddo
    endif

end subroutine write_grid_variables_parallel


subroutine write_2d_variables_parallel(ns, ncid)

    integer, intent(in) :: ns
    integer, intent(in) :: ncid

    integer :: varid
    integer :: n

    do n=1, num_avail_hist_fields_2D
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call check(nf90_inq_varid(ncid, avail_hist_fields(n)%vname, varid), &
                            'inq varid '//avail_hist_fields(n)%vname)
            call put_2d_with_blocks(ncid, varid, avail_hist_fields(n)%vname, &
                                    a2D(:, :, n, :))
        endif
    enddo ! num_avail_hist_fields_2D

end subroutine write_2d_variables_parallel



subroutine write_3d_and_4d_variables_parallel(ns, ncid)

    integer, intent(in) :: ns
    integer, intent(in) :: ncid

    integer :: varid
    integer :: status, n, nn, k, ic

    do n = n2D + 1, n3Dccum
        nn = n - n2D
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then

            call check(nf90_inq_varid(ncid,avail_hist_fields(n)%vname,varid), &
                        'inq varid '//avail_hist_fields(n)%vname)
            call put_3d_with_blocks(ncid, varid, avail_hist_fields(n)%vname, &
                                    ncat_hist, a3Dc(:, :, :, nn, :))
        endif
    enddo ! num_avail_hist_fields_3Dc


    do n = n3Dccum+1, n3Dzcum
        nn = n - n3Dccum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then

            call check(nf90_inq_varid(ncid, avail_hist_fields(n)%vname, varid), &
                       'inq varid '//avail_hist_fields(n)%vname)
            call put_3d_with_blocks(ncid, varid, avail_hist_fields(n)%vname, &
                                    nzilyr, a3Dz(:, :, :, nn, :))
        endif
    enddo ! num_avail_hist_fields_3Dz


    do n = n3Dzcum+1, n3Dbcum
        nn = n - n3Dzcum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call check(nf90_inq_varid(ncid, avail_hist_fields(n)%vname, varid), &
                       'inq varid '//avail_hist_fields(n)%vname)
            call put_3d_with_blocks(ncid, varid, avail_hist_fields(n)%vname, &
                                    nzblyr, a3Db(:, :, :, nn, :))
        endif
    enddo ! num_avail_hist_fields_3Db


    do n = n3Dbcum+1, n4Dicum
        nn = n - n3Dbcum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then

            call check(nf90_inq_varid(ncid, avail_hist_fields(n)%vname, varid), &
                       'inq varid for '//avail_hist_fields(n)%vname)
            call put_4d_with_blocks(ncid, varid, avail_hist_fields(n)%vname, &
                                    nzilyr, ncat_hist, a4Di(:, :, :, :, nn, :))
        endif
    enddo ! num_avail_hist_fields_4Di


    do n = n4Dicum+1, n4Dscum
        nn = n - n4Dicum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call check(nf90_inq_varid(ncid ,avail_hist_fields(n)%vname, varid), &
                       'inq varid for '//avail_hist_fields(n)%vname)
            call put_4d_with_blocks(ncid, varid, avail_hist_fields(n)%vname, &
                                    nzslyr, ncat_hist, a4Ds(:, :, :, :, nn, :))
        endif
    enddo ! num_avail_hist_fields_4Ds

    do n = n4Dscum+1, n4Dbcum
        nn = n - n4Dscum
        if (avail_hist_fields(n)%vhistfreq == histfreq(ns) .or. write_ic) then
            call check(nf90_inq_varid(ncid, avail_hist_fields(n)%vname, varid), &
                       'inq varid '//avail_hist_fields(n)%vname)

            call put_4d_with_blocks(ncid, varid, avail_hist_fields(n)%vname, &
                                    nzblyr, ncat_hist, a4Db(:, :, :, :, nn, :))
        endif
    enddo ! num_avail_hist_fields_4Db

end subroutine write_3d_and_4d_variables_parallel


subroutine put_2d_with_blocks(ncid, varid, var_name, data)

    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: var_name
    real(kind=dbl_kind), dimension(nx_block, ny_block, max_blocks), intent(in) :: data

    integer :: iblk
    integer :: ilo, jlo, ihi, jhi, gilo, gjlo, gihi, gjhi
    integer, dimension(2) :: start, count
    type(block) :: the_block

    do iblk=1, nblocks
        the_block = get_block(blocks_ice(iblk), iblk)
        ilo = the_block%ilo
        jlo = the_block%jlo
        ihi = the_block%ihi
        jhi = the_block%jhi

        gilo = the_block%i_glob(ilo)
        gjlo = the_block%j_glob(jlo)
        gihi = the_block%i_glob(ihi)
        gjhi = the_block%j_glob(jhi)

        start = (/ gilo, gjlo /)
        count = (/ gihi - gilo + 1, gjhi - gjlo + 1 /)
        call check(nf90_put_var(ncid, varid, data(ilo:ihi, jlo:jhi, iblk), &
                                start=start, count=count), &
                    'put_2d_with_blocks put '//trim(var_name))
    enddo

end subroutine put_2d_with_blocks

subroutine put_3d_with_blocks(ncid, varid, var_name, len_3dim, data)

    integer, intent(in) :: ncid, varid, len_3dim
    character(len=*), intent(in) :: var_name
    real(kind=dbl_kind), dimension(nx_block, ny_block, len_3dim, max_blocks), intent(in) :: data

    integer :: iblk
    integer :: ilo, jlo, ihi, jhi, gilo, gjlo, gihi, gjhi
    integer, dimension(3) :: start, count
    type(block) :: the_block

    do iblk=1, nblocks
        the_block = get_block(blocks_ice(iblk), iblk)
        ilo = the_block%ilo
        jlo = the_block%jlo
        ihi = the_block%ihi
        jhi = the_block%jhi

        gilo = the_block%i_glob(ilo)
        gjlo = the_block%j_glob(jlo)
        gihi = the_block%i_glob(ihi)
        gjhi = the_block%j_glob(jhi)

        start = (/ gilo, gjlo, 1 /)
        count = (/ gihi - gilo + 1, gjhi - gjlo + 1, len_3dim /)
        call check(nf90_put_var(ncid, varid, &
                                data(ilo:ihi, jlo:jhi, 1:len_3dim, iblk), &
                                start=start, count=count), &
                    'put_3d_with_blocks put '//trim(var_name))
    enddo

end subroutine put_3d_with_blocks


subroutine put_4d_with_blocks(ncid, varid, var_name, len_3dim, len_4dim, data)

    integer, intent(in) :: ncid, varid, len_3dim, len_4dim
    character(len=*), intent(in) :: var_name
    real(kind=dbl_kind), dimension(nx_block, ny_block, len_3dim, &
                                   len_4dim, max_blocks), intent(in) :: data

    integer :: iblk
    integer :: ilo, jlo, ihi, jhi, gilo, gjlo, gihi, gjhi
    integer, dimension(4) :: start, count
    type(block) :: the_block

    do iblk=1, nblocks
        the_block = get_block(blocks_ice(iblk), iblk)
        ilo = the_block%ilo
        jlo = the_block%jlo
        ihi = the_block%ihi
        jhi = the_block%jhi

        gilo = the_block%i_glob(ilo)
        gjlo = the_block%j_glob(jlo)
        gihi = the_block%i_glob(ihi)
        gjhi = the_block%j_glob(jhi)

        start = (/ gilo, gjlo, 1, 1 /)
        count = (/ gihi - gilo + 1, gjhi - gjlo + 1, len_3dim, len_4dim /)
        call check(nf90_put_var(ncid, varid, &
                                data(ilo:ihi, jlo:jhi, 1:len_3dim, 1:len_4dim, iblk), &
                                start=start, count=count), &
                    'put_4d_with_blocks put '//trim(var_name))
    enddo

end subroutine put_4d_with_blocks


end module ice_history_write
