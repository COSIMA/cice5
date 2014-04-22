module cpl_netcdf_setup

use ice_kinds_mod
use ice_calendar, only: idate, sec

implicit none

include 'netcdf.inc'	!define nf_* 
!use netcdf             !define nf90_*

integer(kind=int_kind) :: pLonDimId, pLatDimId, timeDimId, pDepDimId
contains 

!=========================================================================
subroutine ncheck(status) 

implicit none

integer(kind=int_kind), intent(in) :: status

if (status /= nf_noerr) then
  write(*,'(/a)')   'error - from NetCDF library'
  write(*,'(a/)')   trim(nf_strerror(status))
  stop
end if
end subroutine ncheck

!============================================================================
subroutine read_nc(vout,nx,ny,vname,fname)
!
!read in a 2D (X-Y) array with double precision !*** not needed anymore! *** 
!  instead, the more concise 'ice_read_nc' is use for all nc data reading! 
!
!Note the compiling option ("idbl4" I think) does not allow the models to 
!     read in any single precision variables from a netcdf file such as the
!     A2I_..nc etc. So we must feed the model with dbl precision vars in nc
!     files. Similarly, we need write out dbl precision fields, NOT single
!     precision arrays, into ncfiles! 
!     --- the above note is perhaps not true! ---
!     if we define "real*4 vout", then nf_get_vara_real(...,vout) may work. 

implicit none

real(kind=dbl_kind), dimension(nx,ny), intent(out) :: vout 
integer(kind=int_kind), intent(in) :: nx, ny
character(len=*), intent(in) :: vname,fname

!real, dimension(nx,ny) :: vtmp      !single precision for 'conversion'
integer(kind=int_kind) :: varid, ndim, ncid
integer(kind=int_kind), dimension(:), allocatable :: count, start

!Open file for read access
call ncheck( nf_open(fname, nf_nowrite, ncid) )

!Get variable ID
call ncheck(nf_inq_varid(ncid, vname, varid))

!Get number of dimensions
call ncheck(nf_inq_varndims(ncid, varid, ndim))

!Allocate count and start
allocate (count(ndim), start(ndim))

!Get the 2D array data out of a ndim-D variable
if (ndim == 2) then       !currently only have ndim = 2!
  start = (/ 1, 1 /)
  count = (/ nx, ny /)
else if (ndim == 3) then
  start = (/ 1, 1, 1/)
  count = (/ nx, ny, 1 /)
else  !* ndim = 4 *!
  start = (/ 1, 1, 1, 1 /)
  count = (/ nx, ny, 1, 1 /)
endif

!call ncheck(nf_get_vara_real(ncid, varid, start, count, vtmp))
call ncheck(nf_get_vara_double(ncid, varid, start, count, vout))

!Close file
call ncheck( nf_close(ncid) )

!vout = vtmp

return
end subroutine read_nc 

!===========================================================================
subroutine create_ncfile(ncfile, ncid, ii, jj, kk, ll, ilout)
!
!to create 2, 3,or 4D ncfile, depending on optional args (kk,ll)! 
!

implicit none

integer(kind=int_kind), intent(in) :: ii,jj     	!x, y dimension size
!!!integer(kind=int_kind), optional, intent(in) :: kk, ll !z, t dimension size
!!!integer(kind=int_kind), optional, intent(in) :: ilout  !format io file id
! * 'optional' att can NOT be with 'intent(in) *' 
integer(kind=int_kind), optional :: kk, ll !z, t dimension size
integer(kind=int_kind), optional :: ilout  !format io file id
character(len=*), intent(in) :: ncfile  
integer(kind=int_kind), intent(out) :: ncid

if (present(ilout)) write(ilout,*) 'creating a new netcdf file: ',ncfile

!create a new NetCDF and define the grid:
call ncheck(nf_create(trim(ncfile),nf_write,ncid))

!define the dimensions
if (present(ll)) call ncheck(nf_def_dim(ncid,"time", nf_unlimited,  timeDimId))
if (present(kk)) call ncheck(nf_def_dim(ncid,"nz", kk,  pDepDimId))
call ncheck(nf_def_dim(ncid, "ny", jj,  pLatDimId))
call ncheck(nf_def_dim(ncid, "nx", ii,  pLonDimId))

!end of the definition phase
call ncheck(nf_enddef(ncid))

!close NetCDF file
!call ncheck(nf_close(ncid))
!do NOT close it here!
write(*,'(2a)') 'ncfile created: ',trim(ncfile)

return
end subroutine create_ncfile

!===========================================================================
subroutine open_ncfile(ncfile, ncid, ii, jj, kk, ll, ilout)
!
!to open exsiting ncfile to modify it
!

implicit none

integer(kind=int_kind), intent(in) :: ii,jj             !x, y dimension size
!!!integer(kind=int_kind), optional, intent(in) :: kk, ll !z, t dimension size
!!!integer(kind=int_kind), optional, intent(in) :: ilout  !format io file id
! * 'optional' att can NOT be with 'intent(in) *'
integer(kind=int_kind), optional :: kk, ll !z, t dimension size
integer(kind=int_kind), optional :: ilout  !format io file id
character(len=*), intent(in) :: ncfile
integer(kind=int_kind), intent(out) :: ncid

if (present(ilout)) write(ilout,*) 'opening a existing netcdf file: ',ncfile

!create a new NetCDF and define the grid:
call ncheck(nf_open(trim(ncfile),nf_write,ncid))

!close NetCDF file
!call ncheck(nf_close(ncid))
!do NOT close it here!
write(*,'(2a)') 'ncfile opened: ',trim(ncfile)

return
end subroutine open_ncfile

!===========================================================================
subroutine create_nc(ncfile,ncid,ii,jj)
!
!not needed anymore 'cos its function is covered by 'create_ncfile'
! 

implicit none

integer(kind=int_kind), intent(in) :: ii,jj
character(len=*), intent(in) :: ncfile
integer(kind=int_kind), intent(out) :: ncid
 
print *
print *, 'creating a new netcdf file: ',ncfile

!create a new NetCDF and define the grid:
call ncheck(nf_create(trim(ncfile),nf_write,ncid))

!define the dimensions
call ncheck(nf_def_dim(ncid,"ny", jj,  pLatDimId))
call ncheck(nf_def_dim(ncid,"nx", ii,  pLonDimId))
!B: note the sizes of dimensions ('jj' and 'ii') can't be missing!

!end of the definition phase
call ncheck(nf_enddef(ncid))
 
!close NetCDF file
!call ncheck(nf_close(ncid))
!do NOT close it here!
write(*,'(2a)') 'ncfile created: ',trim(ncfile)

return
end subroutine create_nc

!===========================================================================
subroutine write_nc_1Dtime(vin, nt, vname, ncid)

implicit none

integer(kind=int_kind), intent(in) :: ncid,nt
integer(kind=int_kind) :: varid, ncstatus
integer(kind=int_kind), dimension(1:6) :: adate
!real(kind=dbl_kind), dimension(nt), intent(in) :: vin
real, intent(in) :: vin 	
! NOTE here real is default real*8 (which is actually the same as dbl_kind!)
! somehow the netcdf lib used here takes 'real' as real*4. therefore we need: 
real*4 :: vtmp
character(len=*), intent(in) :: vname
character*80 ctimeatt

vtmp = real(vin)  ! isn't this real here real*8 ?
print *, 'write_nc_1Dtime: time to write field -- ', vtmp, vname
print *, 'write_nc_1Dtime: idate, sec          -- ', idate, sec

ncstatus=nf_inq_varid(ncid,vname,varid)

if (ncstatus/=nf_noerr) then
  adate(1) = idate/10000
  adate(2) = (idate - (idate/10000)*10000)/100
  adate(3) = idate - (idate/100)*100
  adate(4:6) = 0  !OK for 'whole-day' runs               
  call ncheck(nf_redef(ncid))
  call ncheck(nf_def_var(ncid,trim(vname),nf_real, 1, timeDimId, varid))
  write(ctimeatt, &
      '("seconds since ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') adate(:)
  !ctimeatt="hours since 0000-01-01 00:00:00"
  call ncheck(nf_put_att_text(ncid,varid,"units",len_trim(ctimeatt),trim(ctimeatt)))
  call ncheck(nf_enddef(ncid))
end if

!write values into the 1D array
!call ncheck(nf_put_vara_real(ncid,varid,nt,1,vin))
call ncheck(nf_put_vara_real(ncid,varid,nt,1,vtmp))
!B: must indicate the start point and number of the record, ie., nt and 1!  

return
end subroutine write_nc_1Dtime

!===========================================================================
subroutine write_nc2D(ncid, vname, vin, prcn, nx, ny, istep, ilout)
!
!to output a 2D array into a 3D field (with time dimension) 
!with either single or double precisioin depending on argumnet 'prcn'!
!

implicit none

integer(kind=int_kind), intent(in) :: ncid
integer(kind=int_kind), intent(in) :: prcn	!precision choice (1/2: signle/double)
character(len=*), intent(in) :: vname
integer(kind=int_kind), intent(in) :: nx, ny
integer(kind=int_kind), intent(in) :: istep	!position in the time dim (No of record) 
!!!integer(kind=int_kind), intent(in), optional :: ilout 
integer(kind=int_kind), optional :: ilout
real(kind=dbl_kind), dimension(nx,ny), intent(in) :: vin

integer(kind=int_kind) :: varid, ncstatus 
real*4, dimension(nx,ny) :: vtmp   !single precision

if (present(ilout)) write(ilout,*) 'write_nc2D: handling var *** ',vname, ' rec: ', istep

ncstatus=nf_inq_varid(ncid,vname,varid)
if (ncstatus/=nf_noerr) then
  call ncheck(nf_redef(ncid))
  if (prcn == 1) then
    call ncheck(nf_def_var(ncid,trim(vname),nf_real, 3, &
            (/pLonDimId, pLatDimId, timeDimId/),varid))
  else
    call ncheck(nf_def_var(ncid,trim(vname),nf_double, 3, &
            (/pLonDimId, pLatDimId, timeDimId/),varid))
  endif
  call ncheck(nf_enddef(ncid))
  if (present(ilout)) write(ilout,*) 'write_nc2D: defined new var ***', vname 
else
  if (present(ilout)) write(ilout,*) 'write_nc2D: found   old var ***', vname
end if

select case(prcn)
  case (1)
    vtmp = real(vin) !dbl precision to single precision
    call ncheck(nf_put_vara_real(ncid,varid,(/1,1,istep/),(/nx,ny,1/),vtmp))
  case default    !case (2)
    call ncheck(nf_put_vara_double(ncid,varid,(/1,1,istep/),(/nx,ny,1/),vin))
end select

return
end subroutine write_nc2D

!===========================================================================
subroutine modify_nc2D(ncid, vname, vin, prcn, nx, ny, istep, ilout)
!
!to modify a 2D array into an existing 3D field (with time dimension)
!with either single or double precisioin depending on argumnet 'prcn'!
!

implicit none

integer(kind=int_kind), intent(in) :: ncid
integer(kind=int_kind), intent(in) :: prcn      !precision choice (1/2: signle/double)
character(len=*), intent(in) :: vname
integer(kind=int_kind), intent(in) :: nx, ny
integer(kind=int_kind), intent(in) :: istep     !position in the time dim (No of record)
!!!integer(kind=int_kind), intent(in), optional :: ilout
integer(kind=int_kind), optional :: ilout
real(kind=dbl_kind), dimension(nx,ny), intent(in) :: vin

integer(kind=int_kind) :: varid, ncstatus
real*4, dimension(nx,ny) :: vtmp   !single precision

if (present(ilout)) write(ilout,*) 'modify_nc2D: handling var *** ',vname, ' rec: ', istep

ncstatus=nf_inq_varid(ncid,vname,varid)
if (ncstatus/=nf_noerr) then
  if (present(ilout)) write(ilout,*) 'modify_nc2D: Error- not found old var ***', vname
  stop
else
  if (present(ilout)) write(ilout,*) 'modify_nc2D: found old var ***', vname
end if

select case(prcn)
  case (1)
    vtmp = real(vin) !dbl precision to single precision
    !call ncheck(nf_put_vara_real(ncid,varid,(/1,1,istep/),(/nx,ny,1/),vtmp))
    call ncheck(nf_put_vara_real(ncid,varid,(/1,1/),(/nx,ny/),vtmp))
  case default    !case (2)
    !call ncheck(nf_put_vara_double(ncid,varid,(/1,1,istep/),(/nx,ny,1/),vin))
    call ncheck(nf_put_vara_double(ncid,varid,(/1,1/),(/nx,ny/),vin))
end select

return
end subroutine modify_nc2D

!===========================================================================
subroutine write_nc(vin, nx, ny, vname, ncid)
!
! *** its function is covered by write_nc2D and thus not used any more. ***
!

implicit none

integer(kind=int_kind), intent(in) :: nx, ny, ncid
integer(kind=int_kind) :: varid
real(kind=dbl_kind), dimension(nx,ny), intent(in) :: vin
character(len=*), intent(in) :: vname

!switch to define mode
call ncheck(nf_redef(ncid))

!define (2D) variable to be written
call ncheck(nf_def_var(ncid,trim(vname),nf_double, 2, &
            (/pLonDimId, pLatDimId/),varid))
!B: here '2' indicates the dimension of the to-be-written variable 'vname',
!   and (/pLonDimId, pLatDimId/) the dimention.
!   'nf_real' determines the output precision!

!leave define mode
call ncheck(nf_enddef(ncid))

!get varId and write to array
call ncheck(nf_inq_varid(ncid,trim(vname),varid))

!write values into the 2D array
call ncheck(nf_put_vara_double(ncid,varid,(/1,1/),(/nx,ny/),vin))
!B: must indicate the start point and number of the record (ie. count as in
!   routine read_nc above) (/1,1/) and (/nx,ny/) for the to-be-written vin!

return
end subroutine write_nc

!===========================================================================

end module cpl_netcdf_setup
