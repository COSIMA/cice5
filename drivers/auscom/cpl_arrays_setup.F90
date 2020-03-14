!============================================================================
!
module cpl_arrays_setup
!
! It's designed to include the following 'possible' coupling fields at the 
! air-ice-sea surface:
!
! A> atm==>ice 
!                                       
! (1) 2m air temp (K)                      		tair0 
! (2) 10m 'zonal' wind speed (m/s)         		uwnd0 
! (3) 10m 'meridional' wind speed (m/s)     		vwnd0 
! (4) shortwave radiation (down, J/m^2)     		swflx0
! (5) longwave radiation  (down, J/m^2)     		lwflx0
! (6) 2m air specific humidity (kg/kg)      		qair0 
! (7) rainfall rate     (kg/m^2/s)        		rain0 
! (8) snowfall rate     (kg/m^2/s)                      snow0
! (9) pressure            (Pa)                		press0
! (10)runof               (kg/m^2/s)          		runof0 
!
! B> ocn==>ice 
!                          
! (1) sea surface temperature  (K)              	ssto 
! (2) sea surface salinity   (psu)             	 	ssso
! (3) zonal water speed      (m/s)              	ssuo
! (4) meridional water speed (m/s)              	ssvo
! (5) sea surface gradient (zonal)     (m/m)   	 	sslx
! (6) sea surface gradient (meridional)(m/m)    	ssly
! (7) potential ice frm/mlt heatflux (W/m^2)    	pfmice
!
! D> ice==>ocn
!             
! (1) ice-ocean stress, x-direction (kg/m s^2)    	iostrsu
! (2) ice-ocean stress, y-direction (kg/m s^2)    	iostrsv
! (3) fresh water flux to ocean--rain (kg/m^2/s)      	iorain
! (4) fresh water flux to ocean--snow (kg/m^2/s)	iosnow
! (5) salt flux to ocean (kg/m^2/s)               	iostflx
! (6) 'net' heat flux to ocean (W/m^2)            	iohtflx
!     *(note word 'net' is misleading!) it is actually ice
!     *'melt' heatflux into ocean. (ref: ice_coupling.F, 
!     *it says:
!     *'buffs(n,index_i2c_Fioi_melth) = fhnet(i,j) 
!     *                          ! hf from melting'
! (7) shortwave penetrating to ocean (W/m^2)      	ioswflx
!     Also, we let the following 'atmospheric fluxes' 
!     (some maybe calculated in cice) be passed into ocean:
! (8) latent heat flux (W/m^2, positive out of ocean)   ioqflux
! (9) sensible heat flux (W/m^2, postive out of ocean  	ioshflx
!--- note sensible/latent heatfluxes are calculated in cice being 
!    positive into ocean! they must change sign before sent to 
!    mom4 (which requires these 2 item as positive out of ocean!) 
!    this is done in routine" get_i2o_fluxes"
!
!(10) long wave radiation                         	iolwflx
!(11) runoff                                      	iorunof
!(12) pressure                                    	iopress
!(13) ice concentration (fraction)			ioaice
!
! Seperate ice melting/forcation associated water fluxes from the rainfall field:
!
!(14) ice melt waterflux				iomelt
!(15) ice form waterflux				ioform
!(16) land ice waterflux                iolicefw
!(17) land ice heatflux                 iolicefh
!
!
! 18 in, 18 out => thus we set jpfldout=18, jpfldin=18 (in cpl_parameters)! 
!
!----------------------------------------------------------------------------
!Note: 
!   1. runoff forcing (read in from NCEP2 by matm and passed in cice) is NOT 
!      used! instead, the core2 runoff remapped onto AusCOM grid (by Marsmald)
!      should be read in by cice and passed directly into mom4p1! (17/07/2009)
!      * This approach is only 'temporarily suitable' for AusCOM. Eventually
!      * we have to work out a remapping scheme for real a2i coupling. 
!   2. ice coverage is passed into mom4 
!============================================================================

use ice_kinds_mod

implicit none

! Fields in
real(kind=dbl_kind), dimension(:,:,:), allocatable :: &   !from atm
    tair0, swflx0, lwflx0, uwnd0, vwnd0, qair0, rain0, snow0 & !(for ice)
   ,runof0, press0, calv0                                      !(for ocn)    

real(kind=dbl_kind), dimension(:,:,:), allocatable :: runof, calv, press

! CORE runoff remapped onto the AusCOM grid
real(kind=dbl_kind), dimension(:,:,:), allocatable :: & 
    core_runoff

real(kind=dbl_kind), dimension(:,:,:), allocatable :: &   !from ocn
    ssto,  ssso,   ssuo,   ssvo,   sslx,  ssly,  pfmice  
real(kind=dbl_kind), dimension(:,:), allocatable :: gwork
    !global domain work array, used for coupling data passing and global data output. 
real(kind=dbl_kind), dimension(:,:,:), allocatable :: vwork  

real(kind=dbl_kind),dimension(:,:,:), allocatable :: &     !to ocn (time averaged)
    iostrsu, iostrsv, iorain, iosnow, iostflx, iohtflx, ioswflx &
   ,ioqflux, ioshflx, iolwflx, iorunof, iopress, ioaice &
   ,iomelt, ioform, iolicefw, iolicefh
real(kind=dbl_kind),dimension(:,:,:), allocatable :: &     !to ocn (temporary)
    tiostrsu, tiostrsv, tiorain, tiosnow, tiostflx, tiohtflx, tioswflx &
   ,tioqflux, tioshflx, tiolwflx, tiorunof, tiopress, tioaice &
   ,tiomelt, tioform, tiolicefw, tiolicefh

! other stuff 
real(kind=dbl_kind),dimension(:,:,:), allocatable :: & 
    sicemass
! 20100111: for gfdl surface roughness calculation
real(kind=dbl_kind),dimension(:,:,:), allocatable :: &
    u_star0, rough_mom0, rough_heat0, rough_moist0 

!=================================================================================
end module cpl_arrays_setup

