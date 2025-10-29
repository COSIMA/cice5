!=======================================================================
!
!BOP
!
! !MODULE: ice_coupling - contains coupling related routines used by Met Office
!
! !DESCRIPTION:
!
!  Contains routines relating to coupling fields used by Met Office
!
! !REVISION HISTORY:
!  SVN:$Id: 
!
!  authors: Alison McLaren, Met Office
!  Feb 2014: Amended by Alex West for use in CICE 5.0.  
!
! !INTERFACE:
!
      module ice_coupling
!
! !USES:
!
      use ice_constants
      use ice_kinds_mod
!
!EOP
!     
      implicit none
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: sfcflux_to_ocn, set_sfcflux, top_layer_Tandk_init, top_layer_Tandk_run
!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: sfcflux_to_ocn
!
! !DESCRIPTION:
!
! If surface heat fluxes are provided to CICE instead of CICE calculating
! them internally (i.e. .not. calc_Tsfc), then these heat fluxes can 
! be provided at points which do not have ice.  (This is could be due to
! the heat fluxes being calculated on a lower resolution grid or the
! heat fluxes not recalculated at every CICE timestep.)  At ice free points, 
! conserve energy and water by passing these fluxes to the ocean.
!
! !INTERFACE:
!
       subroutine sfcflux_to_ocn(nx_block,   ny_block,     &
                                 tmask,      aice,         &
                                 fsurfn_f,   flatn_f,      &
                                 fresh,      fhocn)
!
! !REVISION HISTORY:
!
! authors: A. McLaren, Met Office
!
! !USES:
!
      use ice_domain_size, only: ncat
!
! !INPUT/OUTPUT PARAMETERS:
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
!
!EOP
!
!#ifdef CICE_IN_NEMO
#ifdef ACCESS
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

! If model is not calculating surface temperature, set the surface
! flux values using values read in from forcing data or supplied via
! coupling (stored in ice_flux).
!
! If CICE is running in NEMO environment, convert fluxes from GBM values 
! to per unit ice area values. If model is not running in NEMO environment, 
! the forcing is supplied as per unit ice area values.
!
! authors Alison McLaren, Met Office

      subroutine set_sfcflux (nx_block,  ny_block, &
                              n,         iblk,     &
                              icells,              & 
                              indxi,     indxj,    &
                              aicen,               &
                              flatn,               &
			      fsensn,              &
                              fsurfn,              &
                              fcondtopn)

      use ice_fileunits, only: nu_diag
      use ice_flux, only: fsurfn_f, fcondtopn_f, flatn_f, fsensn_f

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         n,                  & ! thickness category index
         iblk,               & ! block index
         icells                ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      ! ice state variables
      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aicen           ! concentration of ice

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out):: &
         flatn       , & ! latent heat flux   (W/m^2) 
	 fsensn      , & ! sensible heat flux (W/m^2)
         fsurfn      , & ! net flux to top surface, not including fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij              ! horizontal indices, combine i and j loops

      real (kind=dbl_kind)  :: &
         raicen          ! 1 or 1/aicen

      logical (kind=log_kind) :: &
         extreme_flag    ! flag for extreme forcing values

      logical (kind=log_kind), parameter :: & 
         extreme_test=.false. ! test and write out extreme forcing data

         raicen        = c1
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

!#ifdef CICE_IN_NEMO
#ifdef ACCESS
!----------------------------------------------------------------------
! Convert fluxes from GBM values to per ice area values when 
! running in NEMO environment.  (When in standalone mode, fluxes
! are input as per ice area.)
!----------------------------------------------------------------------
            raicen        = c1 / aicen(i,j)
#endif
            fsurfn(i,j)   = fsurfn_f(i,j,n,iblk)*raicen
            fcondtopn(i,j)= fcondtopn_f(i,j,n,iblk)*raicen
            flatn(i,j)    = flatn_f(i,j,n,iblk)*raicen
	    fsensn(i,j)   = fsensn_f(i,j,n,iblk)*raicen

         enddo

!----------------------------------------------------------------
! Flag up any extreme fluxes
!---------------------------------------------------------------

         if (extreme_test) then
            extreme_flag = .false.

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)         

               if (fcondtopn(i,j) < -100.0_dbl_kind & 
                     .or. fcondtopn(i,j) > 20.0_dbl_kind) then
                  extreme_flag = .true.
               endif

               if (fsurfn(i,j) < -100.0_dbl_kind & 
                    .or. fsurfn(i,j) > 80.0_dbl_kind) then
                  extreme_flag = .true.
               endif

               if (flatn(i,j) < -20.0_dbl_kind & 
                     .or. flatn(i,j) > 20.0_dbl_kind) then
                  extreme_flag = .true.
               endif

            enddo  ! ij

            if (extreme_flag) then
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)         

                  if (fcondtopn(i,j) < -100.0_dbl_kind & 
                       .or. fcondtopn(i,j) > 20.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -100 > fcondtopn > 20'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,fcondtopn = ', & 
                       i,j,n,iblk,aicen(i,j),fcondtopn(i,j)
                  endif

                  if (fsurfn(i,j) < -100.0_dbl_kind & 
                       .or. fsurfn(i,j) > 80.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -100 > fsurfn > 40'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,fsurfn = ', & 
                        i,j,n,iblk,aicen(i,j),fsurfn(i,j)
                  endif

                  if (flatn(i,j) < -20.0_dbl_kind & 
                       .or. flatn(i,j) > 20.0_dbl_kind) then
                     write(nu_diag,*) & 
                       'Extreme forcing: -20 > flatn > 20'
                     write(nu_diag,*) & 
                       'i,j,n,iblk,aicen,flatn = ', & 
                        i,j,n,iblk,aicen(i,j),flatn(i,j)
                  endif

               enddo  ! ij
      
            endif  ! extreme_flag
         endif     ! extreme_test    


      end subroutine set_sfcflux 

!=======================================================================
!=======================================================================
!BOP
!
! !ROUTINE: top_layer_Tandk_init
!
! !DESCRIPTION:
!
! Hacked version to be called upon initialisation (when we're not
! parallelised)
! Calculate the top layer temperature and conductivity for passing
! to atmosphere model or calculating Tsfc explicitly.
!
! This routine is only called if calc_Tsfc = F and heat_capacity = T.
!
! !REVISION HISTORY:
!
! authors: Alison McLaren, Met Office
! Feb 2014: Modified by Alex West to work in CICE 5.0
!
! !INTERFACE:

      subroutine top_layer_Tandk_init
!
! !USES:
!
      use ice_blocks
      use ice_constants
      use ice_domain, only: nblocks
      use ice_domain_size
      use ice_fileunits, only: nu_diag
      use ice_flux, only: Tn_top, keffn_top
      use ice_itd, only: hs_min
      use ice_state, only: aicen, vicen, vsnon, trcrn, nt_qice, nt_sice, nt_qsno
      use ice_therm_mushy, only: liquidus_temperature_mush
      use ice_therm_shared, only: calculate_ki_from_Tin, calculate_Tin_from_qin, ktherm, Tmlt, conduct
      
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      
      integer (kind=int_kind) :: & 
         iblk        , & ! block index 
         n           , & ! thickness category index
         i,j             ! horizontal indices

      real (kind=dbl_kind) ::  &
         rnslyr      , & ! real(nslyr)
         rnilyr      , & ! real(nilyr)
	 hs1         , & ! thickness of top snow layer
	                 ! (so we know whether the top layer is snow or ice)
	 hi1         , & ! thickness of top ice layer
	 Tmlt1       , & ! melting temperature of top ice layer
         ki              ! top ice layer conductivity
	 
	 
      real (kind=dbl_kind) :: &
         ki_hold, &
	 Ti_hold         ! debugging variables
	 
      keffn_top(:,:,:,:) = c0   ! initialise
      Tn_top(:,:,:,:)    = c0   
      rnslyr = real(nslyr,kind=dbl_kind)      
      rnilyr = real(nilyr,kind=dbl_kind)    
              
      do iblk = 1, max_blocks
	 do n = 1, ncat
            do j = 1, ny_block
               do i = 1, nx_block

        	  if (aicen(i,j,n,iblk) > puny) then

        	     hs1 = vsnon(i,j,n,iblk)/(aicen(i,j,n,iblk)*rnslyr)

		     if (hs1 > hs_min/rnslyr) then

                	!snow is top layer
                	Tn_top(i,j,n,iblk)    = (Lfresh + trcrn(i,j,nt_qsno,n,iblk) / rhos)/cp_ice
                	keffn_top(i,j,n,iblk) = c2 * ksno / hs1   

		     else
                	!ice is top layer
                	hi1 = vicen(i,j,n,iblk)/(aicen(i,j,n,iblk)*rnilyr)
        		if (ktherm == 2) then
        		   Tmlt1 = liquidus_temperature_mush(trcrn(i,j,nt_sice,n,iblk))
        		else
        		   Tmlt1 = - trcrn(i,j,nt_sice,n,iblk) * depressT
        		endif

                	Tn_top(i,j,n,iblk)    = & 
                        	  calculate_Tin_from_qin(trcrn(i,j,nt_qice,n,iblk),Tmlt1)
			Ti_hold = calculate_Tin_from_qin(trcrn(i,j,nt_qice,n,iblk),Tmlt1)
			ki_hold = calculate_ki_from_Tin(Tn_top(i,j,n,iblk),trcrn(i,j,nt_sice,n,iblk))
                	ki = calculate_ki_from_Tin(Tn_top(i,j,n,iblk),trcrn(i,j,nt_sice,n,iblk))
			keffn_top(i,j,n,iblk) = c2 * ki / hi1   
		     endif

		  endif     ! aice > puny   
               enddo     ! i
	    enddo     ! i
	 enddo     ! n
      enddo     ! iblk
            
      end subroutine top_layer_Tandk_init

!=======================================================================
!BOP
!
! !ROUTINE: top_layer_Tandk_run
!
! !DESCRIPTION:
!
! Calculate the top layer temperature and conductivity for passing
! to atmosphere model or calculating Tsfc explicitly.
!
! This routine is only called if calc_Tsfc = F and heat_capacity = T.
!
! !REVISION HISTORY:
!
! authors: Alison McLaren, Met Office
! Feb 2014: Modified by Alex West to work in CICE 5.0
!
! !INTERFACE:

      subroutine top_layer_Tandk_run (iblk)
!
! !USES:
!
      use ice_blocks
      use ice_constants
      use ice_domain, only: nblocks
      use ice_domain_size
      use ice_fileunits, only: nu_diag
      use ice_flux, only: Tn_top, keffn_top
      use ice_itd, only: hs_min
      use ice_state, only: aicen, vicen, vsnon, trcrn, nt_qice, nt_sice, nt_qsno
      use ice_therm_mushy, only: liquidus_temperature_mush
      use ice_therm_shared, only: calculate_ki_from_Tin, calculate_Tin_from_qin, ktherm, Tmlt, conduct
      
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      
      integer (kind=int_kind), intent(in) :: & 
         iblk            ! block index 
      integer (kind=int_kind) :: & 
         n           , & ! thickness category index
         i,j             ! horizontal indices

      real (kind=dbl_kind) ::  &
         rnslyr      , & ! real(nslyr)
         rnilyr      , & ! real(nilyr)
	 hs1         , & ! thickness of top snow layer
	                 ! (so we know whether the top layer is snow or ice)
	 hi1         , & ! thickness of top ice layer
	 Tmlt1       , & ! melting temperature of top ice layer
         ki              ! top ice layer conductivity
	 
	 
      real (kind=dbl_kind) :: &
         ki_hold, &
	 Ti_hold         ! debugging variables
	 
      keffn_top(:,:,:,:) = c0   ! initialise
      Tn_top(:,:,:,:)    = c0   
      rnslyr = real(nslyr,kind=dbl_kind)      
      rnilyr = real(nilyr,kind=dbl_kind)    
              
      do n = 1, ncat
         do j = 1, ny_block
            do i = 1, nx_block

               if (aicen(i,j,n,iblk) > puny) then

        	  hs1 = vsnon(i,j,n,iblk)/(aicen(i,j,n,iblk)*rnslyr)

		  if (hs1 > hs_min/rnslyr) then

                     !snow is top layer
                     Tn_top(i,j,n,iblk)    = (Lfresh + trcrn(i,j,nt_qsno,n,iblk) / rhos)/cp_ice
                     keffn_top(i,j,n,iblk) = c2 * ksno / hs1   

		  else
                     !ice is top layer
                     hi1 = vicen(i,j,n,iblk)/(aicen(i,j,n,iblk)*rnilyr)
        	     if (ktherm == 2) then
        		Tmlt1 = liquidus_temperature_mush(trcrn(i,j,nt_sice,n,iblk))
        	     else
        		Tmlt1 = - trcrn(i,j,nt_sice,n,iblk) * depressT
        	     endif

                     Tn_top(i,j,n,iblk)    = & 
                               calculate_Tin_from_qin(trcrn(i,j,nt_qice,n,iblk),Tmlt1)
		     Ti_hold = calculate_Tin_from_qin(trcrn(i,j,nt_qice,n,iblk),Tmlt1)
		     ki_hold = calculate_ki_from_Tin(Tn_top(i,j,n,iblk),trcrn(i,j,nt_sice,n,iblk))
                     ki = calculate_ki_from_Tin(Tn_top(i,j,n,iblk),trcrn(i,j,nt_sice,n,iblk))
		     keffn_top(i,j,n,iblk) = c2 * ki / hi1   
		  endif

	       endif     ! aice > puny   
            enddo     ! i
	 enddo     ! i
      enddo     ! n
            
      end subroutine top_layer_Tandk_run

!======================================================================

      end module ice_coupling

!======================================================================
