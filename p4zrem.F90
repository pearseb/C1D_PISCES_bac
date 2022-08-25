MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zlim
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   REAL(wp), PUBLIC ::   xremikc    !: remineralisation rate of DOC 
   REAL(wp), PUBLIC ::   xremikn    !: remineralisation rate of DON 
   REAL(wp), PUBLIC ::   xremikp    !: remineralisation rate of DOP 
   REAL(wp), PUBLIC ::   xremik     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   nitrif     !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::   xsirem     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsiremlab  !: fast remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsilab     !: fraction of labile biogenic silica 
   REAL(wp), PUBLIC ::   feratb     !: Fe/C quota in bacteria
   REAL(wp), PUBLIC ::   xkferb     !: Half-saturation constant for bacteria Fe/C
   REAL(wp), PUBLIC ::   aoa_pocoef !: Rate of oxygen uptake by AOA via diffusion 
   REAL(wp), PUBLIC ::   aoa_mumax  !: Maximum growth rate of AOA
   REAL(wp), PUBLIC ::   aoa_y_nh4  !: Yield for AOA growth on NH4
   REAL(wp), PUBLIC ::   aoa_y_oxy  !: Yield for AOA growth on O2
   REAL(wp), PUBLIC ::   aoa_CN     !: Carbon:Nitrogen ratio for AOA
   REAL(wp), PUBLIC ::   nob_pocoef !: Rate of oxygen uptake by NOB via diffusion 
   REAL(wp), PUBLIC ::   nob_mumax  !: Maximum growth rate of NOB
   REAL(wp), PUBLIC ::   nob_y_no2  !: Yield for NOB growth on NO2
   REAL(wp), PUBLIC ::   nob_y_oxy  !: Yield for NOB growth on O2
   REAL(wp), PUBLIC ::   nob_CN     !: Carbon:Nitrogen ratio for NOB

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitraoa    !: ammonia oxidiser growth array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitrnob    !: nitrite oxidiser growth array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr       !: full denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitrno3    !: denitrification (NO3-->NO2) array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitrno2    !: denitrification (NO2-->N2) array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zaltrem      !: anaerobic remin without O2, NO3 or NO2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zanammox     !: anammox array (NH4+NO2-->N2)


   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zrem.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremik, zremikc, zremikn, zremikp, zsiremin, zfact 
      REAL(wp) ::   zsatur, zsatur2, znusil, znusil2, zdep, zdepmin, zfactdep
      REAL(wp) ::   zbactfer, zolimit, zonitr, zrfact2
      REAL(wp) ::   zammonic, zoxyremc, zoxyremn, zoxyremp, znitrate2ton
      REAL(wp) ::   zosil, ztem, zdenitnh4, zolimic, zolimin, zolimip, zdenitrn, zdenitrp
      REAL(wp) ::   zuptoxy, zuptnh4, zuptno2, znitr_prodno2
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj    ) :: ztempbac
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdepbac, zolimi, zdepprod, zfacsi, zfacsib, zdepeff, zfebact
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: mu_aoa, mu_nob
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_rem')
      !
      ! Initialisation of arrys
      zdepprod(:,:,:) = 1._wp
      zdepeff (:,:,:) = 0.3_wp
      ztempbac(:,:)   = 0._wp
      zfacsib(:,:,:)  = xsilab / ( 1.0 - xsilab )
      zfebact(:,:,:)  = 0._wp
      zfacsi(:,:,:)   = xsilab
      mu_aoa(:,:,:)   = 0.0
      mu_nob(:,:,:)   = 0.0

      ! Computation of the mean phytoplankton concentration as
      ! a crude estimate of the bacterial biomass
      ! this parameterization has been deduced from a model version
      ! that was modeling explicitely bacteria
      ! -------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep = MAX( hmld(ji,jj), heup(ji,jj) )
               IF( gdept_n(ji,jj,jk) < zdep ) THEN
                  zdepbac(ji,jj,jk) = MIN( 0.7 * ( trb(ji,jj,jk,jpzoo) + 2.* trb(ji,jj,jk,jpmes) ), 4.e-6 )
                  ztempbac(ji,jj)   = zdepbac(ji,jj,jk)
               ELSE
                  zdepmin = MIN( 1., zdep / gdept_n(ji,jj,jk) )
                  zdepbac (ji,jj,jk) = zdepmin**0.683 * ztempbac(ji,jj)
                  zdepprod(ji,jj,jk) = zdepmin**0.273
                  zdepeff (ji,jj,jk) = zdepeff(ji,jj,jk) * zdepmin**0.3
               ENDIF
            END DO
         END DO
      END DO

      IF( ln_p4z ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! DOC ammonification. Depends on depth, phytoplankton biomass
                  ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
                  zremik = xremik * xstep / 1.e-6 * xlimbac(ji,jj,jk) * zdepbac(ji,jj,jk) 
                  zremik = MAX( zremik, 2.74e-4 * xstep )
                  
                  ! Ammonification in oxic waters with oxygen consumption
                  ! -----------------------------------------------------
                  zolimit = zremik * ( 1.- nitrfac(ji,jj,jk) ) * trb(ji,jj,jk,jpdoc) 
                  zolimi(ji,jj,jk) = MIN( ( trb(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimit ) 

                  ! Ammonification in suboxic waters with denitrification
                  ! -------------------------------------------------------
                  zammonic = zremik * nitrfac(ji,jj,jk) * trb(ji,jj,jk,jpdoc)
                  denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )

                  ! Do two-step denitrification
                  ! -------------------------------------------------------
                  znitrate2ton = ( trb(ji,jj,jk,jpno3)+rtrn ) / ( trb(ji,jj,jk,jpno3)+trb(ji,jj,jk,jpno2)+rtrn )
                  denitrno3(ji,jj,jk) = denitr(ji,jj,jk) * znitrate2ton ! NO3 --> NO2 denitrification
                  denitrno2(ji,jj,jk) = denitr(ji,jj,jk) * (1.0 - znitrate2ton) * 0.1 ! NO2 --> N2 denitrification
                  denitrno3(ji,jj,jk)  = MIN((trb(ji,jj,jk,jpno3)-rtrn) / rdenit, denitrno3(ji,jj,jk)) ! only remove available NO3
                  denitrno2(ji,jj,jk)  = MIN((trb(ji,jj,jk,jpno2)-rtrn) / rdenit, denitrno2(ji,jj,jk)) ! only remove available NO2

                  ! anaerobic remineralisation without NO3, NO2 or O2
                  ! -------------------------------------------------------
                  zaltrem(ji,jj,jk) = zammonic - denitr(ji,jj,jk)

                  ! make sure all arrays are positive
                  ! -------------------------------------------------------
                  zolimi(ji,jj,jk) = MAX( 0.e0, zolimi(ji,jj,jk) )
                  denitrno3(ji,jj,jk) = MAX( 0.e0, denitrno3(ji,jj,jk) )
                  denitrno2(ji,jj,jk) = MAX( 0.e0, denitrno2(ji,jj,jk) )
                  zaltrem(ji,jj,jk) = MAX( 0.e0, zaltrem(ji,jj,jk) )

                  ! update the full denitrification array
                  ! -------------------------------------------------------
                  denitr(ji,jj,jk) = denitrno3(ji,jj,jk) + denitrno2(ji,jj,jk)

                  ! Update the tracer arrays
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jpno2) = tra(ji,jj,jk,jpno2) - denitrno2(ji,jj,jk) * rdenit + denitrno3(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - denitrno3(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zolimi(ji,jj,jk) - denitr(ji,jj,jk) - zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zolimi(ji,jj,jk) * o2ut
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) +                                                &
                  &                     rno3 * ( zolimi(ji,jj,jk) + zaltrem(ji,jj,jk + denitr(ji,jj,jk) )    &
                  &                     + rdenit * denitrno2(ji,jj,jk) )  ! removal of NO2-->N2 increases alkalinity
               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! DOC ammonification. Depends on depth, phytoplankton biomass
                  ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
                  ! -----------------------------------------------------------------
                  zremik = xstep / 1.e-6 * MAX(0.01, xlimbac(ji,jj,jk)) * zdepbac(ji,jj,jk) 
                  zremik = MAX( zremik, 2.74e-4 * xstep / xremikc )

                  zremikc = xremikc * zremik
                  zremikn = xremikn / xremikc
                  zremikp = xremikp / xremikc

                  ! Ammonification in oxic waters with oxygen consumption
                  ! -----------------------------------------------------
                  zolimit = zremikc * ( 1.- nitrfac(ji,jj,jk) ) * trb(ji,jj,jk,jpdoc) 
                  zolimic = MAX( 0.e0, MIN( ( trb(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimit ) ) 
                  zolimi(ji,jj,jk) = zolimic
                  zolimin = zremikn * zolimic * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zolimip = zremikp * zolimic * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn ) 

                  ! Ammonification in suboxic waters with denitrification
                  ! -------------------------------------------------------
                  zammonic = zremikc * nitrfac(ji,jj,jk) * trb(ji,jj,jk,jpdoc)
                  denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )
                  denitr(ji,jj,jk)  = MAX(0., MIN(  ( trb(ji,jj,jk,jpno3) - rtrn ) / rdenit, denitr(ji,jj,jk) ) )
                  zoxyremc          = MAX(0., zammonic - denitr(ji,jj,jk))
                  zdenitrn  = zremikn * denitr(ji,jj,jk) * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zdenitrp  = zremikp * denitr(ji,jj,jk) * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zoxyremn  = zremikn * zoxyremc * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zoxyremp  = zremikp * zoxyremc * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn )

                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zolimip + zdenitrp + zoxyremp
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zolimin + zdenitrn + zoxyremn
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - denitr(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zolimic - denitr(ji,jj,jk) - zoxyremc
                  tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) - zolimin - zdenitrn - zoxyremn
                  tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) - zolimip - zdenitrp - zoxyremp
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zolimic * o2ut
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zolimic + denitr(ji,jj,jk) + zoxyremc
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zolimin + zoxyremn + ( rdenit + 1.) * zdenitrn )
               END DO
            END DO
         END DO
         !
      ENDIF


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! NH4 nitrification to NO3. Ceased for oxygen concentrations
               ! below 2 umol/L. Inhibited at strong light 
               ! ----------------------------------------------------------
               !zonitrnh4(ji,jj,jk)  = nitrif * xstep * trb(ji,jj,jk,jpnh4) * (1.- nitrfac(ji,jj,jk) )  &
               !&                      / ( 1.+ emoy(ji,jj,jk) ) * ( 1. +fr_i(ji,jj) * emoy(ji,jj,jk) )
               !zonitrno2(ji,jj,jk)  = nitrif * xstep * trb(ji,jj,jk,jpno2) * (1.- nitrfac(ji,jj,jk) )  &
               !&                      / ( 1.+ emoy(ji,jj,jk) ) * ( 1. +fr_i(ji,jj) * emoy(ji,jj,jk) )

               !~~~~~~~~~~~~~~~~~~~~~~~!
               !!! Ammonia oxidation !!!
               !~~~~~~~~~~~~~~~~~~~~~~~!
                 ! Get uptake rates of oxygen and ammonium
               zuptoxy = aoa_pocoef * trb(ji,jj,jk,jpoxy)*1e6  ! mmol O2 / mmol C / day
               !print*, trb(ji,jj,jk,jpoxy)*1e6, trb(ji,jj,jk,jpnh4), xaoanh4(ji,jj,jk)
               zuptnh4 = ( aoa_mumax / aoa_y_nh4 ) * xaoanh4(ji,jj,jk)  ! mmol NH4 / mmol C / day
                 ! Calculate growth rate (Liebig Law of min)
               mu_aoa(ji,jj,jk) = MAX(0.0, MIN( zuptoxy * aoa_y_oxy, zuptnh4 * aoa_y_nh4 ))  ! /day
                 ! Calculate ammonia oxidiser biomass accumulation rate
               zonitraoa(ji,jj,jk) = mu_aoa(ji,jj,jk) * xstep * trb(ji,jj,jk,jpaoa)  ! mmol C / ts
               
               !~~~~~~~~~~~~~~~~~~~~~~~!
               !!! Nitrite oxidation !!!
               !~~~~~~~~~~~~~~~~~~~~~~~!
                 ! Get uptake rates of oxygen and nitrite
               zuptoxy = nob_pocoef * trb(ji,jj,jk,jpoxy)*1e6  ! mmol O2 / mmol C / day
               zuptno2 = ( nob_mumax / nob_y_no2 ) * xnobno2(ji,jj,jk)  ! mmol NO2 / mmol C / day
                 ! Calculate growth rate (Liebig Law of min)
               mu_nob(ji,jj,jk) = MAX(0.0, MIN( zuptoxy * nob_y_oxy, zuptno2 * nob_y_no2 ))  ! /day
                 ! Calculate nitrite oxidiser biomass accumulation rate
               zonitrnob(ji,jj,jk)  = mu_nob(ji,jj,jk) * xstep * trb(ji,jj,jk,jpnob)  ! mmol C / ts

               ! Loss of NH4 and NO2 due to anammox
               ! ----------------------------------------------------------
               zanammox(ji,jj,jk) = 0.01 * xstep * trb(ji,jj,jk,jpnh4) * nitrfac(ji,jj,jk)
               zanammox(ji,jj,jk) = max(0.0, min(  ( trb(ji,jj,jk,jpno2) - rtrn ) / 1.3, zanammox(ji,jj,jk) ) )
                        ! 1.3 mol of NO2 required per mol of NH4 oxidised by anammox (Brunner et al., 2013 PNAS)
 
               ! Update of the tracers trends (rno3 normalises NH4, NO2 and NO3 back to carbon units)
               ! ----------------------------
               tra(ji,jj,jk,jpaoa) = tra(ji,jj,jk,jpaoa) + zonitraoa(ji,jj,jk) 
               tra(ji,jj,jk,jpnob) = tra(ji,jj,jk,jpnob) + zonitrnob(ji,jj,jk)
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - (zonitraoa(ji,jj,jk)/aoa_y_nh4) / rno3      &
               &                     - (zonitrnob(ji,jj,jk)/nob_CN) / rno3                             &
               &                     - zanammox(ji,jj,jk)
               tra(ji,jj,jk,jpno2) = tra(ji,jj,jk,jpno2)                                               &
               &                     + ( zonitraoa(ji,jj,jk) * (1./aoa_y_nh4 - 1./aoa_CN) ) / rno3     &
               &                     - (zonitrnob(ji,jj,jk)/nob_y_no2) / rno3                          &  
               &                     - zanammox(ji,jj,jk)*1.3
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3)                                               &
               &                     + (zonitrnob(ji,jj,jk)/nob_y_no2) / rno3                          &
               &                     + zanammox(ji,jj,jk)*0.3
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy)                                               &
               &                     - zonitraoa(ji,jj,jk)/aoa_y_oxy                                   &
               &                     - zonitrnob(ji,jj,jk)/nob_y_oxy
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal)                                               &
               &                     - zonitraoa(ji,jj,jk)/aoa_y_nh4                                   &
               &                     - zonitraoa(ji,jj,jk) * (1./aoa_y_nh4 - 1./aoa_CN)
                        ! Wolf-Gladrow et al. (2007) 
                        ! Alkalinity increases by 1 mol for every 1 mol NO3/NO2 removed
                        ! Alkalinity decreases by 1 mol for every 1 mol NH4 removed
                        ! ammonia-oxidation removes two mol Alk because -1 NH4, +1 NO2 = -1 Alk, -1 Alk
                        ! Anammox has no net effect because -1 NH4, -1.3 NO2 and +0.3 NO3 = -1 Alk, +1.3 Alk, -0.3 Alk
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! Bacterial uptake of iron. No iron is available in DOC. So
               ! Bacteries are obliged to take up iron from the water. Some
               ! studies (especially at Papa) have shown this uptake to be significant
               ! ----------------------------------------------------------
               zbactfer = feratb *  rfact2 * 0.6_wp / rday * tgfunc(ji,jj,jk) * xlimbacl(ji,jj,jk)     &
                  &              * trb(ji,jj,jk,jpfer) / ( xkferb + trb(ji,jj,jk,jpfer) )    &
                  &              * zdepprod(ji,jj,jk) * zdepeff(ji,jj,jk) * zdepbac(ji,jj,jk)
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zbactfer*0.33
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zbactfer*0.25
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zbactfer*0.08
               zfebact(ji,jj,jk)   = zbactfer * 0.33
               blim(ji,jj,jk)      = xlimbacl(ji,jj,jk)  * zdepbac(ji,jj,jk) / 1.e-6 * zdepprod(ji,jj,jk)
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      ! Initialization of the array which contains the labile fraction
      ! of bSi. Set to a constant in the upper ocean
      ! ---------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep     = MAX( hmld(ji,jj), heup_01(ji,jj) )
               zsatur   = MAX( rtrn, ( sio3eq(ji,jj,jk) - trb(ji,jj,jk,jpsil) ) / ( sio3eq(ji,jj,jk) + rtrn ) )
               zsatur2  = ( 1. + tsn(ji,jj,jk,jp_tem) / 400.)**37
               znusil   = 0.225  * ( 1. + tsn(ji,jj,jk,jp_tem) / 15.) * zsatur + 0.775 * zsatur2 * zsatur**9.25
               ! Remineralization rate of BSi depedant on T and saturation
               ! ---------------------------------------------------------
               IF ( gdept_n(ji,jj,jk) > zdep ) THEN
                  zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk-1) * EXP( -0.5 * ( xsiremlab - xsirem )  &
                  &                   * znusil * e3t_n(ji,jj,jk) / wsbio4(ji,jj,jk) )
                  zfacsi(ji,jj,jk)  = zfacsib(ji,jj,jk) / ( 1.0 + zfacsib(ji,jj,jk) )
                  zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk) * EXP( -0.5 * ( xsiremlab - xsirem )    &
                  &                   * znusil * e3t_n(ji,jj,jk) / wsbio4(ji,jj,jk) )
               ENDIF
               zsiremin = ( xsiremlab * zfacsi(ji,jj,jk) + xsirem * ( 1. - zfacsi(ji,jj,jk) ) ) * xstep * znusil
               zosil    = zsiremin * trb(ji,jj,jk,jpgsi)
               !
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) - zosil
               tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) + zosil
            END DO
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      IF( knt == nrdttrc ) THEN
          zrfact2 = 1.e3 * rfact2r
          ALLOCATE( zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "REMIN" ) )  THEN
              zw3d(:,:,:) = zolimi(:,:,:) * tmask(:,:,:) * zfact !  Remineralisation rate
              CALL iom_put( "REMIN"  , zw3d )
          ENDIF
          IF( iom_use( "NITRNH4" ) )  THEN
              zw3d(:,:,:) = zonitraoa(:,:,:)/aoa_y_nh4 * tmask(:,:,:) * zfact ! 1st step of nitrification
              CALL iom_put( "NITRNH4"  , zw3d )
          ENDIF
          IF( iom_use( "NITRNO2" ) )  THEN
              zw3d(:,:,:) = zonitrnob(:,:,:)/nob_y_no2 * tmask(:,:,:) * zfact ! 2nd step of nitrification
              CALL iom_put( "NITRNO2"  , zw3d )
          ENDIF
          IF( iom_use( "DENIT" ) )  THEN
              zw3d(:,:,:) = denitr(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENIT"  , zw3d )
          ENDIF
          IF( iom_use( "DENITNO3" ) )  THEN
              zw3d(:,:,:) = denitrno3(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENITNO3"  , zw3d )
          ENDIF
          IF( iom_use( "DENITNO2" ) )  THEN
              zw3d(:,:,:) = denitrno2(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENITNO2"  , zw3d )
          ENDIF
          IF( iom_use( "ANAMMOX" ) )  THEN
              zw3d(:,:,:) = zanammox(:,:,:) * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "ANAMMOX"  , zw3d )
          ENDIF
          IF( iom_use( "ALTREM" ) )  THEN
              zw3d(:,:,:) = zaltrem(:,:,:) * tmask(:,:,:) * zfact ! anaerobic remin (DOC-->NH4) without O2,NO3,NO2
              CALL iom_put( "ALTREM"  , zw3d )
          ENDIF
          IF( iom_use( "BACT" ) )  THEN
               zw3d(:,:,:) = zdepbac(:,:,:) * 1.E6 * tmask(:,:,:)  ! Bacterial biomass
               CALL iom_put( "BACT", zw3d )
          ENDIF
          IF( iom_use( "FEBACT" ) )  THEN
               zw3d(:,:,:) = zfebact(:,:,:) * 1E9 * tmask(:,:,:) * zrfact2   ! Bacterial iron consumption
               CALL iom_put( "FEBACT" , zw3d )
          ENDIF
          !
          DEALLOCATE( zw3d )
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ xremik, nitrif, xsirem, xsiremlab, xsilab, feratb, xkferb, & 
         &                xremikc, xremikn, xremikp, aoa_pocoef, aoa_mumax,          &
         &                aoa_y_nh4, aoa_y_oxy, aoa_CN, nob_pocoef, nob_mumax,       &
         &                nob_y_no2, nob_y_oxy, nob_CN
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of remineralization parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisrem in reference namelist : Pisces remineralization
      READ  ( numnatp_ref, nampisrem, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisrem in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisrem in configuration namelist : Pisces remineralization
      READ  ( numnatp_cfg, nampisrem, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisrem in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist parameters for remineralization, nampisrem'
         IF( ln_p4z ) THEN
            WRITE(numout,*) '      remineralization rate of DOC              xremik    =', xremik
         ELSE
            WRITE(numout,*) '      remineralization rate of DOC              xremikc   =', xremikc
            WRITE(numout,*) '      remineralization rate of DON              xremikn   =', xremikn
            WRITE(numout,*) '      remineralization rate of DOP              xremikp   =', xremikp
         ENDIF
         WRITE(numout,*) '      remineralization rate of Si               xsirem    =', xsirem
         WRITE(numout,*) '      fast remineralization rate of Si          xsiremlab =', xsiremlab
         WRITE(numout,*) '      fraction of labile biogenic silica        xsilab    =', xsilab
         WRITE(numout,*) '      NH4 nitrification rate                    nitrif    =', nitrif
         WRITE(numout,*) '      Bacterial Fe/C ratio                      feratb    =', feratb
         WRITE(numout,*) '      Half-saturation constant for bact. Fe/C   xkferb    =', xkferb
         WRITE(numout,*) '      AOA diffusive uptake limit (m3/mmolC/day) aoa_pocoef =', aoa_pocoef
         WRITE(numout,*) '      AOA maximum growth rate (/day)            aoa_mumax =', aoa_mumax
         WRITE(numout,*) '      AOA biomass yield on NH4 (mol C / mol NH4) aoa_y_nh4 =', aoa_y_nh4
         WRITE(numout,*) '      AOA biomass yield on O2 (mol C / mol O2)  aoa_y_oxy =', aoa_y_oxy
         WRITE(numout,*) '      AOA Carbon:Nitrogen ratio of biomass      aoa_CN    =', aoa_CN
         WRITE(numout,*) '      NOB diffusive uptake limit (m3/mmolC/day) nob_pocoef =', nob_pocoef
         WRITE(numout,*) '      NOB maximum growth rate (/day)            nob_mumax =', nob_mumax
         WRITE(numout,*) '      NOB biomass yield on NO2 (mol C / mol NO2) nob_y_no2 =', nob_y_no2
         WRITE(numout,*) '      NOB biomass yield on O2 (mol C / mol O2)  nob_y_oxy =', nob_y_oxy
         WRITE(numout,*) '      NOB Carbon:Nitrogen ratio of biomass      nob_CN    =', nob_CN
      ENDIF
      !
      zonitraoa(:,:,:) = 0._wp
      zonitrnob(:,:,:) = 0._wp
      denitr(:,:,:) = 0._wp
      denitrno3(:,:,:) = 0._wp
      denitrno2(:,:,:) = 0._wp
      zaltrem(:,:,:) = 0._wp
      zanammox(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( zonitraoa(jpi,jpj,jpk), zonitrnob(jpi,jpj,jpk),denitr(jpi,jpj,jpk),   &
      &         denitrno3(jpi,jpj,jpk), denitrno2(jpi,jpj,jpk),zaltrem(jpi,jpj,jpk),  &
      &         zanammox(jpi,jpj,jpk), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_rem_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_rem_alloc

   !!======================================================================
END MODULE p4zrem
