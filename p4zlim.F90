MODULE p4zlim
   !!======================================================================
   !!                         ***  MODULE p4zlim  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
   !!----------------------------------------------------------------------
   !!   p4z_lim        :   Compute the nutrients limitation terms 
   !!   p4z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE oce_trc         ! Shared ocean-passive tracers variables
   USE trc             ! Tracers defined
   USE sms_pisces      ! PISCES variables
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC p4z_lim    
   PUBLIC p4z_lim_init    
   PUBLIC p4z_lim_alloc

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concnno3    !:  NO3, PO4 half saturation   
   REAL(wp), PUBLIC ::  concdno3    !:  Phosphate half saturation for diatoms  
   REAL(wp), PUBLIC ::  concnnh4    !:  NH4 half saturation for phyto  
   REAL(wp), PUBLIC ::  concdnh4    !:  NH4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concnfer    !:  Iron half saturation for nanophyto 
   REAL(wp), PUBLIC ::  concdfer    !:  Iron half saturation for diatoms  
   REAL(wp), PUBLIC ::  concbno3    !:  NO3 half saturation  for bacteria 
   REAL(wp), PUBLIC ::  concbnh4    !:  NH4 half saturation for bacteria
   REAL(wp), PUBLIC ::  xsizedia    !:  Minimum size criteria for diatoms
   REAL(wp), PUBLIC ::  xsizephy    !:  Minimum size criteria for nanophyto
   REAL(wp), PUBLIC ::  xsizern     !:  Size ratio for nanophytoplankton
   REAL(wp), PUBLIC ::  xsizerd     !:  Size ratio for diatoms
   REAL(wp), PUBLIC ::  xksi1       !:  half saturation constant for Si uptake 
   REAL(wp), PUBLIC ::  xksi2       !:  half saturation constant for Si/C 
   REAL(wp), PUBLIC ::  xkdoc       !:  2nd half-sat. of DOC remineralization  
   REAL(wp), PUBLIC ::  concbfe     !:  Fe half saturation for bacteria 
   REAL(wp), PUBLIC ::  oxymin      !:  half saturation constant for anoxia
   REAL(wp), PUBLIC ::  qnfelim     !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qdfelim     !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  caco3r      !:  mean rainratio 
   REAL(wp), PUBLIC ::  aoa_k_nh4   !:  half saturation constant for AOA NH4 uptake 
   REAL(wp), PUBLIC ::  nob_k_no2   !:  half saturation constant for NOB NO2 uptake 
   REAL(wp), PUBLIC ::  aox_k_nh4   !:  half saturation constant for AOX NH4 uptake 
   REAL(wp), PUBLIC ::  aox_k_no2   !:  half saturation constant for AOX NO2 uptake 
   REAL(wp), PUBLIC ::  nar_k_doc   !:  half saturation constant for NAR DOC uptake 
   REAL(wp), PUBLIC ::  nir_k_doc   !:  half saturation constant for NIR DOC uptake 
   REAL(wp), PUBLIC ::  nar_k_no3   !:  half saturation constant for NAR NO3 uptake 
   REAL(wp), PUBLIC ::  nir_k_no2   !:  half saturation constant for NIR NO2 uptake 

   !!* Phytoplankton limitation terms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanono3   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatno3   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanonh4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatnh4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanopo4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatpo4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphy    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdia    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimsi     !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimbac    !: ??
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimbacl   !: ??
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concnfe    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xaoanh4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnobno2   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xaoxnh4   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xaoxno2   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnardoc   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnirdoc   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnarno3   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnirno2   !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xhetfer   !: ???

   ! Coefficient for iron limitation
   REAL(wp) ::  xcoef1   = 0.0016  / 55.85  
   REAL(wp) ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.625 * 0.5 * 1.5
   REAL(wp) ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.625 * 0.5 

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zlim.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_lim( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!              for the various phytoplankton species
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zno3, zferlim, zton
      REAL(wp) ::   zconcd, zconcd2, zconcn, zconcn2
      REAL(wp) ::   z1_trbdia, z1_trbphy, ztem1, ztem2, zetot1, zetot2
      REAL(wp) ::   zdenom, zratio, zironmin
      REAL(wp) ::   zconc1d, zconc1dnh4, zconc0n, zconc0nnh4   
      REAL(wp) ::   zlimnh4, zlimno3, znutlimtot, zbactnh4, zbactno3
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_lim')
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               
               ! Tuning of the iron concentration to a minimum level that is set to the detection limit
               !-------------------------------------
               zno3    = trb(ji,jj,jk,jpno3) / 40.e-6
               zferlim = MAX( 3e-11 * zno3 * zno3, 5e-12 )
               zferlim = MIN( zferlim, 7e-11 )
               trb(ji,jj,jk,jpfer) = MAX( trb(ji,jj,jk,jpfer), zferlim )

               ! Computation of a variable Ks for iron on diatoms taking into account
               ! that increasing biomass is made of generally bigger cells
               !------------------------------------------------
               zconcd   = MAX( 0.e0 , trb(ji,jj,jk,jpdia) - xsizedia )
               zconcd2  = trb(ji,jj,jk,jpdia) - zconcd
               zconcn   = MAX( 0.e0 , trb(ji,jj,jk,jpphy) - xsizephy )
               zconcn2  = trb(ji,jj,jk,jpphy) - zconcn
               z1_trbphy   = 1. / ( trb(ji,jj,jk,jpphy) + rtrn )
               z1_trbdia   = 1. / ( trb(ji,jj,jk,jpdia) + rtrn )

               concdfe(ji,jj,jk) = MAX( concdfer, ( zconcd2 * concdfer + concdfer * xsizerd * zconcd ) * z1_trbdia )
               zconc1d           = MAX( concdno3, ( zconcd2 * concdno3 + concdno3 * xsizerd * zconcd ) * z1_trbdia )
               zconc1dnh4        = MAX( concdnh4, ( zconcd2 * concdnh4 + concdnh4 * xsizerd * zconcd ) * z1_trbdia )

               concnfe(ji,jj,jk) = MAX( concnfer, ( zconcn2 * concnfer + concnfer * xsizern * zconcn ) * z1_trbphy )
               zconc0n           = MAX( concnno3, ( zconcn2 * concnno3 + concnno3 * xsizern * zconcn ) * z1_trbphy )
               zconc0nnh4        = MAX( concnnh4, ( zconcn2 * concnnh4 + concnnh4 * xsizern * zconcn ) * z1_trbphy )

               ! Total oxidised inorganic nitrogen
               zton = trb(ji,jj,jk,jpno2) + trb(ji,jj,jk,jpno3)

               ! Michaelis-Menten Limitation term for nutrients Small bacteria
               ! -------------------------------------------------------------
               zlimnh4 = trb(ji,jj,jk,jpnh4) / ( concbno3 + trb(ji,jj,jk,jpnh4))
               zlimno3 = zton / ( concbno3 + zton)
               znutlimtot = ( trb(ji,jj,jk,jpnh4) + zton ) / (concbno3 + trb(ji,jj,jk,jpnh4) + zton )

               zbactnh4 = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
               zbactno3 = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
               !
               zlim1    = zbactno3 + zbactnh4
               zlim2    = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + concbnh4 )
               zlim3    = trb(ji,jj,jk,jpfer) / ( concbfe + trb(ji,jj,jk,jpfer) )
               zlim4    = trb(ji,jj,jk,jpdoc) / ( xkdoc   + trb(ji,jj,jk,jpdoc) )
               xlimbacl(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               xlimbac (ji,jj,jk) = MIN( zlim1, zlim2, zlim3 ) * zlim4

               ! Michaelis-Menten Limitation term for nutrients Small flagellates
               ! -----------------------------------------------
               zlimnh4 = trb(ji,jj,jk,jpnh4) / ( zconc0n + trb(ji,jj,jk,jpnh4) )
               zlimno3 = zton / ( zconc0n + zton )
               znutlimtot = ( trb(ji,jj,jk,jpnh4) + zton ) / ( zconc0n + trb(ji,jj,jk,jpnh4) + zton )

               xnanonh4(ji,jj,jk) = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
               xnanono3(ji,jj,jk) = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
               !
               zlim1    = xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk)
               zlim2    = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + zconc0nnh4 )
               zratio   = trb(ji,jj,jk,jpnfe) * z1_trbphy
               zironmin = xcoef1 * trb(ji,jj,jk,jpnch) * z1_trbphy + xcoef2 * zlim1 + xcoef3 * xnanono3(ji,jj,jk)
               zlim3    = MAX( 0.,( zratio - zironmin ) / qnfelim )
               xnanopo4(ji,jj,jk) = zlim2
               xlimnfe (ji,jj,jk) = MIN( 1., zlim3 )
               xlimphy (ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               !
               !   Michaelis-Menten Limitation term for nutrients Diatoms
               !   ----------------------------------------------
               zlimnh4 = trb(ji,jj,jk,jpnh4) / ( zconc1d + trb(ji,jj,jk,jpnh4) )
               zlimno3 = zton / ( zconc1d + zton )
               znutlimtot = ( trb(ji,jj,jk,jpnh4) + zton ) / ( zconc1d + trb(ji,jj,jk,jpnh4) + zton )

               xdiatnh4(ji,jj,jk) = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
               xdiatno3(ji,jj,jk) = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + rtrn )
               !
               zlim1    = xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk)
               zlim2    = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + zconc1dnh4  )
               zlim3    = trb(ji,jj,jk,jpsil) / ( trb(ji,jj,jk,jpsil) + xksi(ji,jj) )
               zratio   = trb(ji,jj,jk,jpdfe) * z1_trbdia
               zironmin = xcoef1 * trb(ji,jj,jk,jpdch) * z1_trbdia + xcoef2 * zlim1 + xcoef3 * xdiatno3(ji,jj,jk)
               zlim4    = MAX( 0., ( zratio - zironmin ) / qdfelim )
               xdiatpo4(ji,jj,jk) = zlim2
               xlimdfe (ji,jj,jk) = MIN( 1., zlim4 )
               xlimdia (ji,jj,jk) = MIN( zlim1, zlim2, zlim3, zlim4 )
               xlimsi  (ji,jj,jk) = MIN( zlim1, zlim2, zlim4 )
               !
               !   Limitation terms for chemoautotrophs
               !   ----------------------------------------------
               xaoanh4(ji,jj,jk) = ( (trb(ji,jj,jk,jpnh4)+rtrn) / ( trb(ji,jj,jk,jpnh4) + aoa_k_nh4/rno3 + rtrn ) )
               xnobno2(ji,jj,jk) = ( (trb(ji,jj,jk,jpno2)+rtrn) / ( trb(ji,jj,jk,jpno2) + nob_k_no2/rno3 + rtrn ) )
               xaoxnh4(ji,jj,jk) = ( (trb(ji,jj,jk,jpnh4)+rtrn) / ( trb(ji,jj,jk,jpnh4) + aox_k_nh4/rno3 + rtrn ) )
               xaoxno2(ji,jj,jk) = ( (trb(ji,jj,jk,jpno2)+rtrn) / ( trb(ji,jj,jk,jpno2) + aox_k_no2/rno3 + rtrn ) )
               !
               !   Limitation terms for heterotrophs
               !   ----------------------------------------------
               xnardoc(ji,jj,jk) = ( (trb(ji,jj,jk,jpdoc)+rtrn) / ( trb(ji,jj,jk,jpdoc) + nar_k_doc + rtrn ) )
               xnirdoc(ji,jj,jk) = ( (trb(ji,jj,jk,jpdoc)+rtrn) / ( trb(ji,jj,jk,jpdoc) + nir_k_doc + rtrn ) )
               xnarno3(ji,jj,jk) = ( (trb(ji,jj,jk,jpno3)+rtrn) / ( trb(ji,jj,jk,jpno3) + nar_k_no3/rno3 + rtrn ) )
               xnirno2(ji,jj,jk) = ( (trb(ji,jj,jk,jpno2)+rtrn) / ( trb(ji,jj,jk,jpno2) + nir_k_no2/rno3 + rtrn ) )
               !print*, jk, xnirno2(ji,jj,jk), trb(ji,jj,jk,jpno2)
               !print*, "DOC", xnirdoc(ji,jj,jk), trb(ji,jj,jk,jpdoc)
               xhetfer(ji,jj,jk) = trb(ji,jj,jk,jpfer) / ( concbfe + trb(ji,jj,jk,jpfer) )
               !
           END DO
         END DO
      END DO

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! --------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zton = trb(ji,jj,jk,jpno2) + trb(ji,jj,jk,jpno3) ! total oxidised nitrogen
               zlim1 =  ( zton * concnnh4 + trb(ji,jj,jk,jpnh4) * concnno3 )    &
                  &   / ( concnno3 * concnnh4 + concnnh4 * zton + concnno3 * trb(ji,jj,jk,jpnh4) )
               zlim2  = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + concnnh4 )
               zlim3  = trb(ji,jj,jk,jpfer) / ( trb(ji,jj,jk,jpfer) +  5.E-11   )
               ztem1  = MAX( 0., tsn(ji,jj,jk,jp_tem) )
               ztem2  = tsn(ji,jj,jk,jp_tem) - 10.
               zetot1 = MAX( 0., etot_ndcy(ji,jj,jk) - 1.) / ( 4. + etot_ndcy(ji,jj,jk) ) 
               zetot2 = 30. / ( 30. + etot_ndcy(ji,jj,jk) ) 

               xfracal(ji,jj,jk) = caco3r * MIN( zlim1, zlim2, zlim3 )                  &
                  &                       * ztem1 / ( 0.1 + ztem1 )                     &
                  &                       * MAX( 1., trb(ji,jj,jk,jpphy) * 1.e6 / 2. )  &
                  &                       * zetot1 * zetot2               &
                  &                       * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )         &
                  &                       * MIN( 1., 50. / ( hmld(ji,jj) + rtrn ) )
               xfracal(ji,jj,jk) = MIN( 0.8 , xfracal(ji,jj,jk) )
               xfracal(ji,jj,jk) = MAX( 0.02, xfracal(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! denitrification factor computed from O2 levels
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.e-6  - trb(ji,jj,jk,jpoxy) )    &
                  &                                / ( oxymin + trb(ji,jj,jk,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
               !
               ! denitrification factor computed from NO3 levels
               zton = trb(ji,jj,jk,jpno2) + trb(ji,jj,jk,jpno3) ! total oxidised nitrogen
               nitrfac2(ji,jj,jk) = MAX( 0.e0,       ( 1.E-6 - zton )  &
                  &                                / ( 1.E-6 + zton ) )
               nitrfac2(ji,jj,jk) = MIN( 1., nitrfac2(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN        ! save output diagnostics
        IF( iom_use( "xfracal" ) )   CALL iom_put( "xfracal", xfracal(:,:,:) * tmask(:,:,:) )  ! euphotic layer deptht
        IF( iom_use( "LNnut"   ) )   CALL iom_put( "LNnut"  , xlimphy(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LDnut"   ) )   CALL iom_put( "LDnut"  , xlimdia(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LNFe"    ) )   CALL iom_put( "LNFe"   , xlimnfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LDFe"    ) )   CALL iom_put( "LDFe"   , xlimdfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_lim')
      !
   END SUBROUTINE p4z_lim


   SUBROUTINE p4z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of nutrient limitation parameters
      !!
      !! ** Method  :   Read the nampislim namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampislim
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zlim/ concnno3, concdno3, concnnh4, concdnh4, concnfer, concdfer, concbfe,   &
         &                concbno3, concbnh4, xsizedia, xsizephy, xsizern, xsizerd,          & 
         &                xksi1, xksi2, xkdoc, qnfelim, qdfelim, caco3r, oxymin,             &
         &                aoa_k_nh4, nob_k_no2, aox_k_nh4, aox_k_no2,                        & 
         &                nar_k_doc, nir_k_doc, nar_k_no3, nir_k_no2
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_lim_init : initialization of nutrient limitations'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampislim in reference namelist : Pisces nutrient limitation parameters
      READ  ( numnatp_ref, namp4zlim, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zlim in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampislim in configuration namelist : Pisces nutrient limitation parameters 
      READ  ( numnatp_cfg, namp4zlim, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zlim in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zlim )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zlim'
         WRITE(numout,*) '      mean rainratio                           caco3r    = ', caco3r
         WRITE(numout,*) '      NO3 half saturation of nanophyto         concnno3  = ', concnno3
         WRITE(numout,*) '      NO3 half saturation of diatoms           concdno3  = ', concdno3
         WRITE(numout,*) '      NH4 half saturation for phyto            concnnh4  = ', concnnh4
         WRITE(numout,*) '      NH4 half saturation for diatoms          concdnh4  = ', concdnh4
         WRITE(numout,*) '      half saturation constant for Si uptake   xksi1     = ', xksi1
         WRITE(numout,*) '      half saturation constant for Si/C        xksi2     = ', xksi2
         WRITE(numout,*) '      half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(numout,*) '      Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(numout,*) '      Iron half saturation for diatoms         concdfer  = ', concdfer
         WRITE(numout,*) '      size ratio for nanophytoplankton         xsizern   = ', xsizern
         WRITE(numout,*) '      size ratio for diatoms                   xsizerd   = ', xsizerd
         WRITE(numout,*) '      NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(numout,*) '      NH4 half saturation for bacteria         concbnh4  = ', concbnh4
         WRITE(numout,*) '      Minimum size criteria for diatoms        xsizedia  = ', xsizedia
         WRITE(numout,*) '      Minimum size criteria for nanophyto      xsizephy  = ', xsizephy
         WRITE(numout,*) '      Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(numout,*) '      halk saturation constant for anoxia       oxymin   =' , oxymin
         WRITE(numout,*) '      optimal Fe quota for nano.               qnfelim   = ', qnfelim
         WRITE(numout,*) '      Optimal Fe quota for diatoms             qdfelim   = ', qdfelim
         WRITE(numout,*) '      NH4 half saturation for AOA              aoa_k_nh4 = ', aoa_k_nh4
         WRITE(numout,*) '      NO2 half saturation for NOB              nob_k_no2 = ', nob_k_no2
         WRITE(numout,*) '      NH4 half saturation for AOX              aox_k_nh4 = ', aox_k_nh4
         WRITE(numout,*) '      NO2 half saturation for AOX              aox_k_no2 = ', aox_k_no2
         WRITE(numout,*) '      DOC half saturation for NAR              nar_k_doc = ', nar_k_doc
         WRITE(numout,*) '      DOC half saturation for NIR              nir_k_doc = ', nir_k_doc
         WRITE(numout,*) '      NO3 half saturation for NAR              nar_k_no3 = ', nar_k_no3
         WRITE(numout,*) '      NO2 half saturation for NIR              nir_k_no2 = ', nir_k_no2
      ENDIF
      !
      nitrfac (:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_lim_init


   INTEGER FUNCTION p4z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim_alloc  ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      !!----------------------------------------------------------------------

      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xnanono3(jpi,jpj,jpk), xdiatno3(jpi,jpj,jpk),       &
         &      xnanonh4(jpi,jpj,jpk), xdiatnh4(jpi,jpj,jpk),       &
         &      xnanopo4(jpi,jpj,jpk), xdiatpo4(jpi,jpj,jpk),       &
         &      xaoanh4(jpi,jpj,jpk),  xnobno2(jpi,jpj,jpk),        &
         &      xaoxnh4(jpi,jpj,jpk),  xaoxno2(jpi,jpj,jpk),        &
         &      xnardoc(jpi,jpj,jpk),  xnirdoc(jpi,jpj,jpk),        &
         &      xnarno3(jpi,jpj,jpk),  xnirno2(jpi,jpj,jpk),        &
         &      xlimphy (jpi,jpj,jpk), xlimdia (jpi,jpj,jpk),       &
         &      xlimnfe (jpi,jpj,jpk), xlimdfe (jpi,jpj,jpk),       &
         &      xlimbac (jpi,jpj,jpk), xlimbacl(jpi,jpj,jpk),       &
         &      concnfe (jpi,jpj,jpk), concdfe (jpi,jpj,jpk),       &
         &      xlimsi  (jpi,jpj,jpk), xhetfer(jpi,jpj,jpk), STAT=p4z_lim_alloc )
      !
      IF( p4z_lim_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_lim_alloc

   !!======================================================================
END MODULE p4zlim
