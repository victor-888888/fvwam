MODULE WAM_SOURCE_MODULE
    use fre_dir_mod
    use mesh_mod
    use mask_mod
    use const_mod
    use wave_propag_mod
    use initial_mod


IMPLICIT NONE

!    2. PARAMETERS FOR COUPLING.

REAL(real_kind), PARAMETER :: RNUAIR  = 1.5E-5     !! KINEMATIC AIR VISCOSITY
REAL(real_kind), PARAMETER :: RNUAIRM = 0.11*RNUAIR !! KINEMATIC AIR VISCOSITY FOR MOMENTUM TRANSFER
REAL(real_kind), PARAMETER :: XEPS = rho_air/rho_sea
REAL(real_kind), PARAMETER :: XINVEPS = 1.0_rk/XEPS

REAL(real_kind)            :: BETAMAX = 1.20       !! PARAMETER FOR WIND INPUT (ECMWF CY45R1).
REAL(real_kind)            :: ZALP    = 0.0080     !! SHIFTS GROWTH CURVE (ECMWF CY45R1).
REAL(real_kind)            :: ALPHA   = 0.0060     !! MINIMUM CHARNOCK CONSTANT (ECMWF CY45R1).
                                                   !! if LE 30 frequencies changed
                                                   !! to 0.0075 in subroutine INITMDL
REAL(real_kind)            :: TAUWSHELTER= 0.0     !! SHELTERING COEFFICIENT  (ECMWF CY45R1)
REAL(real_kind), PARAMETER :: XKAPPA  = 0.40       !! VON KARMAN CONSTANT.
REAL(real_kind), PARAMETER :: XNLEV   = 10.0       !! WINDSPEED REF. LEVEL.

!    4. PARAMETERS FOR SDISSIP_ARD ARDHUIN et al. 2010
REAL(real_kind), PARAMETER :: SDSBR = 9.0E-4
!!     Saturation dissipation coefficient
INTEGER, PARAMETER :: ISDSDTH = 80
INTEGER, PARAMETER :: ISB=2
INTEGER, PARAMETER :: IPSAT=2
REAL(real_kind), PARAMETER :: SSDSC2 = -2.2E-5
REAL(real_kind), PARAMETER :: SSDSC4 = 1.0
REAL(real_kind), PARAMETER :: SSDSC6 = 0.3
REAL(real_kind), PARAMETER :: MICHE = 1.0
REAL(real_kind), PARAMETER :: SSDSC3 = 0.0
REAL(real_kind), PARAMETER :: SSDSBRF1 = 0.5
REAL(real_kind), PARAMETER :: BRKPBCOEF=28.16
REAL(real_kind), PARAMETER :: SSDSC5  = 0.0

REAL(real_kind), PARAMETER :: DKMAX = 40.        !! MAXIMUM VALUE OF DEPTH*WAVENUMBER.

! fre-dir mod
integer:: KL 
integer:: ML
real(real_kind) :: CO = 1.1
REAL(real_kind),    SAVE      :: INV_LOG_CO     !! 1./LOG10(FREQUENCY RATIO).
REAL(real_kind),    SAVE      :: FMIN            !! MINIMUM ENERGY DENSITY = HS = 7CM

REAL(real_kind)      :: EMIN = 1.0E-12    !! REPLACES THE INTRINSIC TINY

INTEGER, PARAMETER :: JUMAX = 300   !! TABLE DIMENSION FOR U10.
REAL(real_kind),    PARAMETER :: UMAX  = 75.   !! MAXIMUM WIND SPEED IN TABLE.
REAL(real_kind)               :: DELU  = 0.25   !! WIND INCREMENT.
REAL(real_kind),    PARAMETER :: EPS1 = 0.00001 !! SMALL NUMBER TO MAKE SURE THAT A
                                      !! SOLUTION IS OBTAINED IN ITERATION
                                     !! WITH TAU>TAUW.
REAL(real_kind), PARAMETER :: FLMIN = 0.000001  !! ABSOLUTE MINIMUM ENERGY IN SPECTRAL BINS
REAL(real_kind),   ALLOCATABLE, DIMENSION(:,:) :: T_TAIL  !! TABLE FOR K**(-3)/
                                              !GROUP VELOCITY.
REAL(real_kind), ALLOCATABLE, DIMENSION(:,:) ::  FLMINFR !! THE MINIMUM VALUE IN SPECTRAL

REAL(real_kind), PARAMETER :: GAMD  = 0.8  !! Parameter of depth limited wave height


!! flux
REAL(real_kind),    ALLOCATABLE :: PHIOC(:)    !! ENERGY FLUX TO OCEAN.
REAL(real_kind),    ALLOCATABLE :: PHIAW(:)    !! ENERGY FLUX FROM WIND TO WAVES.
REAL(real_kind),    ALLOCATABLE :: TAUOC_X(:)  !! MOMENTUM FLUX INTO OCEAN.
REAL(real_kind),    ALLOCATABLE :: TAUOC_Y(:)  !! MOMENTUM FLUX INTO OCEAN.
REAL(real_kind),    ALLOCATABLE :: PHIBOT(:)   !! BOTTOM ENERGY FLUX TO OCEAN.
REAL(real_kind),    ALLOCATABLE :: TAUBOT_X(:) !! BOTTOM MOMENTUM FLUX INTO OCEAN.
REAL(real_kind),    ALLOCATABLE :: TAUBOT_Y(:) !! BOTTOM MOMENTUM FLUX INTO OCEAN.


! logical

LOGICAL :: shallow_run=.true., wave_breaking_run=.true.
INTEGER :: ISNONLIN = 0         !! = 0 OLD DEPTH SCALING FOR SNL,
                             !! ≠ 0 NEW DEPTH SCALING FOR SNL.

! ---------------------------------------------------------------------------- !
!
!    1. INDICES AND WEIGHTS USED TO COMPUTE THE NONLINEAR TRANSFER RATE.
!       ----------------------------------------------------------------

INTEGER, PARAMETER :: NINL = 5  !! SIZE OF INLCOEF
INTEGER, PARAMETER :: NRNL = 25 !! SIZE OF RNLCOEF

INTEGER :: KFRH     !! SIZE OF FRH
INTEGER :: MFRSTLW  !! INDEX OF FIRST EXTRA LOW FREQUENCY FOR SNL
INTEGER :: MLSTHG   !! INDEX OF LAST EXTRA HIGH FREQUENCY FOR SNL

INTEGER,ALLOCATABLE, DIMENSION(:)   :: IKP    !! FREQUENCY INDEX FOR STORING
                                              !! ENERGY TRANSFER INCREMENTS
                                              !! INTO BINS, WAVE NO. 3.
INTEGER,ALLOCATABLE, DIMENSION(:)   :: IKP1   !! IKP+1. 
INTEGER,ALLOCATABLE, DIMENSION(:)   :: IKM    !! FREQUENCY INDEX FOR STORING
                                              !! ENERGY TRANSFER INCREMENTS
                                              !! INTO BINS, WAVE NO. 4. 
INTEGER,ALLOCATABLE, DIMENSION(:)   :: IKM1   !! IKM+1
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: K1W    !! ANGULAR INDEX FOR STORING
                                              !! ENERGY TRANSFER INCREMENTS
                                              !! INTO BINS, WAVE NO. 3.
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: K2W    !! ANGULAR INDEX FOR STORING
                                              !! ENERGY TRANSFER INCREMENTS
                                              !! INTO BINS, WAVE NO. 4.
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: K11W   !! K1W(.,1)-1, K1W(.,2)+1.
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: K21W   !! K2W(.,1)+1, K2W(.,2)-1.
INTEGER,ALLOCATABLE, DIMENSION(:,:) :: INLCOEF  !! STORES ALL FREQUENCY DEPENDENT
                                                !! INDICES FOUND IN SNL

REAL(real_kind),   ALLOCATABLE, DIMENSION(:)   :: AF11   !! WEIGHTS FOR APPROXIMATION
                                              !! OF NONL. TRANSFER (ONE 
                                              !! TERM ONLY SET TO 3000). 
                                              !! MULTIPLIED BY FREQ. **11.
REAL(real_kind),   ALLOCATABLE, DIMENSION(:)   :: FKLAP  !! WEIGHT IN FREQUENCY GRID  
                                              !! FOR INTERPOLATION, WAVE  
                                              !! NO. 3 ("1+LAMBDA" TERM)
REAL(real_kind),   ALLOCATABLE, DIMENSION(:)   :: FKLAP1 !! 1-FKLAP.
REAL(real_kind),   ALLOCATABLE, DIMENSION(:)   :: FKLAM  !! WEIGHT IN FREQUENCY GRID  
                                             !! FOR INTERPOLATION, WAVE 
                                             !! NO. 4 ("1-LAMBDA" TERM).
REAL(real_kind),   ALLOCATABLE, DIMENSION(:)   :: FKLAM1 !! 

REAL(real_kind), ALLOCATABLE, DIMENSION(:)     :: FRH    !! TAIL FREQUENCY RATION **5
REAL(real_kind), ALLOCATABLE, DIMENSION(:)     :: FTRF   !! FRONT TAIL REDUCTION FACTOR 
                                             !! USED TO A SPECTRAL TAIL IN FRONT
                                             !! OF THE FIRST DISCRETISED FREQUENCY
REAL(real_kind), ALLOCATABLE, DIMENSION(:,:)   :: RNLCOEF !! STORES ALL FREQUENCY DEPENDENT
                                              !! COEFFICIENT FOUND IN SNL

REAL(real_kind)    :: ACL1        !! WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
                         !! WAVE NO. 3 ("1+LAMBDA" TERM).
REAL(real_kind)    :: ACL2        !! WEIGHT IN ANGULAR GRID FOR INTERPOLATION,
                         !! WAVE NO. 4 ("1-LAMBDA" TERM).
REAL(real_kind)    :: CL11        !! 1.-ACL1.
REAL(real_kind)    :: CL21        !! 1.-ACL2.
REAL(real_kind)    :: DAL1        !! 1./ACL1.
REAL(real_kind)    :: DAL2        !! 1./ACL2.

! ---------------------------------------------------------------------------- !
!                                                                              !
!    2. NONLINEAR TRANSFER FUNCTION COEFFICIENTS FOR SHALLOW WATER.            !
!       -----------------------------------------------------------            !

REAL(real_kind),   ALLOCATABLE, DIMENSION(:,:) :: ENH

! ---------------------------------------------------------------------------- !
!                                                                              !
!    3. INTEGRATION WEIGHT FOR HIGH FREQ STRESS AND ENERGY (SUB.TAU_PHI_HF).   !
!       --------------------------------------------------------------         !

INTEGER, PARAMETER :: JTOT_TAUHF=19  !! DIMENSION OF WTAUHF. IT MUST BE ODD !!!
REAL(real_kind)               :: WTAUHF(JTOT_TAUHF) !! INTEGRATION WEIGHT FOR TAU_PHI_HF
REAL(real_kind)               :: X0TAUHF        !! LOWEST LIMIT FOR INTEGRATION IN TAU_PHI_HF: X0 *(G/USTAR)

! ---------------------------------------------------------------------------- !
!
!    4. TABLE FOR SINPUT_ARD AND SDISSIP_ARD.
!       -------------------------------------

! FOR SINPUT_ARD
INTEGER(real_kind), PARAMETER :: IAB=200       !! DIMENSION OF SWELLFT.
REAL(real_kind)               :: SWELLFT(IAB)  !! FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS.

! FOR SDISSIP_ARD
INTEGER :: NSDSNTH  !! NUMBER OF DIRECTIONS TO COMPUTE THE SPECTRAL SATURATION
INTEGER :: NDIKCUMUL !! INTEGER DIFFERENCE IN FREQUENCY BANDS.
INTEGER, ALLOCATABLE :: INDICESSAT(:,:)
REAL(real_kind), ALLOCATABLE :: SATWEIGHTS(:,:)
REAL(real_kind), ALLOCATABLE :: CUMULW(:,:,:,:)

REAL(real_kind), ALLOCATABLE :: FL3(:,:,:)
REAL(real_kind), ALLOCATABLE :: USTAR(:)
REAL(real_kind), ALLOCATABLE :: Z0(:)
REAL(real_kind), ALLOCATABLE :: TAUW(:)
REAL(real_kind), ALLOCATABLE :: WSTAR(:)
REAL(real_kind), ALLOCATABLE :: ROAIRN(:)

CONTAINS


SUBROUTINE WAM_SOURCE_init()
implicit none
REAL(real_kind), PARAMETER :: GAMMA = 3.000000
REAL(real_kind), PARAMETER :: SA=7.000000E-02
REAL(real_kind), PARAMETER :: SB=9.000000E-02
INTEGER :: J
REAL(real_kind), ALLOCATABLE :: U10(:), FPK(:), ALPHJO(:)

KL = nDir
ML = nFre
!   PREPARE_SOURCE - ROUTINE TO PREPARE WAM SOURCE MODULE.                     !

FR5(:) = Fre(:)**5
FRM5(:) = 1/FR5(:)

INV_LOG_CO = 1./LOG10(ratioFre)
FMIN = 0.07**2 /(16.*(Fre(nFre)-Fre(1))*PI2)

ALLOCATE (U10(1:JUMAX), FPK(1:JUMAX), ALPHJO(1:JUMAX))
DO J = 1,JUMAX
   U10(J) = REAL(J)*DELU
END DO

ALLOCATE(FLMINFR(1:JUMAX,1:ML))
CALL FETCH_LAW (FETCH, Fre(nFre), U10, ALPHJO, FPK)
ALPHJO(:) = 0.01*ALPHJO

CALL JONSWAP (Fre, ALPHJO, GAMMA, SA, SB, FPK, FLMINFR(:,:))
FLMINFR(:,:) = MAX(FLMINFR(:,:),FLMIN)  !! AVOID TOO SMALL NUMBERS
DEALLOCATE (U10, FPK, ALPHJO)

!     1. WEIGHT OF NON-LINEAR INTERACTION.                                     !
!        ---------------------------------                                     !

CALL NLWEIGT

CALL INIT_SNONLIN

CALL MAKE_SHALLOW_SNL(bottomDepth(1:nCells))

! ---------------------------------------------------------------------------- !
!     2. TABLES AND PRE_COMPUTED CONSTANT FOR WAM PHYSICS.                     !
!        -------------------------------------------------                     !

CALL TABU_SWELLFT

CALL INIT_SDISSP_ARD

CALL INIT_X0TAUHF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ALLOCATE ARRAYS FOR FLUXES.                                           !
!        ---------------------------                                           !

IF (ALLOCATED(PHIOC)) DEALLOCATE(PHIOC)
ALLOCATE (PHIOC(nCells))
PHIOC(:) = 0.0
IF (ALLOCATED(PHIAW)) DEALLOCATE(PHIAW)
ALLOCATE (PHIAW(nCells))
PHIAW(:) = 0.0
IF (ALLOCATED(TAUOC_X)) DEALLOCATE(TAUOC_X)
ALLOCATE (TAUOC_X(nCells))
TAUOC_X(:) = 0.0
IF (ALLOCATED(TAUOC_Y)) DEALLOCATE(TAUOC_Y)
ALLOCATE (TAUOC_Y(nCells))
TAUOC_Y(:) = 0.0
IF (ALLOCATED(PHIBOT)) DEALLOCATE(PHIBOT)
ALLOCATE (PHIBOT(nCells))
PHIBOT(:) = 0.0
IF (ALLOCATED(TAUBOT_X)) DEALLOCATE(TAUBOT_X)
ALLOCATE (TAUBOT_X(nCells))
TAUBOT_X(:) = 0.0
IF (ALLOCATED(TAUBOT_Y)) DEALLOCATE(TAUBOT_Y)
ALLOCATE (TAUBOT_Y(nCells))
TAUBOT_Y(:) = 0.0

ALLOCATE (FL3(nCells, nDir, nFre))
ALLOCATE (USTAR(nCells))
ALLOCATE (Z0(nCells))
ALLOCATE (TAUW(nCells))
ALLOCATE (WSTAR(nCells))
ALLOCATE (ROAIRN(nCells))

USTAR  = 0.0_rk
Z0     = 0.0_rk
TAUW   = 0.0_rk
WSTAR  = 0.0_rk
ROAIRN = rho_air


END SUBROUTINE WAM_SOURCE_init

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE IMPLSCH(FL3, dt, U10, UDIR, TAUW, USTAR, Z0, ROAIRN, WSTAR, DEPTH)

REAL(real_kind),    INTENT(IN)    :: dt             !! Time step of source terms[S].
REAL(real_kind),    INTENT(IN)    :: U10   (:)      !! WIND SPEED [M/S].
REAL(real_kind),    INTENT(IN)    :: UDIR  (:)      !! WIND DIRECTION [RAD].
REAL(real_kind),    INTENT(IN)    :: ROAIRN(:)      !! SURFACE AIR DENSITY [kg/m**3].
REAL(real_kind),    INTENT(IN)    :: WSTAR (:)      !! CONVECTIVE VELOCITY SCALE [m/s]
REAL(real_kind),    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].
REAL(real_kind),    INTENT(INOUT) :: FL3(:,:,:)     !! FREQUENCY SPECTRUM.
REAL(real_kind),    INTENT(INOUT) :: TAUW  (:)      !! WAVE STRESS IN (M/S)**2 
REAL(real_kind),    INTENT(INOUT) :: Z0    (:)      !! ROUGHNESS LENGTH [M].
REAL(real_kind),    INTENT(OUT)   :: USTAR (:)      !! FRICTION VELOCITY [M/S].

! ---------------------------------------------------------------------------- !
!                                                                              !
!   IMPLSCH - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE FUNCTIONS.        !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS             !
!       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.             !
!       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH                  !
!       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.                                 !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS FN+1=FN+DELT*(SN+SN+1)/2.,  !
!       WHERE SN IS THE TOTAL SOURCE FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF   !
!       - ONLY THE DIAGONAL TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE COMPUTED, !
!         THE NONDIAGONAL TERMS ARE NEGLIGIBLE.                                !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !

INTEGER :: K, M, IJ
REAL(real_kind)    :: DELT

REAL(real_kind)    :: FL(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))  !! DIAGONAL MATRIX OF
                                                    !! FUNCTIONAL DERIVATIVE.
REAL(real_kind)    :: SL(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))  !! TOTAL SOURCE FUNCTION.
REAL(real_kind)    :: SL_BOT(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! BOTTOM SOURCE FUNCTION.

INTEGER :: MIJ(SIZE(FL3,1))    !! LAST FREQUENCY INDEX OF PROGNOSTIC PART.
INTEGER :: JU(SIZE(FL3,1))     !! U10 TABLE INDEX.
REAL(real_kind)    :: EMEAN(SIZE(FL3,1))  !! TOTAL ENERGY
REAL(real_kind)    :: FMEAN(SIZE(FL3,1))  !! MEAN FREQUENCY
REAL(real_kind)    :: AKMEAN(SIZE(FL3,1)) !! MEAN WAVENUMBER BASED ON SQRT(1/K)-MOMENT
REAL(real_kind)    :: TEMP(SIZE(FL3,1),SIZE(FL3,3))
REAL(real_kind)    :: DELFL
LOGICAL :: LLWS(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3))

REAL(real_kind)    :: F1MEAN(SIZE(FL3,1))    !! MEAN FREQUENCY BASED ON F-MOMENT
REAL(real_kind)    :: XKMEAN(SIZE(FL3,1))    !! MEAN WAVENUMBER BASED ON SQRT(K)-MOMENT
REAL(real_kind)    :: EMEANWS(SIZE(FL3,1))   !! TOTAL WINDSEA ENERGY
REAL(real_kind)    :: FMEANWS(SIZE(FL3,1))   !! MEAN WINDSEA FREQUENCY
REAL(real_kind)    :: SPRD(SIZE(FL3,1),SIZE(FL3,2)) !! SPREAD FUNCTIONS.

REAL(real_kind)    :: SPOS(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! POSITIVE SINPUT
REAL(real_kind)    :: SMIN(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! NEGATIVE SINPUT
REAL(real_kind)    :: SSOURCE(SIZE(FL3,1),SIZE(FL3,2),SIZE(FL3,3)) !! SOURCE TERMS CONTRIBUTING
                                                        !! SURFACE WAVE FLUXES. 

! MODULATION OF SOURCE TERM BY IMPLICIT FACTOR IN THE CALCULATION
! OF THE SURFACE WAVE FLUXES To THE OCEANS.
LOGICAL, PARAMETER :: LLIMPFLX=.FALSE.


!OPENACC
REAL(real_kind)   :: TEMPS1,TEMPS2,TEMPS3


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. CALCULATE ROUGHNESS LENGTH AND FRICTION VELOCITIES.                   !
!        ---------------------------------------------------                   !


DO  K = 1,KL
   DO IJ = 1, SIZE(UDIR)
      SPRD(IJ,K) = MAX(0.0_rk,COS(Theta(K)-UDIR(IJ)))**2   !! cosine spreading.
   END DO
END DO

DO IJ = 1,SIZE(U10)
   JU(IJ) = MIN(JUMAX, MAX(NINT(U10(IJ)/DELU),1))   !! u10 tabel index.
END DO


! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN PARAMETERS.                                              !
!        ------------------------                                              !

!     2.1 Total wave energy.                                                   !
!CALL TOTAL_ENERGY_B(FL3, EMEAN)
CALL TOTAL_ENERGY_B_OPENACC(FL3, EMEAN)

!     2.2 Reduce wave energy if larger than depth limited wave height.         !

!IF (WAVE_BREAKING_RUN) THEN
  !CALL SDEPTHLIM(DEPTH, EMEAN, FL3) 
  CALL SDEPTHLIM_OPENACC(DEPTH, EMEAN, FL3) 
!ENDIF
 
!     2.3 Mean frequencies and wave numbers.                                   !

!CALL FEMEAN(FL3, EMEAN, FMEAN)
CALL FEMEAN_OPENACC (FL3, EMEAN, FMEAN)

!CALL TM1_TM2_PERIODS_B(FL3, EMEAN, TM1=F1MEAN)
CALL TM1_TM2_PERIODS_B_OPENACC(FL3, EMEAN, TM1=F1MEAN)

DO IJ = 1, SIZE(F1MEAN)
   F1MEAN(IJ) = 1.0_rk/F1MEAN(IJ)
END DO


CALL WM1_WM2_WAVENUMBER_B_OPENACC(FL3, EMEAN, WM1=AKMEAN, WM2=XKMEAN)

! -------------------------------------------------------------------------!
!                                                                          !
!     3. COMPUTATION OF SOURCE FUNCTIONS AND DERIVATIVES.                  !
!        ------------------------------------------------                  !

CALL AIRSEA_OPENACC(U10, TAUW, USTAR, Z0)
 
CALL SINPUT_ARD_OPENACC(FL3, SL, SPOS, FL, USTAR, UDIR, Z0, ROAIRN, WSTAR, LLWS)

CALL TOTAL_ENERGY_B_OPENACC(FL3, EMEANWS, LLWS)
CALL FEMEAN_OPENACC (FL3, EMEANWS, FMEANWS, LLWS)
CALL FRCUTINDEX_OPENACC (FMEAN, FMEANWS, USTAR, MIJ)

CALL STRESSO_OPENACC(FL3, SPOS, USTAR, UDIR, Z0, MIJ, TAUW, PHIAW)

! re-evalute the input
CALL AIRSEA_OPENACC (U10, TAUW, USTAR, Z0)

CALL IMPHFTAIL_OPENACC(MIJ,  FL3)

CALL SINPUT_ARD_OPENACC (FL3, SL, SPOS, FL, USTAR, UDIR, Z0, ROAIRN, WSTAR, LLWS)

CALL STRESSO_OPENACC (FL3, SPOS, USTAR, UDIR, Z0, MIJ, TAUW, PHIAW)

CALL SDISSIP_ARD_OPENACC (FL3, SL, FL, USTAR, UDIR, ROAIRN)

CALL SNONLIN_OPENACC (FL3, SL, FL, DEPTH, AKMEAN)

IF (SHALLOW_RUN) THEN
   DO M = 1,SIZE(SL_BOT,3)
     DO K = 1,SIZE(SL_BOT,2) 
       DO IJ = 1,SIZE(SL_BOT,1)
         SL_BOT(IJ,K,M) = 0.0
       END DO
     END DO
   END DO 
   CALL SBOTTOM_OPENACC (FL3, SL_BOT, FL, DEPTH)
END IF

IF (SHALLOW_RUN) THEN
   DO M = 1,SIZE(SL,3)
      DO K = 1,SIZE(SL,2)
        DO IJ = 1,SIZE(SL,1)
           SL(IJ,K,M) = SL(IJ,K,M) + SL_BOT(IJ,K,M)
        END DO
      END DO
   END DO
END IF

IF (WAVE_BREAKING_RUN) THEN
   CALL SFBRK_OPENACC (FL3, SL, FL, EMEAN, FMEAN, DEPTH)
END IF

!IF (PHILLIPS_RUN) THEN
!  call source_phillips (sl, ustar, udir, depth)
!END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTATION OF NEW SPECTRA.                                           !
!        ---------------------------                                           !
!                                                                              !
!       INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE             !
!       FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.                    !

!     4.2 INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE           !
!         FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.                  !
!         ----------------------------------------------------------           !

DELT = dt
DO M=1,SIZE(SL,3)
   DO K=1,SIZE(SL,2)
      DO IJ = 1,SIZE(SL,1)
         DELFL = 5.0E-07*G/FRE(M)**4*DELT
         TEMPS2 = USTAR(IJ)*DELFL*MAX(FMEANWS(IJ), FMEAN(IJ))
         TEMPS1 = DELT*SL(IJ,K,M)/ MAX(1.0_rk, 1.0_rk-DELT*FL(IJ,K,M))
         TEMPS3 = MIN(ABS(TEMPS1),TEMPS2)
         FL3(IJ,K,M) = FL3(IJ,K,M)+SIGN(TEMPS3,TEMPS1)
         FL3(IJ,K,M) = MAX(FL3(IJ,K,M),FLMINFR(JU(IJ),M)*SPRD(IJ,K))
      END DO
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     5. REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.                 !
!        -----------------------------------------------------                 !

!    5.1 COMPUTE MEAN PARAMETERS.                                              !
!        ------------------------                                              !

CALL TOTAL_ENERGY_B_OPENACC(FL3, EMEAN)
CALL FEMEAN_OPENACC (FL3, EMEAN, FMEAN)
CALL TOTAL_ENERGY_B_OPENACC(FL3, EMEANWS, LLWS)
CALL FEMEAN_OPENACC (FL3, EMEANWS, FMEANWS, LLWS)


!     5.2 COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.         !
!         ------------------------------------------------------------         !

CALL FRCUTINDEX_OPENACC (FMEAN, FMEANWS, USTAR, MIJ)


!     5.3 COMPUTE TAIL ENERGY RATIOS AND MERGE TAIL INTO SPECTRA.              !
!         -------------------------------------------------------              !

CALL IMPHFTAIL_OPENACC(MIJ,  FL3)


END SUBROUTINE IMPLSCH


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE AIRSEA (UTOP, TAUW, USTAR, Z0)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   AIRSEA - - COMPUTATION OF TOTAL STRESS AND ROUGHNESS LENGTH SCALE.         !
!                                                                              !
!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TOTAL STRESS.                                                  !
!                                                                              !
!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)+RNUAIRM/USTAR

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)    :: UTOP(:)       !! WIND SPEED AT REFERENCE LEVEL XNLEV.
REAL(real_kind),    INTENT(IN)    :: TAUW(:)       !! WAVE STRESS.
REAL(real_kind),    INTENT(OUT)   :: USTAR(:)      !! FRICTION VELOCITY.
REAL(real_kind),    INTENT(OUT)   :: Z0(:)         !! ROUGHNESS LENGTH.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER, PARAMETER :: NITER=15

REAL(real_kind), PARAMETER :: TWOXMP1=3.0
REAL(real_kind), PARAMETER :: EPSUS = 1.0E-6

!     *ACD*       COEFFICIENTS FOR SIMPLE CD(U10) RELATION
!     *BCD*       CD = ACD + BCD*U10
REAL(real_kind), PARAMETER :: ACD=8.0E-4
REAL(real_kind), PARAMETER :: BCD=8.0E-5

INTEGER :: IJ, ITER

! REAL :: ALPHA
REAL(real_kind) :: XLOGXL, ALPHAOG, XKUTOP, XOLOGZ0
REAL(real_kind) :: USTOLD, TAUOLD, TAUNEW, X, F, DELF
REAL(real_kind) :: USTM1, Z0TOT, Z0CH, Z0VIS, ZZ

! ----------------------------------------------------------------------

XLOGXL=LOG(XNLEV)
ALPHAOG=ALPHA/G

DO IJ=1,SIZE(USTAR)
   XKUTOP = XKAPPA*UTOP(IJ)
   USTOLD = UTOP(IJ)*SQRT(ACD+BCD*UTOP(IJ))
   TAUOLD = MAX(USTOLD**2,TAUW(IJ)+EPS1)
   USTAR(IJ) = SQRT(TAUOLD)
   USTM1 = 1.0/MAX(USTAR(IJ),EPSUS)

   DO ITER=1,NITER
      X = TAUW(IJ)/TAUOLD
!      Z0CH = ALPHAOG*TAUOLD/SQRT(1.0-X)
      Z0CH = ALPHAOG*TAUOLD/SQRT(MAX(1.0-X,EPS1))
      Z0VIS = RNUAIRM*USTM1
      Z0TOT = Z0CH+Z0VIS

      XOLOGZ0= 1.0/(XLOGXL-LOG(Z0TOT))
      F = USTAR(IJ)-XKUTOP*XOLOGZ0
      ZZ = USTM1*(Z0CH*(2.0-TWOXMP1*X)/(1.0-X)-Z0VIS)/Z0TOT
      DELF= 1.0-XKUTOP*XOLOGZ0**2*ZZ

      USTAR(IJ) = USTAR(IJ)-F/DELF
      TAUNEW = MAX(USTAR(IJ)**2,TAUW(IJ)+EPS1)
      USTAR(IJ) = SQRT(TAUNEW)
      IF (TAUNEW.EQ.TAUOLD) EXIT
      USTM1 = 1.0/MAX(USTAR(IJ),EPSUS)
      TAUOLD = TAUNEW
   ENDDO

   Z0(IJ)=Z0CH

ENDDO

END SUBROUTINE AIRSEA

SUBROUTINE AIRSEA_OPENACC (UTOP, TAUW, USTAR, Z0)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   AIRSEA - - COMPUTATION OF TOTAL STRESS AND ROUGHNESS LENGTH SCALE.         !
!                                                                              !
!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990                             !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TOTAL STRESS.                                                  !
!                                                                              !
!     METHOD.
!     -------

!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH

!                  Z1=Z0/SQRT(1-TAUW/TAU)+RNUAIRM/USTAR

!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)    :: UTOP(:)       !! WIND SPEED AT REFERENCE LEVEL XNLEV.
REAL(real_kind),    INTENT(IN)    :: TAUW(:)       !! WAVE STRESS.
REAL(real_kind),    INTENT(OUT)   :: USTAR(:)      !! FRICTION VELOCITY.
REAL(real_kind),    INTENT(OUT)   :: Z0(:)         !! ROUGHNESS LENGTH.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER, PARAMETER :: NITER=15

REAL(real_kind), PARAMETER :: TWOXMP1=3.0_rk
REAL(real_kind), PARAMETER :: EPSUS = 1.0E-6

!     *ACD*       COEFFICIENTS FOR SIMPLE CD(U10) RELATION
!     *BCD*       CD = ACD + BCD*U10
REAL(real_kind), PARAMETER :: ACD=8.0E-4
REAL(real_kind), PARAMETER :: BCD=8.0E-5

INTEGER :: IJ, ITER

! REAL :: ALPHA
REAL(real_kind) :: XLOGXL, ALPHAOG, XKUTOP, XOLOGZ0
REAL(real_kind) :: USTOLD, TAUOLD, TAUNEW, X, F, DELF
REAL(real_kind) :: USTM1, Z0TOT, Z0CH, Z0VIS, ZZ

! ----------------------------------------------------------------------
XLOGXL=LOG(XNLEV)
ALPHAOG=ALPHA/G

DO IJ=1,SIZE(USTAR)
   XKUTOP = XKAPPA*UTOP(IJ)
   USTOLD = UTOP(IJ)*SQRT(ACD+BCD*UTOP(IJ))
   TAUOLD = MAX(USTOLD**2,TAUW(IJ)+EPS1)
   USTAR(IJ) = SQRT(TAUOLD)
   USTM1 = 1.0/MAX(USTAR(IJ),EPSUS)

   DO ITER=1,NITER
      X = TAUW(IJ)/TAUOLD
!      Z0CH = ALPHAOG*TAUOLD/SQRT(1.0-X)
      Z0CH = ALPHAOG*TAUOLD/SQRT(MAX(1.0-X,EPS1))
      Z0VIS = RNUAIRM*USTM1
      Z0TOT = Z0CH+Z0VIS

      XOLOGZ0= 1.0/(XLOGXL-LOG(Z0TOT))
      F = USTAR(IJ)-XKUTOP*XOLOGZ0
      ZZ = USTM1*(Z0CH*(2.0-TWOXMP1*X)/(1.0-X)-Z0VIS)/Z0TOT
      DELF= 1.0-XKUTOP*XOLOGZ0**2*ZZ

      USTAR(IJ) = USTAR(IJ)-F/DELF
      TAUNEW = MAX(USTAR(IJ)**2,TAUW(IJ)+EPS1)
      USTAR(IJ) = SQRT(TAUNEW)
      
      IF (TAUNEW.EQ.TAUOLD) EXIT
      USTM1 = 1.0/MAX(USTAR(IJ),EPSUS)
      TAUOLD = TAUNEW
   ENDDO

   Z0(IJ)=Z0CH

ENDDO

END SUBROUTINE AIRSEA_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MAKE_SHALLOW_SNL (DEPTH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!      MAKE_SHALLOW_SNL - COMPUTE THE NONLINEAR TRANSFER FUNCTION COEFFICIENTS !
!                         FOR SHALLOW WATER.                                   !
!                                                                              !
!      P. JANSSEN     ECMWF  JUNE 2005                                         !
!      H. GUNTHER     HZG    JANUARY 2015  CYCLE_4.5.4                         !
!                                                                              !
!     REFERENCES.                                                              !

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, M
REAL(real_kind)    :: D, OM, XK

REAL(real_kind), PARAMETER  :: ENH_MAX=10.0_rk !! MAXIMUM COEFFICIENT.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.                   !
!         --------------------------------------------------                   !

IF (.NOT.ALLOCATED(ENH)) ALLOCATE(ENH(SIZE(DEPTH),ML+4))

IF (.NOT.ALLOCATED(ENH)) ALLOCATE(ENH(SIZE(DEPTH),MLSTHG))

IF (ISNONLIN.NE.0) THEN
   DO M = 1,ML
      DO IJ = 1, SIZE(DEPTH)
         D = DEPTH(IJ)
         OM = Pi2*Fre(M)
         XK = wave_numb(OM,D)
         ENH(IJ,M) = MIN(ENH_MAX,TRANSF(XK,D))
      END DO
    END DO
!       NOTE THAT FR IS NOT DEFINED FOR M>ML.
    DO M = ML+1, MLSTHG
       DO IJ = 1, SIZE(DEPTH)
          D = DEPTH(IJ)
          OM = PI2*Fre(ML)*ratioFre**(M-ML)
          XK = wave_numb(OM,D)
          ENH(IJ,M) = MIN(ENH_MAX,TRANSF(XK,D))
       END DO
    END DO
ENDIF

END SUBROUTINE MAKE_SHALLOW_SNL

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!                                                                              !
!     G. PRIVATE MODULE PROCEDURES.                                            !
!                                                                              !
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE NLWEIGT

! ---------------------------------------------------------------------------- !
!                                                                              !
!   NLWEIGT - COMPUTATION OF INDEX ARRAYS AND WEIGHTS FOR THE COMPUTATION OF   !
!             THE NONLINEAR TRANSFER RATE.                                     !
!                                                                              !
!     SUSANNE HASSELMANN JUNE 86.                                              !
!                                                                              !
!     H. GUNTHER   ECMWF/GKSS  DECEMBER 90 - CYCLE_4 MODIFICATIONS.            !
!                                            4 FREQUENCIES ADDED.              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF PARAMETERS USED IN DISCRETE INTERACTION                 !
!       PARAMETERIZATION OF NONLINEAR TRANSFER.                                !
!                                                                              !
!     *PARAMETER*  FOR DISCRETE APPROXIMATION OF NONLINEAR TRANSFER.           !

REAL(real_kind), PARAMETER :: ALAMD   = 0.25     !! LAMBDA
REAL(real_kind), PARAMETER :: CON     = 3000.    !! WEIGHT FOR DISCRETE APPROXIMATION OF
                                      !! NONLINEAR TRANSFER
!REAL, PARAMETER :: DELPHI1 = -11.48   !!
!REAL, PARAMETER :: DELPHI2 = 33.56    !!

INTEGER :: KLP1, IC, KH, KLH, K, KS, ISG, K1, K11, K2, K21
INTEGER :: M, IKN, I, ISP, ISM
INTEGER, ALLOCATABLE :: JA1(:,:)
INTEGER, ALLOCATABLE :: JA2(:,:)

REAL(real_kind)    :: DELPHI1
REAL(real_kind)    :: DELPHI2
REAL(real_kind)    :: DELTHA, CL1, CL2, AL11, AL12, CH, CL1H, CL2H
REAL(real_kind)    :: F1P1, FRG, FLP, FLM, FKP, FKM, XF, COSTH3, COSTH4
REAL(real_kind),    ALLOCATABLE :: FRLON(:)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     0. ALLOCATE ARRAYS.                                                      !
!        ----------------                                                      !

F1P1 = LOG10(ratioFre)
ISP = INT(LOG10(1.+ALAMD)/F1P1+.000001)
ISM = FLOOR(LOG10(1.-ALAMD)/F1P1+.0000001)

MFRSTLW = 1+ISM
MLSTHG = ML-ISM

KFRH=-ISM+ISP+2

ALLOCATE(JA1(KL,2))
ALLOCATE(JA2(KL,2))
ALLOCATE(FRLON(MFRSTLW:ML+KFRH))

IF (.NOT.ALLOCATED (IKP ))  ALLOCATE (IKP (MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (IKP1))  ALLOCATE (IKP1(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (IKM ))  ALLOCATE (IKM (MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (IKM1))  ALLOCATE (IKM1(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (K1W ))  ALLOCATE (K1W (KL,2))
IF (.NOT.ALLOCATED (K2W ))  ALLOCATE (K2W (KL,2))
IF (.NOT.ALLOCATED (K11W))  ALLOCATE (K11W(KL,2))
IF (.NOT.ALLOCATED (K21W))  ALLOCATE (K21W(KL,2))
IF (.NOT.ALLOCATED (AF11))  ALLOCATE (AF11(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FKLAP ))  ALLOCATE (FKLAP(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FKLAP1))  ALLOCATE (FKLAP1(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FKLAM ))  ALLOCATE (FKLAM(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FKLAM1))  ALLOCATE (FKLAM1(MFRSTLW:MLSTHG))
IF (.NOT.ALLOCATED (FRH))  ALLOCATE(FRH(KFRH))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTATION FOR ANGULAR GRID.                                         !
!        -----------------------------                                         !

XF      = ((1.+ALAMD)/(1.-ALAMD))**4
COSTH3  = (1.+2.*ALAMD+2.*ALAMD**3)/(1.+ALAMD)**2
DELPHI1 = -180./PI*ACOS(COSTH3)
COSTH4  = SQRT(1.-XF+XF*COSTH3**2)
DELPHI2 = 180./PI*ACOS(COSTH4)

DELTHA = delTheta*DEG
CL1 = DELPHI1/DELTHA
CL2 = DELPHI2/DELTHA

!     1.1 COMPUTATION OF INDICES OF ANGULAR CELL.                              !
!         ---------------------------------------                              !

KLP1 = KL+1
IC = 1

DO KH = 1,2
   KLH = KL
   IF (KH.EQ.2) KLH=KLP1
   DO K = 1,KLH
      KS = K
      IF (KH.GT.1) KS=KLP1-K+1
      IF (KS.GT.KL) CYCLE
      CH = IC*CL1
      JA1(KS,KH) = JAFU(CH,K,KL)
      CH = IC*CL2
      JA2(KS,KH) = JAFU(CH,K,KL)
   END DO
   IC = -1
END DO

!     1.2 COMPUTATION OF ANGULAR WEIGHTS.                                      !
!         -------------------------------                                      !

CL1  = CL1-INT(CL1)
CL2  = CL2-INT(CL2)
ACL1 = ABS(CL1)
ACL2 = ABS(CL2)
CL11 = 1.-ACL1
CL21 = 1.-ACL2
AL11 = (1.+ALAMD)**4
AL12 = (1.-ALAMD)**4
DAL1 = 1./AL11
DAL2 = 1./AL12

!     1.3 COMPUTATION OF ANGULAR INDICES.                                      !
!         -------------------------------                                      !

ISG = 1
DO KH = 1,2
   CL1H = ISG*CL1
   CL2H = ISG*CL2
   DO K = 1,KL
      KS = K
      IF (KH.EQ.2) KS = KL-K+2
      IF(K.EQ.1) KS = 1
      K1 = JA1(K,KH)
      K1W(KS,KH) = K1
      IF (CL1H.LT.0.) THEN
         K11 = K1-1
         IF (K11.LT.1) K11 = KL
      ELSE
         K11 = K1+1
         IF (K11.GT.KL) K11 = 1
      END IF
      K11W(KS,KH) = K11
      K2 = JA2(K,KH)
      K2W(KS,KH) = K2
      IF (CL2H.LT.0) THEN
         K21 = K2-1
         IF(K21.LT.1) K21 = KL
      ELSE
         K21 = K2+1
         IF (K21.GT.KL) K21 = 1
      END IF
      K21W(KS,KH) = K21
   END DO
   ISG = -1
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTATION FOR FREQUENCY GRID.                                       !
!        -------------------------------                                       !

FRLON(1:ML) = FRE(1:ML)

DO M=0,MFRSTLW,-1
   FRLON(M)=FRLON(M+1)/ratioFre
ENDDO
DO M=ML+1,ML+KFRH
   FRLON(M) = ratioFre*FRLON(M-1)
ENDDO

DO M = MFRSTLW,MLSTHG
   FRG = FRLON(M)
   AF11(M) = CON * FRG**11
   FLP = FRG*(1.+ALAMD)
   FLM = FRG*(1.-ALAMD)
   IKN = M+ISP
   IKP(M) = IKN
   FKP = FRLON(IKP(M))
   IKP1(M) = IKP(M)+1
   FKLAP(M) = (FLP-FKP)/(FRLON(IKP1(M))-FKP)
   FKLAP1(M) = 1.-FKLAP(M)
   IF (FRLON(MFRSTLW).GE.FLM) THEN
      IKM(M) = 1
      IKM1(M) = 1
      FKLAM(M) = 0.
      FKLAM1(M) = 0.
   ELSE
      IKN = M+ISM
      IKM(M) = IKN
      FKM = FRLON(IKM(M))
      IKM1(M) = IKM(M)+1
      FKLAM(M) = (FLM-FKM)/(FRLON(IKM1(M))-FKM)
      FKLAM1(M) = 1.-FKLAM(M)
      IF (IKN.LT.MFRSTLW) THEN
         IKM(M) = 1
!         IKM1(M) = 1
         FKLAM1(M) = 0.
      ENDIF
   END IF

   IF(IKM(M) .gt. nFre) IKM(M) = nFre
   IF(IKM1(M) .gt. nFre) IKM1(M) = nFre
   IF(IKM(M) .lt. 1) IKM(M) = 1
   IF(IKM1(M) .lt. 1) IKM1(M) =1
#ifdef MPI_DEBUG_OUTPUT
   IF (IKM1(M)<-2) THEN
     WRITE(*,*) "out of range IKM1(M)=",IKM1(M),",M=",M,",mpi_rank=",mpi_rank
     WRITE(*,*) "out of range IKM(M)=",IKM(M),",M=",M,",mpi_rank=",mpi_rank
     WRITE(*,*) "out of range MFRSTLW=",MFRSTLW,"mpi_rank=",mpi_rank
     WRITE(*,*) "out of range MLSTHG=",MLSTHG,"mpi_rank=",mpi_rank
     WRITE(*,*) "out of range ISM=",ISM,"mpi_rank=",mpi_rank
   END IF
#endif
END DO


!                                                                              !
!     3. COMPUTE TAIL FREQUENCY RATIOS.                                        !
!        ------------------------------                                        !

DO I=1,KFRH
   M = ML+I-1
   FRH(I) = (FRLON(ML)/FRLON(M))**5
END DO

END SUBROUTINE NLWEIGT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INIT_SNONLIN

! ---------------------------------------------------------------------------- !
!                                                                              !
!    INIT_SNONLIN - INITIALISE ALL FREQUENCY DEPENDENT ARRAYS USED BY SNONLIN  !
!                                                                              !
!     J. BIDLOT   ECMWF  MAY 2012                                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       USED TO BE IN SNONLIN BUT NOW IT IS ONLY COMPUTED ONCE.                !
! ---------------------------------------------------------------------------- !

INTEGER :: ICOUNT, IRCOUNT
INTEGER :: MC, MP, MP1, MM, MM1, IC, IP, IP1, IM , IM1, ITEMP

REAL :: ALPH, FRR
REAL :: FFACP, FFACP1, FFACM, FFACM1, FTAIL, FKLAMP, FKLAMP1
REAL :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2
REAL :: FKLAP12, FKLAP22, FKLAMM, FKLAMM1, FKLAMMA, FKLAMMB
REAL :: FKLAMM2, FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
REAL :: GW1, GW2, GW3, GW4, GW5, GW6, GW7, GW8

!     INLINE FUNCTION (PIERSON-MOSKOWITZ SMOOTH CUT-OFF)
!     X == FR(1)/FREQUENCY
REAL :: EPMMA, X
EPMMA(X) = EXP(-MIN(1.25*X**4,50.))*(X**5)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. FRONT SPECTRAL TAIL REDUCTION COEFFICIENTS

IF(.NOT.ALLOCATED(FTRF)) ALLOCATE(FTRF(MFRSTLW:1))
ALPH = 1./EPMMA(1.)
FRR = 1.
DO MC=1,MFRSTLW,-1
   FTRF(MC)=ALPH*EPMMA(FRR)
   FRR=FRR*ratioFre
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. WORK ARRAYS STORING THE DIFFERENT INDICES AND COEFFICIENTS

IF(.NOT.ALLOCATED(INLCOEF)) ALLOCATE(INLCOEF(NINL,1:MLSTHG))
IF(.NOT.ALLOCATED(RNLCOEF)) ALLOCATE(RNLCOEF(NRNL,1:MLSTHG))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. FREQUENCY LOOP.
!        ---------------

DO MC=1,MLSTHG
   MP  = IKP (MC)
   MP1 = IKP1(MC)
   MM  = IKM (MC)
   MM1 = IKM1(MC)
   FFACP  = 1.
   FFACP1 = 1.
   FFACM  = 1.
   FFACM1 = 1.
   FTAIL  = 1.
   IC  = MC
   IP  = MP
   IP1 = MP1
   IM  = MM
   IM1 = MM1
!       LOW FREQUENCY FRONT TAIL
   IF (IM.LT.1) THEN
      FFACM = FTRF(IM)
      IM = 1
      IF (IM1.LT.1) THEN
         FFACM1 = FTRF(IM1)
         IM1 = 1
      ENDIF
   ENDIF
!       HIGH FREQUENCY TAIL
   IF (IP1.GT.ML) THEN
! Quick fix from Deborah
      ITEMP=IP1-ML+1
      IF(ITEMP .GT. SIZE(FRH))THEN
         ITEMP=SIZE(FRH)
      ENDIF
!         FFACP1 = FRH(IP1-ML+1)
      FFACP1 = FRH(ITEMP)

      IP1 = ML
      IF (IP .GT.ML) THEN
         FFACP  = FRH(IP -ML+1)
         IP  = ML
         IF (IC .GT.ML) THEN
            FTAIL  = FRH(IC -ML+1)
            IC  = ML
            IF (IM1.GT.ML) THEN
               FFACM1 = FRH(IM1-ML+1)
               IM1 = ML
            ENDIF
         ENDIF
      ENDIF
   ENDIF

   ICOUNT=1
   INLCOEF(ICOUNT,MC) = IC
   ICOUNT=ICOUNT+1
   INLCOEF(ICOUNT,MC) = IP
   ICOUNT=ICOUNT+1
   INLCOEF(ICOUNT,MC) = IP1
   ICOUNT=ICOUNT+1
   INLCOEF(ICOUNT,MC) = IM
   ICOUNT=ICOUNT+1
   INLCOEF(ICOUNT,MC) = IM1

   FKLAMP  = FKLAP(MC)
   FKLAMP1 = FKLAP1(MC)
   GW2 = FKLAMP1*FFACP*DAL1
   GW1 = GW2*CL11
   GW2 = GW2*ACL1
   GW4 = FKLAMP*FFACP1*DAL1
   GW3 = GW4*CL11
   GW4 = GW4*ACL1
   FKLAMPA = FKLAMP*CL11
   FKLAMPB = FKLAMP*ACL1
   FKLAMP2 = FKLAMP1*ACL1
   FKLAMP1 = FKLAMP1*CL11
   FKLAPA2 = FKLAMPA**2
   FKLAPB2 = FKLAMPB**2
   FKLAP12 = FKLAMP1**2
   FKLAP22 = FKLAMP2**2
   IRCOUNT=1
   RNLCOEF(IRCOUNT,MC) = FTAIL
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW1
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW3
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW4
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMPA
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMPB
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMP2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMP1
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAPA2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAPB2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAP12
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAP22

   FKLAMM  = FKLAM(MC)
   FKLAMM1 = FKLAM1(MC)
   GW6 = FKLAMM1*FFACM*DAL2
   GW5 = GW6*CL21
   GW6 = GW6*ACL2
   GW8 = FKLAMM*FFACM1*DAL2
   GW7 = GW8*CL21
   GW8 = GW8*ACL2
   FKLAMMA = FKLAMM*CL21
   FKLAMMB = FKLAMM*ACL2
   FKLAMM2 = FKLAMM1*ACL2
   FKLAMM1 = FKLAMM1*CL21
   FKLAMA2 = FKLAMMA**2
   FKLAMB2 = FKLAMMB**2
   FKLAM12 = FKLAMM1**2
   FKLAM22 = FKLAMM2**2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW5
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW6
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW7
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = GW8
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMMA
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMMB
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMM2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMM1
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMA2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAMB2
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAM12
   IRCOUNT=IRCOUNT+1
   RNLCOEF(IRCOUNT,MC) = FKLAM22

ENDDO

IF(ICOUNT.NE.NINL) THEN
   print*, 'ERROR IN INISNONLIN : ICOUNT NE NINL'
   print*, 'ICOUNT= ',ICOUNT
   print*, 'NINL= ',NINL
   print*, '*************************************'
   stop
ENDIF
IF(IRCOUNT.NE.NRNL) THEN
   print*, '*************************************'
   print*, 'ERROR IN INISNONLIN : IRCOUNT NE NRNL'
   print*, 'IRCOUNT= ',IRCOUNT
   print*, 'NRNL= ',NRNL
   print*, '*************************************'
   stop
ENDIF

END SUBROUTINE INIT_SNONLIN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SBOTTOM (F, SL, FL, DEPTH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SBOTTOM - COMPUTATION OF BOTTOM FRICTION.                                  !
!                                                                              !
REAL(real_kind), INTENT(IN)     :: F (:, :, :)  !! SPECTRUM.
REAL(real_kind), INTENT(INOUT)  :: SL(:, :, :)  !! TOTAL SOURCE FUNCTION ARRAY
REAL(real_kind), INTENT(INOUT)  :: FL(:, :, :)  !! DIAGONAL MATRIX OF FUNC. DERIVATIVE.
REAL(real_kind), INTENT(IN)     :: DEPTH(:)     !! WATER DEPTH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL(real_kind), PARAMETER  :: CONST = -2.0*0.038/G

INTEGER :: M, K, IJ
REAL(real_kind)    :: WAV(SIZE(F,1))
REAL(real_kind)    :: SBO(SIZE(F,1))

! ---------------------------------------------------------------------------- !

FRE: DO M = 1,SIZE(F,3)
    do IJ = 1, SIZE(F,1)
       WAV(IJ) = wave_numb(sigma(m), DEPTH(IJ))
       SBO(IJ) = MIN (2.* DEPTH(IJ)*WAV(IJ) ,50.0_rk)
       SBO(IJ) = CONST*WAV(IJ)/SINH(SBO(IJ))
    END DO

    DIR: DO K = 1,SIZE(F,2)
       SL(:,K,M) = SL(:,K,M) + SBO*F(:,K,M)
       FL(:,K,M) = FL(:,K,M) + SBO
   END DO DIR
END DO FRE

END SUBROUTINE SBOTTOM

SUBROUTINE SBOTTOM_OPENACC (F, SL, FL, DEPTH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SBOTTOM - COMPUTATION OF BOTTOM FRICTION.                                  !
!                                                                              !
!     G.J.KOMEN AND Q.D.GAO                                                    !
!     OPTIMIZED BY L.F. ZAMBRESKY                                              !
!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90                            !
!     E. MYKLEBUST        FEBRUARY 2005       OPTIMIZATION                     !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF BOTTOM FRICTION DISSIPATION                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       SEE REFERENCES.                                                        !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)                 !
!       BOUWS AND KOMEN, JPO 13(1983)1653-1658                                 !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind), INTENT(IN)     :: F (:, :, :)  !! SPECTRUM.
REAL(real_kind), INTENT(INOUT)  :: SL(:, :, :)  !! TOTAL SOURCE FUNCTION ARRAY
REAL(real_kind), INTENT(INOUT)  :: FL(:, :, :)  !! DIAGONAL MATRIX OF FUNC. DERIVATIVE.
REAL(real_kind), INTENT(IN)     :: DEPTH(:)     !! WATER DEPTH
! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL(real_kind), PARAMETER  :: CONST = -2.0*0.038/G

INTEGER :: M, K, IJ, I
REAL(real_kind)    :: WAV
REAL(real_kind)    :: SBO

! ---------------------------------------------------------------------------- !

DO M = 1,SIZE(F,3)
   DO K = 1,SIZE(F,2)
      DO IJ = 1, SIZE(F,1)
         I = cellIdxTab(IJ)
         WAV = kwave(I,M)
         SBO = MIN (2.* DEPTH(IJ)*WAV ,50.0_rk)
         SBO = CONST*WAV/SINH(SBO)
         SL(IJ,K,M) = SL(IJ,K,M) + SBO*F(IJ,K,M)
         FL(IJ,K,M) = FL(IJ,K,M) + SBO
      END DO
   END DO
END DO

END SUBROUTINE SBOTTOM_OPENACC

SUBROUTINE SDISSIP_ARD_OPENACC (F, SL, FL, USTAR, UDIR, ROAIRN)

! ---------------------------------------------------------------------------- !
!**** *SDISSIP_ARD* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     LOTFI AOUF       METEO FRANCE 2013
!     FABRICE ARDHUIN  IFREMER  2013


!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       ARDHUIN et AL. JPO DOI:10.1175/20110JPO4324.1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)    :: F (:, :, :)    !! SPECTRUM.
REAL(real_kind),    INTENT(OUT)   :: SL(:, :, :)    !! TOTAL SOURCE FUNCTION ARRAY
REAL(real_kind),    INTENT(OUT)   :: FL(:, :, :)    !! DIAGONAL MATRIX OF FUNCTIONAL
                                                    !! DERIVATIVE
REAL(real_kind),    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL(real_kind),    INTENT(IN)    :: UDIR (:)       !! WIND DIRECTION.
REAL(real_kind),    INTENT(IN)    :: ROAIRN(:)      !! AIR DENSITY

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------

INTEGER :: IJ, K, M, M2, K2, KK, KLD, I
INTEGER, DIMENSION(SIZE(F,2),SIZE(F,2)) :: KKD

REAL(real_kind) :: XK(SIZE(F,1),SIZE(F,3))
REAL(real_kind) :: TPIINV, TPIINVH, TMP01, TMP03
REAL(real_kind) :: EPSR
REAL(real_kind) :: ROG
REAL(real_kind) :: SSDSC6M1
REAL(real_kind) :: SIG(SIZE(F,3))
REAL(real_kind) :: SSDSC2_SIG(SIZE(F,3))
REAL(real_kind) :: FACTURB(SIZE(F,1))
REAL(real_kind) :: FACSAT(SIZE(F,1),SIZE(F,3))
REAL(real_kind) :: FACWTRB(SIZE(F,1),SIZE(F,3))
REAL(real_kind) :: TEMP1(SIZE(F,1),SIZE(F,3))
REAL(real_kind) :: BTH0(SIZE(F,1),SIZE(F,3))  !saturation spectrum 
REAL(real_kind) :: BTH(SIZE(F,1),SIZE(F,2),SIZE(F,3))  !saturation spectrum 

REAL(real_kind), DIMENSION(SIZE(F,1),SIZE(F,2),SIZE(F,3)) :: D

!!! the following 3 arrays are only used when  LLSSDSC3
!!! is not used, there should be a way to save the memory
REAL(real_kind) :: SCUMUL 
REAL(real_kind) :: RENEWALFREQ
REAL(real_kind) :: WCUMUL

LOGICAL :: LLSSDSC3,  LLSSDSC5

INTEGER :: TEMPC



! INITIALISATION

EPSR=SQRT(SDSBR)

TPIINV = 1.0/PI2
TPIINVH= 0.5*TPIINV

KLD=SIZE(F,2)/2

ROG = rho_sea*G

LLSSDSC3=(SSDSC3.NE.0.0)

TMP03 = 1.0/(SDSBR*MICHE)

LLSSDSC5=(SSDSC5.NE.0.0)

SSDSC6M1=1.0-SSDSC6

DO M = 1,SIZE(F,3)
  SIG(M) = PI2*FRE(M)
  SSDSC2_SIG(M)=SSDSC2*SIG(M)
END DO

IF (SHALLOW_RUN) THEN
  DO M = 1,SIZE(F,3)
    DO IJ = 1,SIZE(F,1)
       I = cellIdxTab(IJ)
       XK(IJ,M) = kwave(I,M)
       FACSAT(IJ,M) = XK(IJ,M)**3*TPIINV*Cg(I,M)
    ENDDO
  ENDDO
ELSE
  DO M= 1,SIZE(F,3)
    DO IJ = 1,SIZE(F,1)
      XK(IJ,M) = (SIG(M)**2)/G
      FACSAT(IJ,M) = XK(IJ,M)**3*TPIINVH*G/SIG(M)
    ENDDO
  ENDDO
ENDIF


! COMPUTE SATURATION SPECTRUM
DO M = 1,SIZE(F,3)
   DO IJ = 1,SIZE(F,1)
      BTH0(IJ,M) = 0.0
   ENDDO
ENDDO

DO M = 1,SIZE(F,3)
  DO K = 1,SIZE(F,2)
    DO IJ = 1,SIZE(F,1)
      BTH(IJ,K,M)= 0.0
      ! integrates in directional sector
      DO K2 = 1,NSDSNTH*2+1
        KK=INDICESSAT(K,K2)
        BTH(IJ,K,M) = BTH(IJ,K,M) + SATWEIGHTS(K,K2)*F(IJ,KK,M)
      ENDDO
      BTH(IJ,K,M)=BTH(IJ,K,M)*FACSAT(IJ,M)
    ENDDO
  ENDDO
ENDDO

DO M = 1,SIZE(F,3)
  DO IJ = 1,SIZE(F,1)
    DO K = 1,SIZE(F,2)
      BTH0(IJ,M)=MAX(BTH0(IJ,M),BTH(IJ,K,M))
    ENDDO
  ENDDO
ENDDO


! SATURATION TERM

DO  M = 1,SIZE(F,3)
  DO IJ = 1,SIZE(F,1)
    TEMP1(IJ,M)=SSDSC6*(MAX(0.0_rk,BTH0(IJ,M)*TMP03-SSDSC4))**IPSAT
  ENDDO
ENDDO

DO  M = 1,SIZE(F,3)
  DO K = 1,SIZE(F,2)
    DO IJ = 1,SIZE(F,1)
      D(IJ,K,M)= SSDSC2_SIG(M)*(TEMP1(IJ,M)+SSDSC6M1*(MAX(0.0_rk,BTH(IJ,K,M)*TMP03-SSDSC4))**IPSAT)
    ENDDO
  ENDDO
ENDDO

!!  ! CUMULATIVE TERM
!!  IF (LLSSDSC3) THEN
!!  
!!    DO M2 = 1,SIZE(F,3)-NDIKCUMUL
!!      DO IJ = 1,SIZE(F,1)
!!        IF(BTH0(IJ,M2).GT.SDSBR) THEN
!!          TEMP1(IJ,M2)=1.0
!!        ELSE
!!          TEMP1(IJ,M2)=0.0
!!        ENDIF
!!      ENDDO
!!    ENDDO
!!  
!!  !  DO M2 = 1,SIZE(F,3)-NDIKCUMUL
!!  !    DO K2 = 1,SIZE(F,2)
!!  !      DO IJ = 1,SIZE(F,1)
!!  !        SCUMUL(IJ,K2,M2)=TEMP1(IJ,M2)*(MAX(SQRT(BTH(IJ,K2,M2))-EPSR,0.0))**2
!!  !      ENDDO
!!  !    ENDDO
!!  !  ENDDO
!!  
!!  
!!          DO K = 1,SIZE(F,2)
!!            DO K2 = 1,SIZE(F,2)
!!              KKD(K2,K)=ABS(K2-K)
!!              IF(KKD(K2,K).GT.KLD) KKD(K2,K)=KKD(K2,K)-KLD
!!            ENDDO
!!          ENDDO
!!  
!!    DO M = NDIKCUMUL+1,SIZE(F,3)
!!      DO K = 1,SIZE(F,2)
!!        DO IJ = 1,SIZE(F,1)
!!  
!!          RENEWALFREQ=0.0
!!          ! OPENACC: to avoid accelerator restriction : live out variable
!!          TEMPC = M
!!  
!!        ! Correction of saturation level for shallow-water kinematics
!!        ! Cumulative effect based on lambda   (breaking probability is
!!        ! the expected rate of sweeping by larger breaking waves)
!!  
!!          DO M2 = 1,TEMPC-NDIKCUMUL
!!            !DO KK = 0,KLD
!!            !  WCUMUL(IJ,KK,M2)=CUMULW(INDEP(IJ),KK,M2,M)
!!            !ENDDO
!!            DO K2 = 1,SIZE(F,2)
!!              KK=KKD(K2,K)
!!              ! Integrates over frequencies M2 and directions K2 to 
!!              ! Integration is performed from M2=1 to a frequency lower than M: IK-NDIKCUMUL
!!              WCUMUL=CUMULW(INDEP(IJ),KK,M2,M)
!!              SCUMUL=TEMP1(IJ,M2)*(MAX(SQRT(BTH(IJ,K2,M2))-EPSR,0.0))**2
!!              RENEWALFREQ=RENEWALFREQ + WCUMUL*SCUMUL
!!            ENDDO
!!          ENDDO
!!  
!!          D(IJ,K,M)= D(IJ,K,M) + RENEWALFREQ
!!  
!!        ENDDO
!!      ENDDO
!!    ENDDO
!!  
!!  ENDIF  ! LLSSDSC3


!     WAVE-TURBULENCE INTERACTION TERM
IF (LLSSDSC5) THEN
  TMP01 = 2.*SSDSC5/ROG
  DO IJ = 1,SIZE(F,1)
    FACTURB(IJ) = TMP01*ROAIRN(IJ)*USTAR(IJ)*USTAR(IJ)
  ENDDO
  DO M= 1, SIZE(F,3)
    DO IJ = 1,SIZE(F,1)
      FACWTRB(IJ,M) = SIG(M)*XK(IJ,M)*FACTURB(IJ)
    ENDDO
  ENDDO 
  DO M= 1, SIZE(F,3)
    DO K = 1,SIZE(F,2)
      DO IJ = 1,SIZE(F,1)
        D(IJ,K,M)= D(IJ,K,M)- FACWTRB(IJ,M)*COS(UDIR(IJ)-THETA(K))
      ENDDO
    ENDDO
  ENDDO
ENDIF


! ADD ALL CONTRIBUTIONS TO SOURCE TERM
DO  M= 1, SIZE(F,3)
  DO K= 1, SIZE(F,2)
    DO IJ = 1,SIZE(F,1)
      SL(IJ,K,M) = SL(IJ,K,M)+D(IJ,K,M)*F(IJ,K,M)
      FL(IJ,K,M) = FL(IJ,K,M)+D(IJ,K,M)
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE SDISSIP_ARD_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SFBRK (F, SL, FL, EMEAN, FMEAN, DEPTH)

! ---------------------------------------------------------------------------- !
!
!     WEIMIN LUO, POL, MAY 1996, COMPUTATION OF WAVE BREAKING
!
!     PURPOSE
!     -------
!
!     COMPUTE DISSIPATION DUE TO DEPTH-INDUCED WAVE BREAKING
!
!     METHOD.
!     -------
!
!       SEE REFERENCES.
!
!     REFERENCES.
!     -----------
!
!     BATTJES & JANSSEN (1978)
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind), INTENT(IN)     :: F (:, :, :)  !! SPECTRUM.
REAL(real_kind), INTENT(INOUT)  :: SL(:, :, :)  !! TOTAL SOURCE FUNCTION ARRAY
REAL(real_kind), INTENT(INOUT)  :: FL(:, :, :)  !! DIAGONAL MATRIX OF FUNC. DERIVATIVE.
REAL(real_kind), INTENT(IN)     :: EMEAN (:)    !! TOTAL ENERGY
REAL(real_kind), INTENT(IN)     :: FMEAN(:)     !! MEAN FREQUENCY
REAL(real_kind), INTENT(IN)     :: DEPTH(:)     !! WATER DEPTH

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL(real_kind), PARAMETER :: ALPHA = 1.0

INTEGER :: M, K
REAL(real_kind) :: QB(SIZE(F,1)), BB(SIZE(F,1))
REAL(real_kind) :: SBR(SIZE(F,1)),DSBR(SIZE(F,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!   1. compute total dissipation rate according to Battjes-Janssen
!      -----------------------------------------------------------

BB = 8.*EMEAN/(GAMD*DEPTH)**2 !! (Hrms / Hmax)**2

CALL CMPQB (BB, QB)           !! fraction of breaking waves

QB = MIN(1.0_rk,QB)
SBR = -ALPHA*2.*FMEAN

WHERE (BB.LE.1.) SBR = SBR*QB/BB

WHERE (BB .LT. 1. .AND. ABS(BB - QB) .GT. 0.)
   DSBR = SBR * (1. - QB) / (BB - QB)
ELSEWHERE
   DSBR = 0.
ENDWHERE

DO M = 1,SIZE(F,3)
  DO K = 1,SIZE(F,2)
     SL(:,K,M) = SL(:,K,M) + SBR*F(:,K,M)
     FL(:,K,M) = FL(:,K,M) + DSBR
   END DO
END DO

END SUBROUTINE SFBRK

SUBROUTINE SFBRK_OPENACC (F, SL, FL, EMEAN, FMEAN, DEPTH)


REAL(real_kind), INTENT(IN)     :: F (:, :, :)  !! SPECTRUM.
REAL(real_kind), INTENT(INOUT)  :: SL(:, :, :)  !! TOTAL SOURCE FUNCTION ARRAY
REAL(real_kind), INTENT(INOUT)  :: FL(:, :, :)  !! DIAGONAL MATRIX OF FUNC. DERIVATIVE.
REAL(real_kind), INTENT(IN)     :: EMEAN (:)    !! TOTAL ENERGY
REAL(real_kind), INTENT(IN)     :: FMEAN(:)     !! MEAN FREQUENCY
REAL(real_kind), INTENT(IN)     :: DEPTH(:)     !! WATER DEPTH

REAL(real_kind), PARAMETER :: ALPHA = 1.0

INTEGER :: M, K, IJ
REAL(real_kind) :: QB, BB
REAL(real_kind) :: SBR(SIZE(F,1)),DSBR(SIZE(F,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!   1. compute total dissipation rate according to Battjes-Janssen
!      -----------------------------------------------------------

DO IJ = 1,SIZE(F,1)
   BB = 8.*EMEAN(IJ)/(GAMD*DEPTH(IJ))**2 !! (Hrms / Hmax)**2
   CALL CMPQB_OPENACC (BB, QB)           !! fraction of breaking waves
   QB = MIN(1.0_rk, QB)
   SBR(IJ) = -ALPHA*2.*FMEAN(IJ)
   IF (BB.LE.1.) SBR(IJ) = SBR(IJ)*QB/BB
   IF (BB.LT.1. .AND. ABS(BB - QB) .GT. 0.0) THEN
       DSBR(IJ) = SBR(IJ) * (1. - QB) / (BB - QB)
   ELSE
       DSBR(IJ) = 0.0
   END IF
END DO

DO M = 1,SIZE(F,3)
  DO K = 1,SIZE(F,2)
    DO IJ = 1,SIZE(F,1)
      SL(IJ,K,M) = SL(IJ,K,M) + SBR(IJ)*F(IJ,K,M)
      FL(IJ,K,M) = FL(IJ,K,M) + DSBR(IJ)
    END DO
  END DO
END DO

END SUBROUTINE SFBRK_OPENACC
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE CMPQB (BB, QB)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   CMPQB - 
!
!     ADDED BY WEIMIN LUO, POL, MAY 1996
!     BASED ON THE CODE OF G. Ph. van Vledder, Delft Hydraulics
!
!  1. Purpose
!
!     Compute fraction of breaking waves for use in
!     SFBRK wave breaking dissipation function
!
!  2. Method
!
!     Newton-Raphson implementation of
!
!     1 - QB
!     ------ = - (HRMS/HMAX)^2
!     ln(QB)
!
!  3. Parameter list
!
!  4. Subroutines used
!
!  5. Error messages
!
!  6. Remarks
!
!     If HRMS > HMAX --> QB = 1
!
!  7. Structure
!
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

IMPLICIT NONE

REAL(real_kind)    , INTENT(IN)  :: BB(:)     !! (RMS wave height/Maximum wave height)**2 !
REAL(real_kind)    , INTENT(OUT) :: QB(:)     !! Fraction of breaking waves

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL(real_kind)    :: Q0(SIZE(BB))

! ---------------------------------------------------------------------------- !
!                                                                              !

Q0 = 0.
WHERE (BB.GE.0.25) Q0 = (2.*SQRT(BB)-1.)**2

! ---------------------------------------------------------------------------- !
!                                                                              !

WHERE (BB.LT.1.)
   QB = Q0 - BB*(Q0 - exp((Q0-1.)/BB))/(BB-exp((Q0-1.)/BB))
ELSEWHERE
    QB = 1.
END WHERE

END SUBROUTINE CMPQB

SUBROUTINE CMPQB_OPENACC (BB, QB)

IMPLICIT NONE

REAL(real_kind)    , INTENT(IN)  :: BB     !! (RMS wave height/Maximum wave height)**2 !
REAL(real_kind)    , INTENT(OUT) :: QB     !! Fraction of breaking waves

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !


REAL(real_kind)    :: Q0

! ---------------------------------------------------------------------------- !
!                                                                              !

Q0 = 0.0_rk
IF (BB.GE.0.25) Q0 = (2.*SQRT(BB)-1.)**2

! ---------------------------------------------------------------------------- !
!                                                                              !

IF (BB.LT.1.0) THEN
   QB = Q0 - BB*(Q0 - exp((Q0-1.)/BB))/(BB-exp((Q0-1.)/BB))
ELSE
    QB = 1.0_rk
END IF

END SUBROUTINE CMPQB_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SINPUT_ARD_OPENACC (F,SL,SPOS,FL,USTAR,UDIR,Z0,ROAIRN,WSTAR,LLWS)

! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY :
!       J-R BIDLOT NOVEMBER 1995
!       J-R BIDLOT FEBRUARY 1996-97
!       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
!       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
!                                  USING NEW STRESS AND ROUGHNESS.
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
!                                 DENSITY AND STABILITY-DEPENDENT
!                                 WIND GUSTINESS
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE
!                                      RUNNING FASTER THAN THE WIND.
!       J-R BIDLOT JANUARY 2013: SHALLOW WATER FORMULATION.


!       L. AOUF    March 2011 : USE OF NEW DISSIPATION DEVELOPED BY ARDHUIN ET AL.2010

!       JEAN BIDLOT : ADAPTED TO ECWAM AND WAM.

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!     METHOD.                                                                  !
!     ------- 

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       WSIGSTAR.

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991
!       ARDHUIN ET AL. 2010

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)    :: F (:, :, :)    !! SPECTRUM.
REAL(real_kind),    INTENT(OUT)   :: SL(:, :, :)    !! TOTAL SOURCE FUNCTION ARRAY
REAL(real_kind),    INTENT(OUT)   :: SPOS(:, :, :)  !! POSITIVE SOURCE FUNCTION ARRAY
REAL(real_kind),    INTENT(OUT)   :: FL(:, :, :)    !! DIAGONAL MATRIX OF FUNCTIONAL
                                                    !! DERIVATIVE
REAL(real_kind),    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL(real_kind),    INTENT(IN)    :: UDIR (:)       !! WIND DIRECTION.
REAL(real_kind),    INTENT(IN)    :: Z0   (:)       !! ROUGHNESS LENGTH.
REAL(real_kind),    INTENT(IN)    :: ROAIRN(:)      !! AIR DENSITY
REAL(real_kind),    INTENT(IN)    :: WSTAR(:)       !! FREE CONVECTION VELOCITY SCLALE 
LOGICAL, INTENT(OUT)   :: LLWS(:, :, :)  !! TRUE WHERE SINPUT IS POSITIVE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER :: IJ,K,M,IND,IGST,I

INTEGER , PARAMETER :: NGST=2

!     TEST471 (modified SWELLF):
REAL(real_kind), PARAMETER :: SWELLF = 0.66 ! controls the turbulent swell dissipation
REAL(real_kind), PARAMETER :: SWELLF2 = -0.018
REAL(real_kind), PARAMETER :: SWELLF3 = 0.022
REAL(real_kind), PARAMETER :: SWELLF4 = 1.5E05
REAL(real_kind), PARAMETER :: SWELLF5 = 1.2  ! controls the viscous swell dissipation
REAL(real_kind), PARAMETER :: SWELLF6 = 1.0
REAL(real_kind), PARAMETER :: SWELLF7 = 3.6E05
REAL(real_kind), PARAMETER :: Z0RAT = 0.04
REAL(real_kind), PARAMETER :: Z0TUBMAX = 0.0005

REAL(real_kind), PARAMETER :: ABMIN = 0.3
REAL(real_kind), PARAMETER :: ABMAX = 8.0 


REAL(real_kind) :: ROG
REAL(real_kind) :: AVG_GST, ABS_TAUWSHELTER 
REAL(real_kind) :: CONST1
REAL(real_kind) :: ZLOG, ZLOG2X
REAL(real_kind) :: XI,X,DELI1
REAL(real_kind) :: FU,FUD,NU_AIR,SMOOTH, HFTSWELLF6,Z0TUB
REAL(real_kind) :: FAC_NU_AIR,FACM1_NU_AIR
REAL(real_kind) :: DELABM1
REAL(real_kind) :: TAUPX,TAUPY
REAL(real_kind) :: DSTAB2
REAL(real_kind), DIMENSION(SIZE(F,3)) :: CONST, SIG, SIGM1, SIG2, COEF, COEF5, DFIM_SIG2
REAL(real_kind), DIMENSION(SIZE(F,1)) :: PVISC, PTURB
REAL(real_kind) :: CONSTF
REAL(real_kind) :: Z0VIS, FWW
REAL(real_kind) :: ZCN
REAL(real_kind) :: TEMP, ZORB
REAL(real_kind), DIMENSION(SIZE(F,1)) :: SIG_N, UORBT, AORB, RE, RE_C
REAL(real_kind) :: CNSN
REAL(real_kind) :: FLP_AVG, SLP_AVG
REAL(real_kind) :: ROGOROAIR, ROAIRN_PVISC
REAL(real_kind) :: XSTRESS, YSTRESS
REAL(real_kind) :: USTP, USDIRP, UCN, UCNZALPD, FLP,SLP
REAL(real_kind) :: USG2
REAL(real_kind), DIMENSION(SIZE(F,1),NGST) :: TAUX, TAUY
REAL(real_kind) :: CM, XK, SH
REAL(real_kind), DIMENSION(SIZE(F,1)) :: TEMP2
REAL(real_kind) :: DSTAB1, TEMP1, COSLP
REAL(real_kind) :: UFAC,DSTAB


      ROG = rho_sea*G
      AVG_GST = 1.0/NGST
      CONST1  = BETAMAX/XKAPPA**2 /rho_sea
      NU_AIR = RNUAIR
      FAC_NU_AIR= RNUAIRM
      FACM1_NU_AIR=4.0/NU_AIR

      FU=ABS(SWELLF3)
      FUD=SWELLF2
      DELABM1= REAL(IAB)/(ABMAX-ABMIN)

      ABS_TAUWSHELTER=ABS(TAUWSHELTER)

!     ESTIMATE THE STANDARD DEVIATION OF GUSTINESS.
      CALL WSIGSTAR_OPENACC (USTAR, Z0, WSTAR, SIG_N)

! ----------------------------------------------------------------------

      DO M = 1,SIZE(F,3)
        SIG(M) = PI2*FRE(M)
        SIGM1(M) = 1.0/SIG(M)
        SIG2(M) = SIG(M)**2
        DFIM_SIG2(M) =DFIM(M)*SIG2(M)
        CONST(M)=SIG(M)*CONST1
        COEF(M) =-SWELLF*16.*SIG2(M)/(G*rho_sea)
        COEF5(M)=-SWELLF5*2.*SQRT(2.*NU_AIR*SIG(M))/rho_sea
      ENDDO
      DO IJ = 1,SIZE(F,1)
        UORBT(IJ) = TINY(1.0) 
        AORB(IJ) = TINY(1.0) 
        DO M = 1,SIZE(F,3)
          TEMP = F(IJ,1,M)   !K=1
          DO K = 2,SIZE(F,2)
            TEMP = TEMP+F(IJ,K,M)
          ENDDO
          UORBT(IJ) = UORBT(IJ)+DFIM_SIG2(M)*TEMP
          AORB(IJ) = AORB(IJ)+DFIM(M)*TEMP
        ENDDO
      ENDDO

      DO IJ = 1,SIZE(F,1)
!!!!
        UORBT(IJ) = 2.0*SQRT(UORBT(IJ))  ! this is the significant orbital amplitude
        AORB(IJ)  = 2.0*SQRT(AORB(IJ))   ! this 1/2 Hs
        RE(IJ)    = FACM1_NU_AIR*UORBT(IJ)*AORB(IJ) ! this is the Reynolds number 
        Z0VIS = FAC_NU_AIR/MAX(USTAR(IJ),0.0001_rk)
        Z0TUB = Z0RAT*MIN(Z0TUBMAX,Z0(IJ))
        !Z0NOZ = MAX(Z0VIS,Z0TUB)
        ZORB  = AORB(IJ)/MAX(Z0VIS,Z0TUB)

! conpute fww
        XI=(LOG10(MAX(ZORB,3.0_rk))-ABMIN)*DELABM1
        IND  = MIN (IAB-1, INT(XI))
        DELI1= MIN (1.0_rk ,XI-FLOAT(IND))
        !DELI2= 1.0 - DELI1
        FWW = SWELLFT(IND)*(1.0 - DELI1)+SWELLFT(IND+1)*DELI1
        TEMP2(IJ) = FWW*UORBT(IJ)
      ENDDO

! Define the critical Reynolds number
      IF( SWELLF6 .EQ. 1.0) THEN
        DO IJ = 1,SIZE(F,1)
          RE_C(IJ) = SWELLF4
        ENDDO
      ELSE
        HFTSWELLF6=1.0-SWELLF6
        DO IJ = 1,SIZE(F,1)
          RE_C(IJ) = SWELLF4*(2.0/AORB(IJ))**HFTSWELLF6
        ENDDO
      ENDIF

! Swell damping weight between viscous and turbulent boundary layer
      IF (SWELLF7 .GT. 0.0) THEN
        DO IJ = 1,SIZE(F,1)
          SMOOTH=0.5*TANH((RE(IJ)-RE_C(IJ))/SWELLF7)
          PTURB(IJ)=0.5+SMOOTH
          PVISC(IJ)=0.5-SMOOTH
        ENDDO
      ELSE
        DO IJ = 1,SIZE(F,1)
          IF (RE(IJ).LE.RE_C(IJ)) THEN
            PTURB(IJ)=0.0_rk
            PVISC(IJ)=0.5_rk
          ELSE
            PTURB(IJ)=0.5_rk
            PVISC(IJ)=0.0_rk
          ENDIF
        ENDDO
      ENDIF


! Initialisation

      DO IGST=1,NGST
        DO IJ = 1,SIZE(F,1)
          !XSTRESS(IJ,IGST)=0.0
          !YSTRESS(IJ,IGST)=0.0
          
          !Wind gustiness
          IF (IGST.EQ.1) THEN
            USG2=(USTAR(IJ)*(1.0+SIG_N(IJ)))**2
          ELSE
            USG2=(USTAR(IJ)*(1.0-SIG_N(IJ)))**2
          ENDIF
          TAUX(IJ,IGST)=USG2*SIN(UDIR(IJ))
          TAUY(IJ,IGST)=USG2*COS(UDIR(IJ))
        ENDDO
      ENDDO

      DO M = 1,SIZE(F,3)
        DO K = 1,SIZE(F,2)
          DO IJ = 1,SIZE(F,1)
            LLWS(IJ,K,M)=.FALSE.
          ENDDO
        ENDDO
      ENDDO

!*    2. MAIN LOOP OVER FREQUENCIES.
!        ---------------------------

      DO IGST = 1,NGST    !NGST can only be 2
        DO IJ = 1,SIZE(F,1)
      
          ROGOROAIR = ROG/MAX(ROAIRN(IJ),1.0_rk)
          ROAIRN_PVISC = ROAIRN(IJ)*PVISC(IJ)

          XSTRESS = 0.0_rk
          YSTRESS = 0.0_rk

!there is iteration dependence on XSTRESS an YSTRESS
          DO M = 1,SIZE(F,3)

!           INVERSE OF PHASE VELOCITIES AND WAVE NUMBER.
            IF (SHALLOW_RUN) THEN
                  I = cellIdxTab(IJ)
                  XK = kwave(I,M)
                  CM = XK*SIGM1(M)
!!                SH(IJ,M) = SIG2(M)/(G*XK(IJ,M)) ! tanh(kh)
                  SH = 1.0
            ELSE
                  CM = SIG(M)/G
                  XK = SIG2(M)/G
                  SH = 1.0_rk
            ENDIF

            TAUPX=TAUX(IJ,IGST)-ABS_TAUWSHELTER*XSTRESS
            TAUPY=TAUY(IJ,IGST)-ABS_TAUWSHELTER*YSTRESS
            USTP=(TAUPX**2+TAUPY**2)**0.25
            USDIRP=ATAN2(TAUPX,TAUPY)

            CONSTF = ROGOROAIR*CM*DFIM(M)


         ! DO IGST=1,NGST
         !   DO K = 1,SIZE(F,2)
         !     COSLP(IJ,K,IGST) = COS(TH(K)-USDIRP(IJ,IGST))
         !   ENDDO
         ! ENDDO

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

            UCN = USTP*CM
            UCNZALPD = XKAPPA/(UCN + ZALP)
            ZCN = LOG(XK*Z0(IJ))
            CNSN = CONST(M)*SH*ROAIRN(IJ)

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

            DSTAB1 = COEF5(M)*ROAIRN_PVISC*XK
            TEMP1 = COEF(M)*ROAIRN(IJ)

            DO K = 1,SIZE(F,2)
              COSLP = COS(THETA(K)-USDIRP)
              IF (COSLP.GT.0.01) THEN
                X    = COSLP*UCN
                ZLOG = ZCN + UCNZALPD/COSLP
                IF (ZLOG.LT.0.0) THEN
                  ZLOG2X=ZLOG*ZLOG*X
                  UFAC = EXP(ZLOG)*ZLOG2X*ZLOG2X
                  LLWS(IJ,K,M)=.TRUE.
                ELSE
                  UFAC=0.0
                ENDIF
              ELSE
                UFAC=0.0
              ENDIF

!       SWELL DAMPING:

              DSTAB2 = TEMP1*(TEMP2(IJ)+(FU+FUD*( COS(THETA(K)-USDIRP) ))*USTP)
              DSTAB = DSTAB1+PTURB(IJ)*DSTAB2



!*    2.2 UPDATE THE SHELTERING STRESS,
!         AND THEN ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ---------------------------------------------------------

            !CONST11=CONSTF*SINTH(K)
            !CONST22=CONSTF*COSTH(K)

              ! SLP: only the positive contributions
              SLP = CNSN*UFAC
              FLP = SLP+DSTAB

              SLP = SLP*F(IJ,K,M)
              XSTRESS=XSTRESS+SLP*CONSTF*SINTHETA(K)
              YSTRESS=YSTRESS+SLP*CONSTF*COSTHETA(K)

           
              IF (IGST==1) THEN
                SPOS(IJ,K,M) = SLP
                FL(IJ,K,M) = FLP
              ELSEIF (IGST==2) THEN
                SPOS(IJ,K,M) = AVG_GST*(SPOS(IJ,K,M)+SLP)
                FL(IJ,K,M) = AVG_GST*(FL(IJ,K,M)+FLP)
                SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
              ELSE
                if(.true.) print * , 'NGST only can be 2 for OPENACC VERSION'
                !CALL ABORT1
              ENDIF



            ENDDO    ! K 

          ENDDO ! END LOOP OVER FREQUENCIES
        ENDDO   ! END LOOP OVER IJ
      ENDDO   !  IGST



END SUBROUTINE SINPUT_ARD_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SNONLIN (F, SL, FL, DEPTH, AKMEAN)

! ---------------------------------------------------------------------------- !
!                                                                              !
! *** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS               !
! ***             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND              !
! ***             ADDITION TO CORRESPONDING NET EXPRESSIONS.                   !
!                                                                              !

REAL(real_kind), INTENT(IN)    :: F (:,:,:)    !! SPECTRA.
REAL(real_kind), INTENT(INOUT) :: SL(:,:,:)    !! TOTAL SOURCE FUNCTION ARRAY.
REAL(real_kind), INTENT(INOUT) :: FL(:,:,:)    !! DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
REAL(real_kind), INTENT(IN)    :: DEPTH (:)    !! WATER DEPTH.
REAL(real_kind), INTENT(IN)    :: AKMEAN (:)   !! MEAN WAVE NUMBER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: MC, MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1, KH, K, K1, K2, K11, K21
INTEGER :: IJ
INTEGER :: MFR1STFR, MFRLSTFR

REAL(real_kind)    :: FTAIL, FKLAMP, FKLAMP1, GW1, GW2, GW3, GW4
REAL(real_kind)    :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22
REAL(real_kind)    :: FKLAMM, FKLAMM1, GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
REAL(real_kind)    :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
REAL(real_kind)    :: SAP, SAM, FIJ, FAD1, FAD2, FCEN

REAL(real_kind), DIMENSION(SIZE(F,1)) :: FTEMP
REAL(real_kind), DIMENSION(SIZE(F,1)) :: ENHFR
REAL(real_kind), DIMENSION(SIZE(F,1)) :: AD
REAL(real_kind), DIMENSION(SIZE(F,1)) :: DELAD
REAL(real_kind), DIMENSION(SIZE(F,1)) :: DELAP
REAL(real_kind), DIMENSION(SIZE(F,1)) :: DELAM

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. SHALLOW WATER INITIALISATION.                                         !
!        -----------------------------                                         !

   IF (ISNONLIN.EQ.0) THEN
      DO IJ = 1,SIZE(F,1)
         ENHFR(IJ) = MAX(0.75*DEPTH(IJ)*AKMEAN(IJ) , 0.5_rk)
         ENHFR(IJ) = 1. + (5.5/ENHFR(IJ)) * (1.-.833*ENHFR(IJ))               &
&                    * EXP(-1.25*ENHFR(IJ))
      END DO
      DO MC=1,MLSTHG
         ENH(:,MC) = ENHFR(:)
      END DO
   END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. FREQUENCY LOOP.                                                       !
!        ---------------                                                       !

MFR1STFR = -MFRSTLW+1
MFRLSTFR = ML-KFRH+MFR1STFR

FRE: DO MC = 1,MLSTHG
   MP  = IKP (MC)
   MP1 = IKP1(MC)
   MM  = IKM (MC)
   MM1 = IKM1(MC)
   IC  = INLCOEF(1,MC)
   IP  = INLCOEF(2,MC)
   IP1 = INLCOEF(3,MC)
   IM  = INLCOEF(4,MC)
   IM1 = INLCOEF(5,MC)

   FTAIL  = RNLCOEF(1,MC)


   FKLAMP  = FKLAP(MC)
   FKLAMP1 = FKLAP1(MC)
   GW1 = RNLCOEF(2,MC)
   GW2 = RNLCOEF(3,MC)
   GW3 = RNLCOEF(4,MC)
   GW4 = RNLCOEF(5,MC)
   FKLAMPA = RNLCOEF(6,MC)
   FKLAMPB = RNLCOEF(7,MC)
   FKLAMP2 = RNLCOEF(8,MC)
   FKLAMP1 = RNLCOEF(9,MC)
   FKLAPA2 = RNLCOEF(10,MC)
   FKLAPB2 = RNLCOEF(11,MC)
   FKLAP12 = RNLCOEF(12,MC)
   FKLAP22 = RNLCOEF(13,MC)

   FKLAMM  = FKLAM(MC)
   FKLAMM1 = FKLAM1(MC)
   GW5 = RNLCOEF(14,MC)
   GW6 = RNLCOEF(15,MC)
   GW7 = RNLCOEF(16,MC)
   GW8 = RNLCOEF(17,MC)
   FKLAMMA = RNLCOEF(18,MC)
   FKLAMMB = RNLCOEF(19,MC)
   FKLAMM2 = RNLCOEF(20,MC)
   FKLAMM1 = RNLCOEF(21,MC)
   FKLAMA2 = RNLCOEF(22,MC)
   FKLAMB2 = RNLCOEF(23,MC)
   FKLAM12 = RNLCOEF(24,MC)
   FKLAM22 = RNLCOEF(25,MC)

   FTEMP(:) = AF11(MC)*ENH(:,MC)

   IF (MC.GT.MFR1STFR .AND. MC.LT.MFRLSTFR ) THEN
!       MC is within the fully resolved spectral domain

      DO KH=1,2
         DO K=1,KL
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

!*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
!*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!             ----------------------------------------------

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )                     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )                     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
!!!! not needed ftail always=1.                FIJ = F(IJ,K  ,IC )*FTAIL
               FIJ = F(IJ,K  ,IC )
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FCEN
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FCEN
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FCEN
            ENDDO

            SL(:,K  ,MC ) = SL(:,K  ,MC ) - 2.*AD(:)
            FL(:,K  ,MC ) = FL(:,K  ,MC ) - 2.*DELAD(:)
            SL(:,K2 ,MM ) = SL(:,K2 ,MM ) + AD(:)*FKLAMM1
            FL(:,K2 ,MM ) = FL(:,K2 ,MM ) + DELAM(:)*FKLAM12
            SL(:,K21,MM ) = SL(:,K21,MM ) + AD(:)*FKLAMM2
            FL(:,K21,MM ) = FL(:,K21,MM ) + DELAM(:)*FKLAM22
            SL(:,K2 ,MM1) = SL(:,K2 ,MM1) + AD(:)*FKLAMMA
            FL(:,K2 ,MM1) = FL(:,K2 ,MM1) + DELAM(:)*FKLAMA2
            SL(:,K21,MM1) = SL(:,K21,MM1) + AD(:)*FKLAMMB
            FL(:,K21,MM1) = FL(:,K21,MM1) + DELAM(:)*FKLAMB2
            SL(:,K1 ,MP ) = SL(:,K1 ,MP ) + AD(:)*FKLAMP1
            FL(:,K1 ,MP ) = FL(:,K1 ,MP ) + DELAP(:)*FKLAP12
            SL(:,K11,MP ) = SL(:,K11,MP ) + AD(:)*FKLAMP2
            FL(:,K11,MP ) = FL(:,K11,MP ) + DELAP(:)*FKLAP22
            SL(:,K1 ,MP1) = SL(:,K1 ,MP1) + AD(:)*FKLAMPA
            FL(:,K1 ,MP1) = FL(:,K1 ,MP1) + DELAP(:)*FKLAPA2
            SL(:,K11,MP1) = SL(:,K11,MP1) + AD(:)*FKLAMPB
            FL(:,K11,MP1) = FL(:,K11,MP1) + DELAP(:)*FKLAPB2
         ENDDO
      ENDDO

   ELSEIF (MC.GE.MFRLSTFR ) THEN
      DO KH=1,2
         DO K=1,KL
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )                     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )                     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
               FIJ = F(IJ,K  ,IC )*FTAIL
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FCEN
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FCEN
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FCEN
            ENDDO

            SL(:,K2 ,MM ) = SL(:,K2 ,MM ) + AD(:)*FKLAMM1
            FL(:,K2 ,MM ) = FL(:,K2 ,MM ) + DELAM(:)*FKLAM12
            SL(:,K21,MM ) = SL(:,K21,MM ) + AD(:)*FKLAMM2
            FL(:,K21,MM ) = FL(:,K21,MM ) + DELAM(:)*FKLAM22

            IF (MM1.LE.ML) THEN
               SL(:,K2 ,MM1) = SL(:,K2 ,MM1) + AD(:)*FKLAMMA
               FL(:,K2 ,MM1) = FL(:,K2 ,MM1) + DELAM(:)*FKLAMA2
               SL(:,K21,MM1) = SL(:,K21,MM1) + AD(:)*FKLAMMB
               FL(:,K21,MM1) = FL(:,K21,MM1) + DELAM(:)*FKLAMB2

               IF (MC .LE.ML) THEN
                  SL(:,K  ,MC ) = SL(:,K  ,MC ) - 2.*AD(:)
                  FL(:,K  ,MC ) = FL(:,K  ,MC ) - 2.*DELAD(:)

                  IF (MP .LE.ML) THEN
                     SL(:,K1 ,MP ) = SL(:,K1 ,MP ) + AD(:)*FKLAMP1
                     FL(:,K1 ,MP ) = FL(:,K1 ,MP ) + DELAP(:)*FKLAP12
                     SL(:,K11,MP ) = SL(:,K11,MP ) + AD(:)*FKLAMP2
                     FL(:,K11,MP ) = FL(:,K11,MP ) + DELAP(:)*FKLAP22

                     IF (MP1.LE.ML) THEN
                        SL(:,K1 ,MP1) = SL(:,K1 ,MP1) + AD(:)*FKLAMPA
                        FL(:,K1 ,MP1) = FL(:,K1 ,MP1) + DELAP(:)*FKLAPA2
                        SL(:,K11,MP1) = SL(:,K11,MP1) + AD(:)*FKLAMPB
                        FL(:,K11,MP1) = FL(:,K11,MP1) + DELAP(:)*FKLAPB2
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   ELSE

      DO KH=1,2
         DO K=1,KL
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
               FIJ = F(IJ,K  ,IC )*FTAIL
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FCEN
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FCEN
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FCEN
            ENDDO

            IF (MM1.GE.1) THEN
               SL(:,K2 ,MM1) = SL(:,K2 ,MM1) + AD(:)*FKLAMMA
               FL(:,K2 ,MM1) = FL(:,K2 ,MM1) + DELAM(:)*FKLAMA2
               SL(:,K21,MM1) = SL(:,K21,MM1) + AD(:)*FKLAMMB
               FL(:,K21,MM1) = FL(:,K21,MM1) + DELAM(:)*FKLAMB2
            ENDIF

            SL(:,K  ,MC ) = SL(:,K  ,MC ) - 2.*AD(:)
            FL(:,K  ,MC ) = FL(:,K  ,MC ) - 2.*DELAD(:)
            SL(:,K1 ,MP ) = SL(:,K1 ,MP ) + AD(:)*FKLAMP1
            FL(:,K1 ,MP ) = FL(:,K1 ,MP ) + DELAP(:)*FKLAP12
            SL(:,K11,MP ) = SL(:,K11,MP ) + AD(:)*FKLAMP2
            FL(:,K11,MP ) = FL(:,K11,MP ) + DELAP(:)*FKLAP22
            SL(:,K1 ,MP1) = SL(:,K1 ,MP1) + AD(:)*FKLAMPA
            FL(:,K1 ,MP1) = FL(:,K1 ,MP1) + DELAP(:)*FKLAPA2
            SL(:,K11,MP1) = SL(:,K11,MP1) + AD(:)*FKLAMPB
            FL(:,K11,MP1) = FL(:,K11,MP1) + DELAP(:)*FKLAPB2
         ENDDO
      ENDDO

   ENDIF
END DO FRE                  !! BRANCH BACK FOR NEXT FREQUENCY.

END SUBROUTINE SNONLIN

SUBROUTINE SNONLIN_OPENACC (F, SL, FL, DEPTH, AKMEAN)

! ---------------------------------------------------------------------------- !
!                                                                              !
! *** *SNONLIN* - COMPUTATION OF NONLINEAR TRANSFER RATE AND ITS               !
! ***             FUNCTIONAL DERIVATIVE (DIAGONAL TERMS ONLY) AND              !
! ***             ADDITION TO CORRESPONDING NET EXPRESSIONS.                   !
!                                                                              !
!     S.D. HASSELMANN.  MPI                                                    !
!                                                                              !
!     G. KOMEN, P. JANSSEN   KNMI             MODIFIED TO SHALLOW WATER        !
!     H. GUENTHER, L. ZAMBRESKY               OPTIMIZED                        !
!     H. GUENTHER       GKSS/ECMWF  JUNE 1991 INTERACTIONS BETWEEN DIAG-       !
!                                             AND PROGNOSTIC PART.             !
!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90                            !
!     E. MYKLEBUST        FEBRUARY 2005       OPTIMIZATION                     !
!     P. JANSSEN  ECMWF  JUNE 2005       IMPROVED SCALING IN SHALLOW WATER     !
!     J. BIDLOT   ECMWF  AUGUST 2006     KEEP THE OLD FORMULATION              !
!                                        UNDER A SWITCH (ISNONLIN = 0 for OLD  !
!                                                                 = 1 for NEW  !
!                                        BE AWARE THAT THE OLD FORMULATION     !
!                                        REQUIRES THE MEAN WAVE NUMBER AKMEAN. !
!     J. BIDLOT   ECMWF  JANUARY 2012    ADD EXTENSION TO LOW FREQUENCIES      !
!                                        OPTIMISATION FOR IBM.                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       SEE ABOVE.                                                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind), INTENT(IN)    :: F (:,:,:)    !! SPECTRA.
REAL(real_kind), INTENT(INOUT) :: SL(:,:,:)    !! TOTAL SOURCE FUNCTION ARRAY.
REAL(real_kind), INTENT(INOUT) :: FL(:,:,:)    !! DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
REAL(real_kind), INTENT(IN)    :: DEPTH (:)    !! WATER DEPTH.
REAL(real_kind), INTENT(IN)    :: AKMEAN (:)   !! MEAN WAVE NUMBER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: MC, MP, MP1, MM, MM1, IC, IP, IP1, IM, IM1, KH, K, K1, K2, K11, K21
INTEGER :: IJ
INTEGER :: MFR1STFR, MFRLSTFR

REAL(real_kind)    :: FTAIL, FKLAMP, FKLAMP1, GW1, GW2, GW3, GW4
REAL(real_kind)    :: FKLAMPA, FKLAMPB, FKLAMP2, FKLAPA2, FKLAPB2, FKLAP12, FKLAP22
REAL(real_kind)    :: FKLAMM, FKLAMM1, GW5, GW6, GW7, GW8, FKLAMMA, FKLAMMB, FKLAMM2
REAL(real_kind)    :: FKLAMA2, FKLAMB2, FKLAM12, FKLAM22
REAL(real_kind)    :: SAP, SAM, FIJ, FAD1, FAD2, FCEN

REAL(real_kind), DIMENSION(SIZE(F,1)) :: FTEMP
REAL(real_kind), DIMENSION(SIZE(F,1)) :: ENHFR
REAL(real_kind), DIMENSION(SIZE(F,1)) :: AD
REAL(real_kind), DIMENSION(SIZE(F,1)) :: DELAD
REAL(real_kind), DIMENSION(SIZE(F,1)) :: DELAP
REAL(real_kind), DIMENSION(SIZE(F,1)) :: DELAM


!OPENACC ---

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. SHALLOW WATER INITIALISATION.                                         !
!        -----------------------------                                         !


IF (SHALLOW_RUN) THEN
   IF (ISNONLIN.EQ.0) THEN
      DO IJ = 1,SIZE(DEPTH)
         ENHFR(IJ) = MAX(0.75*DEPTH(IJ)*AKMEAN(IJ) , 0.5_rk)
         ENHFR(IJ) = 1. + (5.5/ENHFR(IJ)) * (1.-.833*ENHFR(IJ))               &
&                    * EXP(-1.25*ENHFR(IJ))
      END DO
      DO MC=1,MLSTHG
         DO IJ = 1,SIZE(ENH,1)
            ENH(IJ,MC) = ENHFR(IJ)
         END DO
      END DO
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. FREQUENCY LOOP.                                                       !
!        ---------------                                                       !

MFR1STFR = -MFRSTLW+1
MFRLSTFR = ML-KFRH+MFR1STFR
#ifdef MPI_DEBUG_OUTPUT
DO MC=1,MLSTHG
  IF(IKM(MC)<-3) THEN
    WRITE(*,*)"out of range,IKM(MC)=",IKM(MC),",MC=",MC,",mpirank=",mpi_rank
  END IF
  IF(IKM1(MC)<-3) THEN
    WRITE(*,*)"out of range,IKM1(MC)=",IKM1(MC),",MC=",MC,",mpirank=",mpi_rank
  END IF
END DO
#endif
DO MC = 1,MLSTHG
   MP  = IKP (MC)
   MP1 = IKP1(MC)
   MM  = IKM (MC)
   MM1 = IKM1(MC)
   IC  = INLCOEF(1,MC)
   IP  = INLCOEF(2,MC)
   IP1 = INLCOEF(3,MC)
   IM  = INLCOEF(4,MC)
   IM1 = INLCOEF(5,MC)

   FTAIL  = RNLCOEF(1,MC)

   FKLAMP  = FKLAP(MC)
   FKLAMP1 = FKLAP1(MC)
   GW1 = RNLCOEF(2,MC)
   GW2 = RNLCOEF(3,MC)
   GW3 = RNLCOEF(4,MC)
   GW4 = RNLCOEF(5,MC)
   FKLAMPA = RNLCOEF(6,MC)
   FKLAMPB = RNLCOEF(7,MC)
   FKLAMP2 = RNLCOEF(8,MC)
   FKLAMP1 = RNLCOEF(9,MC)
   FKLAPA2 = RNLCOEF(10,MC)
   FKLAPB2 = RNLCOEF(11,MC)
   FKLAP12 = RNLCOEF(12,MC)
   FKLAP22 = RNLCOEF(13,MC)

   FKLAMM  = FKLAM(MC)
   FKLAMM1 = FKLAM1(MC)
   GW5 = RNLCOEF(14,MC)
   GW6 = RNLCOEF(15,MC)
   GW7 = RNLCOEF(16,MC)
   GW8 = RNLCOEF(17,MC)
   FKLAMMA = RNLCOEF(18,MC)
   FKLAMMB = RNLCOEF(19,MC)
   FKLAMM2 = RNLCOEF(20,MC)
   FKLAMM1 = RNLCOEF(21,MC)
   FKLAMA2 = RNLCOEF(22,MC)
   FKLAMB2 = RNLCOEF(23,MC)
   FKLAM12 = RNLCOEF(24,MC)
   FKLAM22 = RNLCOEF(25,MC)

   IF (SHALLOW_RUN) THEN
      DO IJ = 1,SIZE(FTEMP)
         FTEMP(IJ) = AF11(MC)*ENH(IJ,MC)
      END DO
   ELSE
      DO IJ = 1,SIZE(FTEMP)
         FTEMP(IJ) = AF11(MC)
      END DO
   ENDIF


   IF (MC.GT.MFR1STFR .AND. MC.LT.MFRLSTFR ) THEN
!       MC is within the fully resolved spectral domain

      DO KH=1,2
         DO K=1,SIZE(F,2)
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

!*    2.1.1.1 LOOP OVER GRIDPOINTS.. NONLINEAR TRANSFER AND
!*            DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!             ----------------------------------------------

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )                     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )                     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
!!!! not needed ftail always=1.                FIJ = F(IJ,K  ,IC )*FTAIL
               FIJ = F(IJ,K  ,IC )
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               !FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FTEMP(IJ)*FIJ
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FTEMP(IJ)*FIJ
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FTEMP(IJ)*FIJ
            ENDDO

            !add YUAN, pay attention, SL?FL?,dimension K & M is not indep
            !considering using openacc cache clause
            DO IJ = 1,SIZE(AD)
               SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD(IJ)
               FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD(IJ)
               SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD(IJ)*FKLAMM1
               FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM(IJ)*FKLAM12
               SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD(IJ)*FKLAMM2
               FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM(IJ)*FKLAM22
               SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD(IJ)*FKLAMMA
               FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM(IJ)*FKLAMA2
               SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD(IJ)*FKLAMMB
               FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM(IJ)*FKLAMB2
               SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
               FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP ) + DELAP(IJ)*FKLAP12
               SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
               FL(IJ,K11,MP ) = FL(IJ,K11,MP ) + DELAP(IJ)*FKLAP22
               SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD(IJ)*FKLAMPA
               FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1) + DELAP(IJ)*FKLAPA2
               SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD(IJ)*FKLAMPB
               FL(IJ,K11,MP1) = FL(IJ,K11,MP1) + DELAP(IJ)*FKLAPB2
            END DO

         ENDDO
      ENDDO

   ELSEIF (MC.GE.MFRLSTFR ) THEN
      DO KH=1,2
         DO K=1,KL
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )                     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )                     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
               FIJ = F(IJ,K  ,IC )*FTAIL
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               !FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FTEMP(IJ)*FIJ
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FTEMP(IJ)*FIJ
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FTEMP(IJ)*FIJ
            ENDDO

            DO IJ = 1,SIZE(AD)
               SL(IJ,K2 ,MM ) = SL(IJ,K2 ,MM ) + AD(IJ)*FKLAMM1
               FL(IJ,K2 ,MM ) = FL(IJ,K2 ,MM ) + DELAM(IJ)*FKLAM12
               SL(IJ,K21,MM ) = SL(IJ,K21,MM ) + AD(IJ)*FKLAMM2
               FL(IJ,K21,MM ) = FL(IJ,K21,MM ) + DELAM(IJ)*FKLAM22
            END DO

#ifdef MPI_DEBUG_OUTPUT
                  IF(K2<0 .OR. K2>nDir) THEN
                    WRITE(*,*) "K2=",K2," is out of range"
                  END IF
                  IF(K21<0 .OR. K21>nDir) THEN
                    WRITE(*,*) "K21=",K21," is out of range"
                  END IF
                  IF(MM1<0 .OR. MM1>nFre) THEN
                    WRITE(*,*) "MM1=",MM1," is out of range,mpirank=",mpi_rank
                  END IF
#endif

            IF (MM1.LE.ML) THEN
#ifdef MPI_DEBUG_OUTPUT
               IF (SIZE(AD)>1000000 .OR. SIZE(AD)<=0) THEN
                 WRITE(*,*) "SIZE(AD)=",SIZE(AD),",mpi_rank=",mpi_rank
               ENDIF
#endif
               DO IJ = 1, SIZE(AD)
#ifdef MPI_DEBUG_OUTPUT
                  IF(K2<0 .OR. K2>nDir) THEN
                    WRITE(*,*) "K2=",K2," is out of range"
                  END IF
                  IF(K21<0 .OR. K21>nDir) THEN
                    WRITE(*,*) "K21=",K21," is out of range"
                  END IF
                  IF(MM1<0 .OR. MM1>nFre) THEN
                    WRITE(*,*) "MM1=",MM1," is out of range,mpirank=",mpi_rank
                  END IF
                  IF(IJ>nCells) THEN
                    WRITE(*,*) "IJ=",IJ," is out of range"
                  END IF
#endif
                  SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD(IJ)*FKLAMMA
                  FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM(IJ)*FKLAMA2
                  SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD(IJ)*FKLAMMB
                  FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM(IJ)*FKLAMB2
               END DO

               IF (MC .LE.ML) THEN
                  DO IJ = 1,SIZE(AD)
                     SL(IJ,K,MC ) = SL(IJ,K  ,MC ) - 2.*AD(IJ)
                     FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD(IJ)
                  END DO 

                  IF (MP .LE.ML) THEN
                     DO IJ = 1,SIZE(SL,1)
                        SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
                        FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP ) + DELAP(IJ)*FKLAP12
                        SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
                        FL(IJ,K11,MP ) = FL(IJ,K11,MP ) + DELAP(IJ)*FKLAP22
                     END DO

                     IF (MP1.LE.ML) THEN
                        DO IJ = 1,SIZE(AD)
                           SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD(IJ)*FKLAMPA
                           FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1) + DELAP(IJ)*FKLAPA2
                           SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD(IJ)*FKLAMPB
                           FL(IJ,K11,MP1) = FL(IJ,K11,MP1) + DELAP(IJ)*FKLAPB2
                        END DO
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

   ELSE

      DO KH=1,2
         DO K=1,SIZE(F,2)
            K1  = K1W (K,KH)
            K2  = K2W (K,KH)
            K11 = K11W(K,KH)
            K21 = K21W(K,KH)

            DO IJ = 1,SIZE(F,1)
               SAP = GW1*F(IJ,K1 ,IP ) + GW2*F(IJ,K11,IP )     &
&                  + GW3*F(IJ,K1 ,IP1) + GW4*F(IJ,K11,IP1)
               SAM = GW5*F(IJ,K2 ,IM ) + GW6*F(IJ,K21,IM )     &
&                  + GW7*F(IJ,K2 ,IM1) + GW8*F(IJ,K21,IM1)
               FIJ = F(IJ,K  ,IC )*FTAIL
               FAD1 = FIJ*(SAP+SAM)
               FAD2 = FAD1-2.*SAP*SAM
               FAD1 = FAD1+FAD2
               !FCEN = FTEMP(IJ)*FIJ
               AD(IJ) = FAD2*FTEMP(IJ)*FIJ
               DELAD(IJ) = FAD1*FTEMP(IJ)
               DELAP(IJ) = (FIJ-2.*SAM)*DAL1*FTEMP(IJ)*FIJ
               DELAM(IJ) = (FIJ-2.*SAP)*DAL2*FTEMP(IJ)*FIJ
            ENDDO

            IF (MM1.GE.1) THEN
               DO IJ = 1,SIZE(AD)
                  SL(IJ,K2 ,MM1) = SL(IJ,K2 ,MM1) + AD(IJ)*FKLAMMA
                  FL(IJ,K2 ,MM1) = FL(IJ,K2 ,MM1) + DELAM(IJ)*FKLAMA2
                  SL(IJ,K21,MM1) = SL(IJ,K21,MM1) + AD(IJ)*FKLAMMB
                  FL(IJ,K21,MM1) = FL(IJ,K21,MM1) + DELAM(IJ)*FKLAMB2
               END DO
            ENDIF

            DO IJ = 1,SIZE(AD)
               SL(IJ,K  ,MC ) = SL(IJ,K  ,MC ) - 2.*AD(IJ)
               FL(IJ,K  ,MC ) = FL(IJ,K  ,MC ) - 2.*DELAD(IJ)
               SL(IJ,K1 ,MP ) = SL(IJ,K1 ,MP ) + AD(IJ)*FKLAMP1
               FL(IJ,K1 ,MP ) = FL(IJ,K1 ,MP ) + DELAP(IJ)*FKLAP12
               SL(IJ,K11,MP ) = SL(IJ,K11,MP ) + AD(IJ)*FKLAMP2
               FL(IJ,K11,MP ) = FL(IJ,K11,MP ) + DELAP(IJ)*FKLAP22
               SL(IJ,K1 ,MP1) = SL(IJ,K1 ,MP1) + AD(IJ)*FKLAMPA
               FL(IJ,K1 ,MP1) = FL(IJ,K1 ,MP1) + DELAP(IJ)*FKLAPA2
               SL(IJ,K11,MP1) = SL(IJ,K11,MP1) + AD(IJ)*FKLAMPB
               FL(IJ,K11,MP1) = FL(IJ,K11,MP1) + DELAP(IJ)*FKLAPB2
            END DO
         ENDDO
      ENDDO

   ENDIF
END DO                   !! BRANCH BACK FOR NEXT FREQUENCY.


END SUBROUTINE SNONLIN_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SOURCE_PHILLIPS (SL, USTAR, UDIR, DEPTH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   SOURCE_PHILLIPS - COMPUTATION OF PHILLIPS INPUT.                           !
!                                                                              !
!     ROOP LALBEHARRY          ARMN/MSC        NOVEMBER 2003
!                                                                              !
!*    PURPOSE.                                                                 !
!     ---------                                                                !
!                                                                              !
!       COMPUTE PHILLIPS INPUT SOURCE FUNCTION FACTOR.                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       SEE REFERENCE.                                                         !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       O. PHILLIPS 1966, CAVALERI AND RIZZOLI 1981, H. TOLMAN 1992            !
!                                                                              !
! ---------------------------------------------------------------------------- !
!
!*    INTERFACE VARIABLE

REAL(real_kind),   INTENT(INOUT)        :: SL(:,:,:)  !! TOTAL SOURCE FUNCTION.
REAL(real_kind),   INTENT(IN)           :: USTAR(:)   !! FRICTION VELOCITY
REAL(real_kind),   INTENT(IN)           :: UDIR(:)    !! WIND DIRECTION
REAL(real_kind),   INTENT(IN)           :: DEPTH(:)   !! DEPTH

! ---------------------------------------------------------------------------- !
!
!*    LOCAL VARIABLE

INTEGER :: K, M, IJ
 
REAL(real_kind)    :: CONST1
REAL(real_kind)    :: TEMP(1:SIZE(SL,1)), FPM(1:SIZE(SL,1)), KD
REAL(real_kind)    :: TPHOLD(1:SIZE(SL,1),1:SIZE(SL,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    0. CONSTANTS.                                                            !
!        ----------                                                            !

CONST1 = 28/(0.13*G)
DO IJ = 1,SIZE(SL,1)
   FPM(IJ) = MAX(0.0000001_rk, CONST1*USTAR(IJ))
END DO

CONST1 = (80.*16.0*XEPS**2)/(0.5*3.0*G**2)

! ---------------------------------------------------------------------------- !
!                                                                              !
!*    1. PRECALCULATED ANGULAR DEPENDENCE.                                     !
!        ---------------------------------                                     !

DO K = 1,SIZE(SL,2)
   DO IJ = 1,SIZE(SL,1)
      TPHOLD(IJ,K) = MAX(0.0_rk,COS(THETA(K)-UDIR(IJ)))
      TPHOLD(IJ,K) = (USTAR(IJ)*TPHOLD(IJ,K))**4
   END DO
END DO


! ---------------------------------------------------------------------------- !
!                                                                              !
!*    2. SHALLOW WATER                                                         !
!        -------------                                                         !

   DO M = 1,SIZE(SL,3)
      DO IJ = 1,SIZE(SL,1)
        KD = MIN(wave_numb(sigma(M), DEPTH(IJ))*DEPTH(IJ),40.0_rk)
        TEMP(IJ) = EXP(-(FRE(M)*FPM(IJ))**(-4))
        TEMP(IJ)=CONST1*TEMP(IJ)/(1.0+2.0*KD/SINH(2.*KD))
      END DO

      DO K = 1,SIZE(SL,2)
         DO IJ = 1,SIZE(SL,1)
            SL(IJ,K,M) = SL(IJ,K,M) + TEMP(IJ) * TPHOLD(IJ,K)
         END DO
      END DO
   END DO

END SUBROUTINE SOURCE_PHILLIPS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE STRESSO_OPENACC (F, SL, USTAR, UDIR, Z0, MIJ, TAUW, PHIAW)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    STRESSO - COMPUTATION OF WAVE STRESS.                                     !
!                                                                              !
!     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT.       !
!     P.A.E.M. JANSSEN      KNMI  AUGUST    1990                               !
!     J. BIDLOT             ECMWF FEBRUARY  1996-97                            !
!     H. GUENTHER   GKSS  FEBRUARY 2002       FT 90                            !
!     J. BIDLOT             ECMWF           2007  ADD MIJ                      !
!     P.A.E.M. JANSSEN     ECMWF            2011  ADD FLUX CALULATIONS         !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY AND DIRECTIONS. !
!       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE HAS TO BE STORED IN  !
!       *SL* (CALL FIRST SINPUT, THEN STRESSO, AND THEN THE REST OF THE SOURCE !
!       FUNCTIONS)                                                             !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       R SNYDER ET AL,1981.                                                   !
!       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.                  !
!       P. JANSSEN, JPO, 1985                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)    :: F(:,:,:)       !! WAVE SPECTRUM.
REAL(real_kind),    INTENT(IN)    :: SL(:,:,:)      !! INPUT SOURCE FUNCTION.
REAL(real_kind),    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY.
REAL(real_kind),    INTENT(IN)    :: UDIR(:)        !! WIND DIRECTION.
REAL(real_kind),    INTENT(IN)    :: Z0(:)          !! ROUGHNESS LENGTH.
INTEGER, INTENT(IN)    :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE
                                         !! PROGNOSTIC RANGE.
REAL(real_kind),    INTENT(OUT)   :: TAUW(:)        !! WAVE STRESS.
REAL(real_kind),    INTENT(OUT)   :: PHIAW(:)       !! ENERGY FLUX FROM WIND INTO WAVES INTEGRATED
                                         !! OVER THE FULL FREQUENCY RANGE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: K, M, KL, ML, IJ, I

REAL(real_kind) :: GM1, COSW, CONST, SINPLUS
REAL(real_kind),DIMENSION(SIZE(F,1)) :: TAU1,PHI1,XLEVTAIL
REAL(real_kind) :: SIG, SIGM1
REAL(real_kind) :: TEMPS1,TEMPS2,TEMPS3
REAL(real_kind) :: SUMT, SUMX, SUMY
REAL(real_kind) :: CM, RHOWGDFTH,CMRHOWGDFTH
REAL(real_kind),DIMENSION(SIZE(F,1)) :: XSTRESS, YSTRESS
REAL(real_kind) :: TAUHF,PHIHF
INTEGER :: MIJ_MAX

!

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. PRECOMPUTE FREQUENCY SCALING.                                         !
!        -----------------------------                                         !

KL = SIZE(F,2)
ML = SIZE(F,3)

GM1 = 1.0/G
CONST = delTheta*(PI2)**4*GM1

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE WAVE STRESS OF BLOCK.                                         !
!        -----------------------------                                         !
!                                                                              !
!     2.1 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.       !
!         --------------------------------------------------------------       !

MIJ_MAX = 0
        
DO IJ = 1,SIZE(F,1)
   MIJ_MAX = max(MIJ_MAX,MIJ(IJ))
END DO

DO IJ = 1,SIZE(F,1)
   XLEVTAIL(IJ) = 0.0_rk

   TEMPS1 = 0.0_rk
   TEMPS2 = 0.0_rk
   TEMPS3 = 0.0_rk
   
   DO M=1,MIJ_MAX !! THE INTEGRATION ONLY UP TO FR=MIJ SINCE RHOWGDFTH=0 FOR FR>MIJ
!     INVERSE OF PHASE VELOCITIES
      IF (SHALLOW_RUN) THEN
         SIGM1 = 1.0/(PI2*FRE(M))
         I = cellIdxTab(IJ)
         CM = kwave(I,M)*SIGM1
      ELSE
         SIG = PI2*FRE(M)
         CM = SIG*GM1
      ENDIF
!
      RHOWGDFTH = RHOWG_DFIM(M)
      IF (M.EQ.MIJ(IJ) .AND. MIJ(IJ).NE.ML) RHOWGDFTH=0.5*RHOWGDFTH
      IF(M.GT.MIJ(IJ)) RHOWGDFTH = 0.0 
!
      SINPLUS = MAX (SL(IJ,1,M),0.0_rk)
      SUMT = SINPLUS
      SUMX = SINPLUS*SINTHETA(1)
      SUMY = SINPLUS*COSTHETA(1)

      DO K = 2,SIZE(F,2)
         SUMT = SUMT + MAX (SL(IJ,K,M),0.0_rk)
         SUMX = SUMX + MAX (SL(IJ,K,M),0.0_rk)*SINTHETA(K)
         SUMY = SUMY + MAX (SL(IJ,K,M),0.0_rk)*COSTHETA(K)
      ENDDO

      TEMPS1   =  TEMPS1 + SUMT*RHOWGDFTH
      CMRHOWGDFTH = CM*RHOWGDFTH
      TEMPS2 = TEMPS2 + SUMX*CMRHOWGDFTH
      TEMPS3 = TEMPS3 + SUMY*CMRHOWGDFTH
   ENDDO


   PHIAW(IJ)   =  TEMPS1
   XSTRESS(IJ) = TEMPS2/MAX(rho_AIR,1.0_rk)
   YSTRESS(IJ) = TEMPS3/MAX(rho_AIR,1.0_rk)
END DO   ! END DO for IJ


   CALL TAU_PHI_HF_OPENACC(MIJ, USTAR, Z0, XLEVTAIL, TAU1, PHI1)


DO IJ = 1,SIZE(F,1)

!     2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!     ----------------------------------------------------

   COSW     = MAX(COS(THETA(1)-UDIR(IJ)),0.0_rk)
   TEMPS2 = F(IJ,1,MIJ(IJ))*COSW**2
   TEMPS1 = TEMPS2*COSW
   DO K=2,SIZE(F,2)
      TEMPS1 = TEMPS1+F(IJ,K,MIJ(IJ))*MAX(COS(THETA(K)-UDIR(IJ)),0.0_rk)**3
      TEMPS2 = TEMPS2+F(IJ,K,MIJ(IJ))*MAX(COS(THETA(K)-UDIR(IJ)),0.0_rk)**2
   END DO


   TAUHF = CONST*FR5(MIJ(IJ))*GM1*TEMPS1*TAU1(IJ)
   PHIHF = rho_AIR*CONST*FR5(MIJ(IJ))*TEMPS2*PHI1(IJ)
!
   PHIAW(IJ)   = PHIAW(IJ)   + PHIHF
   XSTRESS(IJ) = XSTRESS(IJ) + TAUHF*SIN(UDIR(IJ))
   YSTRESS(IJ) = YSTRESS(IJ) + TAUHF*COS(UDIR(IJ))
   TAUW(IJ) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)


   TAUW(IJ) = MIN(TAUW(IJ),USTAR(IJ)**2-EPS1)
   TAUW(IJ) = MAX(TAUW(IJ),0.0_rk)

ENDDO

END SUBROUTINE STRESSO_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TABU_SWELLFT

!**** *TABU_SWELLFT* - FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS

!     FABRICE ARDHUIN  IFREMER  2013

!*    PURPOSE.
!     --------
!     TO ESTIMATE FRICTION COEFFICIENTS IN OSCILLATORY BOUNDARY LAYERS

!     METHOD.
!     -------
!       TABULATION ON KELVIN FUNCTIONS.

!     EXTERNALS.
!     -----------

!     KERKEI  (zeroth order Kelvin function Ker and Kei)

! ----------------------------------------------------------------------

      INTEGER, PARAMETER :: NITER=100
      REAL(real_kind),    PARAMETER :: ABMIN=0.3
      REAL(real_kind),    PARAMETER :: ABMAX=8.0, KAPPA=0.40
!     VARIABLE.   TYPE.     PURPOSE.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
! ----------------------------------------------------------------------
      INTEGER :: I,ITER
      REAL(real_kind) :: DELAB
      REAL(real_kind) :: KER, KEI
      REAL(real_kind) :: ABR,ABRLOG,L10,FACT,FSUBW,FSUBWMEMO,DZETA0,DZETA0MEMO

! ----------------------------------------------------------------------

      DZETA0 = 0.0
!
      DELAB   = (ABMAX-ABMIN)/REAL(IAB)
      L10=ALOG(10.0)
      DO I=1,IAB
         ABRLOG=ABMIN+REAL(I)*DELAB
         ABR=EXP(ABRLOG*L10)
         FACT=1/ABR/(21.2*KAPPA)
         FSUBW=0.05
         DO ITER=1,NITER
            FSUBWMEMO=FSUBW
            DZETA0MEMO=DZETA0
            DZETA0=FACT*FSUBW**(-0.5)
            CALL KERKEI(2.0*SQRT(DZETA0),KER,KEI)
            FSUBW=0.08/(KER**2+KEI**2)
            FSUBW=0.5*(FSUBWMEMO+FSUBW)
            DZETA0=0.5*(DZETA0MEMO+DZETA0)
         ENDDO   
         SWELLFT(I)=FSUBW
      ENDDO

END SUBROUTINE TABU_SWELLFT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INIT_SDISSP_ARD

! ----------------------------------------------------------------------

!**** *INIT_SDISSP_ARD* - INITIALISATION FOR SDISS_ARD


 INTEGER :: JD, K, M, I_INT, J_INT, M2, KK, KLD

 REAL(real_kind) :: TPIINV, TMP02
 REAL(real_kind) :: DELTH_TRUNC, DELTH_LOC

 REAL(real_kind) :: XLOGDFRTH
 REAL(real_kind) :: BRLAMBDA
 REAL(real_kind), DIMENSION(0:KL/2) :: COSDTH
 REAL(real_kind), DIMENSION(ML) :: SIG, C_, C_C, CGM1, DSIP, TRPZ_DSIP 

! ----------------------------------------------------------------------

      TPIINV = 1.0/Pi2

      KLD=KL/2

      XLOGDFRTH=LOG(CO)*delTheta

!     COMPUTE SATWEIGHTS

!     l(k,th)=1/(2*pi²)= the breaking crest density
      BRLAMBDA=BRKPBCOEF/(2.0*Pi2**2)

      TMP02 = SSDSC3*BRLAMBDA

      NSDSNTH  = MIN(NINT(ISDSDTH*RAD/(delTheta)),KLD-1)
      DELTH_TRUNC=(THETA(1)+ISDSDTH*RAD)-(THETA(1+NSDSNTH)-0.5*delTHeta)
      DELTH_TRUNC=MAX(0.0_rk,MIN(DELTH_TRUNC,delTheta))

      IF(ALLOCATED(INDICESSAT)) DEALLOCATE(INDICESSAT)
      ALLOCATE(INDICESSAT(KL,NSDSNTH*2+1))
      IF(ALLOCATED(SATWEIGHTS)) DEALLOCATE(SATWEIGHTS)
      ALLOCATE(SATWEIGHTS(KL,NSDSNTH*2+1))

      DO K=1,KL
        DO I_INT=K-NSDSNTH, K+NSDSNTH
          J_INT=I_INT
          IF (I_INT.LT.1)  J_INT=I_INT+KL
          IF (I_INT.GT.KL) J_INT=I_INT-KL
          INDICESSAT(K,I_INT-(K-NSDSNTH)+1)=J_INT

          IF(I_INT.EQ.K-NSDSNTH .OR. I_INT.EQ.K+NSDSNTH) THEN
            DELTH_LOC=DELTH_TRUNC
          ELSE
            DELTH_LOC=delTheta
          ENDIF
          SATWEIGHTS(K,I_INT-(K-NSDSNTH)+1)=DELTH_LOC*COS(THETA(K)-THETA(J_INT))**ISB
        END DO
      END DO

!!!!     COMPUTE CUMULW (only if needed)
!!!      IF (SSDSC3.NE.0.0) THEN
!!!        IF(ALLOCATED(CUMULW)) DEALLOCATE(CUMULW)
!!!        ALLOCATE(CUMULW(NDEPTH,0:KLD,ML,ML))
!!!
!!!!       NDIKCUMUL is the  integer difference in frequency bands
!!!!       between the "large breakers" and short "wiped-out waves"
!!!!!! wrong !!???        NDIKCUMUL = NINT(SSDSBRF1/(CO-1.))
!!!        NDIKCUMUL = NINT(-LOG(SSDSBRF1)/LOG(ratioFre))
!!!
!!!        DO KK=0,KLD
!!!          COSDTH(KK)=COS(KK*DELTH)
!!!        ENDDO
!!!
!!!        DO M=1,ML
!!!          SIG(M) = Pi2*Fre(M)
!!!        ENDDO
!!!
!!!        DO JD=1,NDEPTH
!!!
!!!          IF (SHALLOW_RUN) THEN
!!!            DO M=1,ML
!!!              C_(M)=SIG(M)/TFAK(JD,M)
!!!              CGM1(M)=1.0/TCGOND(JD,M)
!!!            ENDDO
!!!          ELSE
!!!            DO M=1,ML
!!!              C_(M)=G/SIG(M)  ! Valid in deep water only !!!!!!!!!!!!
!!!              CGM1(M)=2.0/C_(M) ! deep water !
!!!            ENDDO
!!!          ENDIF
!!!
!!!          DO M=1,ML
!!!            C_C(M)=C_(M)*C_(M)
!!!            DSIP(M)=TMP02*SIG(M)*XLOGDFRTH*CGM1(M) !  coef*dtheta*dk = coef*dtheta*dsigma/cg
!!!          ENDDO
!!!
!!!          DO M=NDIKCUMUL+1,ML
!!!
!!!            IF(M-NDIKCUMUL.GE.3) THEN
!!!              TRPZ_DSIP(1)=0.5*DSIP(1)
!!!              DO M2=2,M-NDIKCUMUL-1
!!!                TRPZ_DSIP(M2)=DSIP(M2)
!!!              ENDDO
!!!              TRPZ_DSIP(M-NDIKCUMUL)=0.5*DSIP(M-NDIKCUMUL)
!!!            ELSE
!!!              DO M2=1,M-NDIKCUMUL
!!!                TRPZ_DSIP(M2)=DSIP(M2)
!!!              ENDDO
!!!            ENDIF
!!!
!!!            DO M2=1,M-NDIKCUMUL
!!!              DO KK=0,KLD
!!!                CUMULW(JD,KK,M2,M)=SQRT(ABS(C_C(M)+C_C(M2)-2.0*C_(M)*C_(M2)*COSDTH(KK)))*TRPZ_DSIP(M2)
!!!              ENDDO 
!!!            ENDDO
!!!          ENDDO
!!!
!!!        ENDDO ! JD
!!!
!!!      ENDIF

 END SUBROUTINE INIT_SDISSP_ARD

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FRCUTINDEX (FM, FMWS, USTAR, MIJ)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FRCUTINDEX - RETURNS THE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!
!     METHOD.
!     -------
!                                                                              !
!     COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.

!!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
!!! re-activated !!!
!                                                                              !
!                                                                              !
!     EXTERNALS.
!     ---------
!                                                                              !
!     REFERENCE.
!     ----------
!                                                                              !
! ----------------------------------------------------------------------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.
!     --------------------
!                                                                              !
REAL(real_kind),    INTENT(IN)    :: FM(:)          !! MEAN FREQUENCY
REAL(real_kind),    INTENT(IN)    :: FMWS(:)        !! MEAN FREQUENCY OF WINDSEA
REAL(real_kind),    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY
INTEGER, INTENT(OUT)   :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !
!                                                                              !
INTEGER :: IJ

REAL(real_kind), PARAMETER :: EPSUS = 1.0E-6
REAL(real_kind), PARAMETER :: FRIC = 28.0
REAL(real_kind), PARAMETER :: TAILFACTOR = 2.5
REAL(real_kind), PARAMETER :: TAILFACTOR_PM = 3.0

REAL(real_kind) :: FPMH, FPPM, FM2, FPM, FPM4

! -----------------------------------------------------------------------------!
!                                                                              !
!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
!*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
!*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
!     ------------------------------------------------------------

FPMH = TAILFACTOR/FRE(1)
FPPM = TAILFACTOR_PM*G/(FRIC*Pi2*Fre(1))

DO IJ = 1,SIZE(USTAR)
    FM2 = MAX(FMWS(IJ),FM(IJ))*FPMH
    FPM = FPPM/MAX(USTAR(IJ),EPSUS)
    FPM4 = MAX(FM2,FPM)
    MIJ(IJ) = NINT(LOG10(FPM4)*INV_LOG_CO)+1
    MIJ(IJ) = MIN(MAX(1,MIJ(IJ)),ML)
ENDDO

END SUBROUTINE FRCUTINDEX

SUBROUTINE FRCUTINDEX_OPENACC (FM, FMWS, USTAR, MIJ)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FRCUTINDEX - RETURNS THE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!
!     METHOD.
!     -------
!                                                                              !
!     COMPUTES LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.

!!! be aware that if this is NOT used, for iphys=1, the cumulative dissipation has to be
!!! re-activated !!!
!                                                                              !
!                                                                              !
!     EXTERNALS.
!     ---------
!                                                                              !
!     REFERENCE.
!     ----------
!                                                                              !
! ----------------------------------------------------------------------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.
!     --------------------
!                                                                              !
REAL(real_kind),    INTENT(IN)    :: FM(:)          !! MEAN FREQUENCY
REAL(real_kind),    INTENT(IN)    :: FMWS(:)        !! MEAN FREQUENCY OF WINDSEA
REAL(real_kind),    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY
INTEGER, INTENT(OUT)   :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !
!                                                                              !
INTEGER :: IJ

REAL(real_kind), PARAMETER :: EPSUS = 1.0E-6
REAL(real_kind), PARAMETER :: FRIC = 28.0
REAL(real_kind), PARAMETER :: TAILFACTOR = 2.5
REAL(real_kind), PARAMETER :: TAILFACTOR_PM = 3.0

REAL(real_kind) :: FPMH, FPPM, FM2, FPM, FPM4

! -----------------------------------------------------------------------------!
!                                                                              !
!*    COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
!*    FREQUENCIES LE MAX(TAILFACTOR*MAX(FMNWS,FM),TAILFACTOR_PM*FPM),
!*    WHERE FPM IS THE PIERSON-MOSKOWITZ FREQUENCY BASED ON FRICTION
!*    VELOCITY. (FPM=G/(FRIC*ZPI*USTAR))
!     ------------------------------------------------------------


FPMH = TAILFACTOR/FRE(1)
FPPM = TAILFACTOR_PM*G/(FRIC*PI2*FRE(1))

DO IJ = 1,SIZE(USTAR)
    FM2 = MAX(FMWS(IJ),FM(IJ))*FPMH
    FPM = FPPM/MAX(USTAR(IJ),EPSUS)
    FPM4 = MAX(FM2,FPM)
    MIJ(IJ) = NINT(LOG10(FPM4)*INV_LOG_CO)+1
    MIJ(IJ) = MIN(MAX(1,MIJ(IJ)),ML)
ENDDO

END SUBROUTINE FRCUTINDEX_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE IMPHFTAIL_OPENACC (MIJ, FL3) 

! ----------------------------------------------------------------------

!**** *IMPHFTAIL* - IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM

!*    PURPOSE.
!     --------

!     IMPOSE A HIGH FREQUENCY TAIL TO THE SPECTRUM ABOVE FREQUENCY INDEX MIJ

!     METHOD.
!     -------

!     EXTERNALS.
!     ---------

!     REFERENCE.
!     ----------

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

INTEGER, INTENT(IN)    :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE
                                         !! PROGNOSTIC RANGE.
REAL(real_kind),    INTENT(INOUT) :: FL3 (:, :, :)  !! SPECTRUM.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER :: NFRE
INTEGER :: IJ, K, M, I

REAL :: AKM1, TFAC
REAL, DIMENSION(SIZE(FL3,1),SIZE(FL3,3)) :: TEMP2

! ----------------------------------------------------------------------

!*    DIAGNOSTIC TAIL.
!     ----------------

NFRE = SIZE(FL3,3)


IF (SHALLOW_RUN) THEN
  DO IJ = 1,SIZE(FL3,1)
    DO M=MIJ(IJ),NFRE
      I = cellIdxTab(IJ)
      AKM1 = 1./kwave(I,M)
      TEMP2(IJ,M) = AKM1**3/Cg(I,M)
    ENDDO
  ENDDO
ELSE
  DO IJ = 1,SIZE(FL3,1)
    DO M=MIJ(IJ),NFRE
      TEMP2(IJ,M) = FRM5(M)
    ENDDO
  ENDDO
ENDIF

DO IJ = 1,SIZE(FL3,1)
  DO M=MIJ(IJ)+1,NFRE
    TEMP2(IJ,M) = TEMP2(IJ,M)/TEMP2(IJ,MIJ(IJ))
  ENDDO
ENDDO

!*    MERGE TAIL INTO SPECTRA.
!     ------------------------

DO K=1,SIZE(FL3,2)
  DO IJ = 1,SIZE(FL3,1)
    DO M=MIJ(IJ)+1,NFRE
      TFAC = FL3(IJ,K,MIJ(IJ))
      FL3(IJ,K,M) = TEMP2(IJ,M)*TFAC
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE IMPHFTAIL_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE SDEPTHLIM(DEPTH, EMEAN, FL3)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------
!     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
!     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].
REAL(real_kind),    INTENT(INOUT) :: EMEAN (:)      !! SPECTRAL VARIANCE 
REAL(real_kind),    INTENT(INOUT) :: FL3(:,:,:)     !! SPECTRUM.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !

INTEGER :: K, M
REAL(real_kind), DIMENSION(SIZE(FL3,1)) :: EM

! ----------------------------------------------------------------------

EM(:)=MIN((0.25*GAMD*DEPTH(:))**2/EMEAN(:),1.0_rk)
DO M = 1,SIZE(FL3,3)
  DO K = 1,SIZE(FL3,2)
     FL3(:,K,M) = FL3(:,K,M)*EM(:)
  ENDDO
ENDDO
EMEAN(:) = EM(:)*EMEAN(:)

END SUBROUTINE SDEPTHLIM

SUBROUTINE SDEPTHLIM_OPENACC(DEPTH, EMEAN, FL3)

! ----------------------------------------------------------------------

!*    PURPOSE.
!     --------
!     LIMITS THE SPECTRAL VARIANCE SUCH THAT THE TOTAL VARIANCE
!     DOES NOT EXCEED THE MAXIMUM WAVE VARIANCE ALLOWED FOR A GIVEN DEPTH

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)    :: DEPTH (:)      !! WATER DEPTH [M].
REAL(real_kind),    INTENT(INOUT) :: EMEAN (:)      !! SPECTRAL VARIANCE 
REAL(real_kind),    INTENT(INOUT) :: FL3(:,:,:)     !! SPECTRUM.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !

INTEGER :: K, M, IJ
REAL(real_kind), DIMENSION(SIZE(FL3,1)) :: EM

! ----------------------------------------------------------------------
DO IJ = 1,SIZE(FL3,1)
  EM(IJ)=MIN((0.25*GAMD*DEPTH(IJ))**2/EMEAN(IJ),1.0_rk)
  EMEAN(IJ) = EM(IJ)*EMEAN(IJ)
END DO
DO M = 1,SIZE(FL3,3)
  DO K = 1,SIZE(FL3,2)
     DO IJ = 1,SIZE(FL3,1)
        FL3(IJ,K,M) = FL3(IJ,K,M)*EM(IJ)
     END DO
  ENDDO
ENDDO

END SUBROUTINE SDEPTHLIM_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INIT_X0TAUHF

! ----------------------------------------------------------------------

!**** *INIT_X0TAUHF* -

!*    PURPOSE.
!     ---------

!     INITIALISATION FOR TAU_PHI_HF


!**   INTERFACE.
!     ----------

!       *CALL* *INIT_X0TAUHF

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------   

INTEGER :: J

REAL(real_kind) :: CONST1, X0, FF, F, DF

! ---------------------------------------------------------------------------- !


!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

!     find lowest limit for integration X0 *(G/USTAR)
!     ALPHA*X0**2*EXP(XKAPPA/(X0+ZALP))=1
      X0=0.005
      DO J=1,30
        FF=EXP(XKAPPA/(X0+ZALP))
        F=ALPHA*X0**2*FF-1.0
        IF (F.EQ.0.0) EXIT
        DF=ALPHA*FF*(2.0*X0-XKAPPA*(X0/(X0+ZALP))**2)
        X0=X0-F/DF
      ENDDO
      X0TAUHF=X0

      CONST1 = (BETAMAX/XKAPPA**2)/3.0

      ! Simpson Integration weights (JTOT_TAUHF must be odd) !
      WTAUHF(1)=CONST1
      DO J=2,JTOT_TAUHF-1,2
        WTAUHF(J)=4.0*CONST1
        WTAUHF(J+1)=2.0*CONST1
      ENDDO
      WTAUHF(JTOT_TAUHF)=CONST1

END SUBROUTINE INIT_X0TAUHF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TAU_PHI_HF(MIJ, USTAR, Z0, XLEVTAIL, TAUHF, PHIHF)

! ----------------------------------------------------------------------

!**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
!                                   HIGH-FREQUENCY ENERGY FLUX.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
!     JEAN BIDLOT  ECMWF  JANUARY 2017

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX
!       FOR BOTH ECMWF PHYSICS AND METREO FRANCE PHYSICS.


!     METHOD.
!     -------

!       IT NEEDS A CALL TO INIT_X0TAUHF TO INITIALISE
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ---------------------------------------------------------------------------- !

!
!     INTERFACE VARIABLES.
!
!     --------------------
!
INTEGER, INTENT(IN)    :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE
                                         !! PROGNOSTIC RANGE.
REAL(real_kind),    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY
REAL(real_kind),    INTENT(IN)    :: Z0(:)          !! ROUGHNESS LENGTH.
REAL(real_kind),    INTENT(IN)    :: XLEVTAIL(:)    !! TAIL LEVEL
REAL(real_kind),    INTENT(OUT)   :: TAUHF(:)       !! HIGH-FREQUENCY STRESS 
REAL(real_kind),    INTENT(OUT)   :: PHIHF(:)       !! HIGH-FREQUENCY ENERGY FLUX INTO OCEAN 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: J, IJ

REAL(real_kind), PARAMETER :: ZSUP = 0.0  !  LOG(1.)
! REAL :: BETAMAX, ZALP, ALPHA
REAL(real_kind) :: OMEGA, OMEGAC, OMEGACC
REAL(real_kind) :: X0G, UST, UST0, TAUW, TAUW0
REAL(real_kind) :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
REAL(real_kind) :: DELZ, ZINF
REAL(real_kind) :: FNC2, SQRTZ0OG, SQRTGZ0, GM1, GZ0, XLOGGZ0

! ----------------------------------------------------------------------

      GM1= 1.0/G

!     See INIT_X0TAUHF
      X0G=X0TAUHF*G

!*    COMPUTE THE INTEGRALS
!     ---------------------

      DO IJ=1,SIZE(USTAR)
        OMEGAC    = Pi2*Fre(MIJ(IJ))
        UST0      = USTAR(IJ)
        TAUW0     = UST0**2
        GZ0       = G*Z0(IJ)
        XLOGGZ0   = LOG(GZ0)
        OMEGACC   = MAX(OMEGAC,X0G/UST0)

        SQRTZ0OG  = SQRT(Z0(IJ)*GM1)
        SQRTGZ0   = 1.0/SQRTZ0OG
        YC        = OMEGACC*SQRTZ0OG
        ZINF      = LOG(YC)
        DELZ      = MAX((ZSUP-ZINF)/REAL(JTOT_TAUHF-1),0.0_rk)

        TAUHF(IJ)= 0.0_rk
        PHIHF(IJ)= 0.0_rk

        TAUW     = TAUW0
        UST      = UST0
        ! Intergrals are integrated following a change of variable : Z=LOG(Y)
        DO J=1,JTOT_TAUHF
          Y         = EXP(ZINF+REAL(J-1)*DELZ)
          OMEGA     = Y*SQRTGZ0
          CM1       = OMEGA*GM1
          ZX        = UST*CM1 +ZALP
          ZARG      = XKAPPA/ZX
          ZLOG      = XLOGGZ0+2.0*LOG(CM1)+ZARG 
          ZLOG      = MIN(ZLOG,0.0_rk)
          ZBETA     = EXP(ZLOG)*ZLOG**4
          FNC2      = ZBETA*TAUW*WTAUHF(J)*DELZ
          TAUW      = MAX(TAUW-XLEVTAIL(IJ)*FNC2,0.0_rk)
          UST       = SQRT(TAUW)
          TAUHF(IJ) = TAUHF(IJ) + FNC2 
          PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
        ENDDO
        TAUHF(IJ) = TAUHF(IJ)
        PHIHF(IJ) = SQRTZ0OG*PHIHF(IJ)

      ENDDO

END SUBROUTINE TAU_PHI_HF

SUBROUTINE TAU_PHI_HF_OPENACC(MIJ, USTAR, Z0, XLEVTAIL, TAUHF, PHIHF)

! ----------------------------------------------------------------------

!**** *TAU_PHI_HF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
!                                   HIGH-FREQUENCY ENERGY FLUX.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
!     JEAN BIDLOT  ECMWF  JANUARY 2017

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS AND ENERGY FLUX
!       FOR BOTH ECMWF PHYSICS AND METREO FRANCE PHYSICS.


!     METHOD.
!     -------

!       IT NEEDS A CALL TO INIT_X0TAUHF TO INITIALISE
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ---------------------------------------------------------------------------- !

!
!     INTERFACE VARIABLES.
!
!     --------------------
!
INTEGER, INTENT(IN)    :: MIJ(:)         !! LAST FREQUENCY INDEX OF THE
                                         !! PROGNOSTIC RANGE.
REAL(real_kind),    INTENT(IN)    :: USTAR(:)       !! FRICTION VELOCITY
REAL(real_kind),    INTENT(IN)    :: Z0(:)          !! ROUGHNESS LENGTH.
REAL(real_kind),    INTENT(IN)    :: XLEVTAIL(:)    !! TAIL LEVEL
REAL(real_kind),    INTENT(OUT)   :: TAUHF(:)       !! HIGH-FREQUENCY STRESS 
REAL(real_kind),    INTENT(OUT)   :: PHIHF(:)       !! HIGH-FREQUENCY ENERGY FLUX INTO OCEAN 

! ---------------------------------------------------------------------------- !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: J, IJ

REAL(real_kind), PARAMETER :: ZSUP = 0.0_rk  !  LOG(1.)
!! REAL :: BETAMAX, ZALP, ALPHA
!REAL :: OMEGA, OMEGAC, OMEGACC
!REAL :: X0G, UST, UST0, TAUW, TAUW0
!REAL :: YC, Y, CM1, ZX, ZARG, ZLOG, ZBETA
!REAL :: DELZ, ZINF
!REAL :: FNC2, SQRTZ0OG, SQRTGZ0, GM1, GZ0, XLOGGZ0

REAL(real_kind) :: OMEGA, OMEGAC, OMEGACC
REAL(real_kind) :: X0G, UST, UST0, TAUW, TAUW0
REAL(real_kind) :: Y, CM1, ZX, ZLOG
REAL(real_kind) :: DELZ, ZINF
REAL(real_kind) :: FNC2, SQRTZ0OG, GM1, GZ0, XLOGGZ0

! ----------------------------------------------------------------------

      GM1= 1.0/G

!     See INIT_X0TAUHF
      X0G=X0TAUHF*G

!*    COMPUTE THE INTEGRALS
!     ---------------------

      DO IJ=1,SIZE(USTAR)
        OMEGAC    = PI2*FRE(MIJ(IJ))
        UST0      = USTAR(IJ)
        TAUW0     = UST0**2
        GZ0       = G*Z0(IJ)
        XLOGGZ0   = LOG(GZ0)
        OMEGACC   = MAX(OMEGAC,X0G/UST0)

        SQRTZ0OG  = SQRT(Z0(IJ)*GM1)
        !SQRTGZ0   = 1.0/SQRTZ0OG
        !YC        = OMEGACC*SQRTZ0OG
        ZINF      = LOG(OMEGACC*SQRTZ0OG)
        DELZ      = MAX((ZSUP-ZINF)/REAL(JTOT_TAUHF-1),0.0_rk)

        TAUHF(IJ)= 0.0_rk
        PHIHF(IJ)= 0.0_rk

        TAUW     = TAUW0
        UST      = UST0
        ! Intergrals are integrated following a change of variable : Z=LOG(Y)
        !add YUAN, TAUW is not cycle-independent,sequential
        DO J=1,JTOT_TAUHF
          Y         = EXP(ZINF+REAL(J-1)*DELZ)
          OMEGA     = Y/SQRTZ0OG
          CM1       = OMEGA*GM1
          ZX        = UST*CM1 +ZALP
          !ZARG      = XKAPPA/ZX
          ZLOG      = XLOGGZ0+2.0*LOG(CM1) + XKAPPA/ZX 
          ZLOG      = MIN(ZLOG,0.0_rk)
          !ZBETA     = EXP(ZLOG)*ZLOG**4
          FNC2      = EXP(ZLOG)*ZLOG**4 * TAUW*WTAUHF(J)*DELZ
          TAUW      = MAX(TAUW-XLEVTAIL(IJ)*FNC2,0.0_rk)
          UST       = SQRT(TAUW)
          TAUHF(IJ) = TAUHF(IJ) + FNC2 
          PHIHF(IJ) = PHIHF(IJ) + FNC2/Y
        ENDDO
        TAUHF(IJ) = TAUHF(IJ)
        PHIHF(IJ) = SQRTZ0OG*PHIHF(IJ)

      ENDDO

END SUBROUTINE TAU_PHI_HF_OPENACC
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE KERKEI (X, KER, KEI)

!**********************************************************************
! Computes the values of the zeroth order Kelvin function Ker and Kei
! These functions are used to determine the friction factor fw as a
! function of the bottom roughness length assuming a linear profile
! of eddy viscosity (See Grant and Madsen, 1979)
!**********************************************************************

      REAL(real_kind) :: ZR, ZI, CYR, CYI, CYR1, CYI1

      REAL(real_kind) :: X, KER, KEI

      ZR = X*0.50*SQRT(2.0)
      ZI = ZR
      CALL KZEONE (ZR, ZI, CYR, CYI, CYR1, CYI1)
      KER = CYR/EXP(ZR)
      KEI = CYI/EXP(ZR)
END SUBROUTINE KERKEI

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE KZEONE(X, Y, RE0, IM0, RE1, IM1)
!  June 1999 adaptation to CRESTb, all tests on range of (x,y) have been
!  bypassed, we implicitly expect X to be positive or |x,y| non zero
! 
! THE VARIABLES X AND Y ARE THE REAL AND IMAGINARY PARTS OF
! THE ARGUMENT OF THE FIRST TWO MODIFIED BESSEL FUNCTIONS
! OF THE SECOND KIND,K0 AND K1.  RE0,IM0,RE1 AND IM1 GIVE
! THE REAL AND IMAGINARY PARTS OF EXP(X)*K0 AND EXP(X)*K1,
! RESPECTIVELY.  ALTHOUGH THE REAL NOTATION USED IN THIS
! SUBROUTINE MAY SEEM INELEGANT WHEN COMPARED WITH THE
! COMPLEX NOTATION THAT FORTRAN ALLOWS, THIS VERSION RUNS
! ABOUT 30 PERCENT FASTER THAN ONE WRITTEN USING COMPLEX
! VARIABLES.
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IMPLICIT NONE
      REAL(real_kind) :: X, Y, X2, Y2, RE0, IM0, RE1, IM1,               &
     &        R1, R2, T1, T2, P1, P2, RTERM, ITERM, L

      REAL(real_kind), PARAMETER, DIMENSION(8) :: EXSQ =                 &
     &   (/ 0.5641003087264E0,  0.4120286874989E0,            &
     &      0.1584889157959E0,  0.3078003387255E-1,           &
     &      0.2778068842913E-2, 0.1000044412325E-3,           &
     &      0.1059115547711E-5, 0.1522475804254E-8 /)

      REAL(real_kind), PARAMETER, DIMENSION(8) :: TSQ =                  &
     &   (/ 0.0E0,              3.19303633920635E-1,          &
     &      1.29075862295915E0, 2.95837445869665E0,           &
     &      5.40903159724444E0, 8.80407957805676E0,           &
     &      1.34685357432515E1, 2.02499163658709E1 /)

      INTEGER :: N, M, K, LL
! THE ARRAYS TSQ AND EXSQ CONTAIN THE SQUARE OF THE
! ABSCISSAS AND THE WEIGHT FACTORS USED IN THE GAUSS-
! HERMITE QUADRATURE.
      R2 = X*X + Y*Y
      IF (R2.GE.1.96E2) GO TO 50
      IF (R2.GE.1.849E1) GO TO 30
! THIS SECTION CALCULATES THE FUNCTIONS USING THE SERIES
! EXPANSIONS
      X2 = X/2.0E0
      Y2 = Y/2.0E0
      P1 = X2*X2
      P2 = Y2*Y2
      T1 = -(LOG(P1+P2)/2.0E0+0.5772156649015329E0)
! THE CONSTANT IN THE PRECEDING STATEMENT IS EULER*S
! CONSTANT
      T2 = -ATAN2(Y,X)
      X2 = P1 - P2
      Y2 = X*Y2
      RTERM = 1.0E0
      ITERM = 0.0E0
      RE0 = T1
      IM0 = T2
      T1 = T1 + 0.5E0
      RE1 = T1
      IM1 = T2
      P2 = SQRT(R2)
      L = 2.106E0*P2 + 4.4E0
      IF (P2.LT.8.0E-1) L = 2.129E0*P2 + 4.0E0
      LL=NINT(L)
      DO N=1,LL
        P1 = N
        P2 = N*N
        R1 = RTERM
        RTERM = (R1*X2-ITERM*Y2)/P2
        ITERM = (R1*Y2+ITERM*X2)/P2
        T1 = T1 + 0.5E0/P1
        RE0 = RE0 + T1*RTERM - T2*ITERM
        IM0 = IM0 + T1*ITERM + T2*RTERM
        P1 = P1 + 1.0E0
        T1 = T1 + 0.5E0/P1
        RE1 = RE1 + (T1*RTERM-T2*ITERM)/P1
        IM1 = IM1 + (T1*ITERM+T2*RTERM)/P1
      END DO
      R1 = X/R2 - 0.5E0*(X*RE1-Y*IM1)
      R2 = -Y/R2 - 0.5E0*(X*IM1+Y*RE1)
      P1 = EXP(X)
      RE0 = P1*RE0
      IM0 = P1*IM0
      RE1 = P1*R1
      IM1 = P1*R2
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE INTEGRAL
! REPRESENTATION, EQN 3, EVALUATED WITH 15 POINT GAUSS-
! HERMITE QUADRATURE
   30 X2 = 2.0E0*X
      Y2 = 2.0E0*Y
      R1 = Y2*Y2
      P1 = SQRT(X2*X2+R1)
      P2 = SQRT(P1+X2)
      T1 = EXSQ(1)/(2.0E0*P1)
      RE0 = T1*P2
      IM0 = T1/P2
      RE1 = 0.0E0
      IM1 = 0.0E0
      DO N=2,8
        T2 = X2 + TSQ(N)
        P1 = SQRT(T2*T2+R1)
        P2 = SQRT(P1+T2)
        T1 = EXSQ(N)/P1
        RE0 = RE0 + T1*P2
        IM0 = IM0 + T1/P2
        T1 = EXSQ(N)*TSQ(N)
        RE1 = RE1 + T1*P2
        IM1 = IM1 + T1/P2
      END DO
      T2 = -Y2*IM0
      RE1 = RE1/R2
      R2 = Y2*IM1/R2
      RTERM = 1.41421356237309E0*COS(Y)
      ITERM = -1.41421356237309E0*SIN(Y)
! THE CONSTANT IN THE PREVIOUS STATEMENTS IS,OF COURSE,
! SQRT(2.0).
      IM0 = RE0*ITERM + T2*RTERM
      RE0 = RE0*RTERM - T2*ITERM
      T1 = RE1*RTERM - R2*ITERM
      T2 = RE1*ITERM + R2*RTERM
      RE1 = T1*X + T2*Y
      IM1 = -T1*Y + T2*X
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE
! ASYMPTOTIC EXPANSIONS
   50 RTERM = 1.0E0
      ITERM = 0.0E0
      RE0 = 1.0E0
      IM0 = 0.0E0
      RE1 = 1.0E0
      IM1 = 0.0E0
      P1 = 8.0E0*R2
      P2 = SQRT(R2)
      L = 3.91E0+8.12E1/P2
      LL=NINT(L)
      R1 = 1.0E0
      R2 = 1.0E0
      M = -8
      K = 3
      DO N=1,LL
        M = M + 8
        K = K - M
        R1 = FLOAT(K-4)*R1
        R2 = FLOAT(K)*R2
        T1 = FLOAT(N)*P1
        T2 = RTERM
        RTERM = (T2*X+ITERM*Y)/T1
        ITERM = (-T2*Y+ITERM*X)/T1
        RE0 = RE0 + R1*RTERM
        IM0 = IM0 + R1*ITERM
        RE1 = RE1 + R2*RTERM
        IM1 = IM1 + R2*ITERM
      END DO
      T1 = SQRT(P2+X)
      T2 = -Y/T1
      P1 = 8.86226925452758E-1/P2
! THIS CONSTANT IS SQRT(PI)/2.0, WITH PI=3.14159...
      RTERM = P1*COS(Y)
      ITERM = -P1*SIN(Y)
      R1 = RE0*RTERM - IM0*ITERM
      R2 = RE0*ITERM + IM0*RTERM
      RE0 = T1*R1 - T2*R2
      IM0 = T1*R2 + T2*R1
      R1 = RE1*RTERM - IM1*ITERM
      R2 = RE1*ITERM + IM1*RTERM
      RE1 = T1*R1 - T2*R2
      IM1 = T1*R2 + T2*R1

END SUBROUTINE KZEONE

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WSIGSTAR (USTAR, Z0, WSTAR, SIG_N)

! ----------------------------------------------------------------------

!**** *WSIGSTAR* - COMPUTATION OF THE RELATIVE STANDARD DEVIATION OF USTAR.

!*    PURPOSE.
!     ---------

!     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS
!     RELATIVE TO USTAR

!     METHOD.
!     -------

!     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
!     USTAR AND (Zi/L) THE MONIN-OBOKHOV LENGTH (Zi THE INVERSION HEIGHT).
!     (but with the background gustiness set to 0.)
!     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
!     WITH CD=A+B*U10 (see below).

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------

!     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.
!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)

! ---------------------------------------------------------------------------- !

!
!     INTERFACE VARIABLES.
!
!     --------------------
!

REAL(real_kind),    INTENT(IN)    :: USTAR(:) !! FRICTION VELOCITY
REAL(real_kind),    INTENT(IN)    :: Z0(:)    !! ROUGHNESS LENGTH.
REAL(real_kind),    INTENT(IN)    :: WSTAR(:) !! FREE CONVECTION VELOCITY SCALE.
REAL(real_kind),    INTENT(OUT)   :: SIG_N(:) !! RELATIVE STANDARD DEVIATION OF USTAR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL(real_kind), PARAMETER :: WSPMIN = 1.0 ! MINIMUM WIND SPEED
REAL(real_kind), PARAMETER :: BG_GUST = 0.0 ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
REAL(real_kind), PARAMETER :: ONETHIRD = 1.0/3.0
REAL(real_kind), PARAMETER :: LOG10 = LOG(10.0)
REAL(real_kind), PARAMETER :: C1 = 1.03E-3
REAL(real_kind), PARAMETER :: C2 = 0.04E-3
REAL(real_kind), PARAMETER :: P1 = 1.48
REAL(real_kind), PARAMETER :: P2 = -0.21

INTEGER :: IJ

REAL(real_kind) :: U10, C_D, DC_DDU, SIG_CONV
REAL(real_kind) :: XKAPPAD
REAL(real_kind) :: U10M1, C2U10P1, U10P2

! ---------------------------------------------------------------------------- !

      XKAPPAD=1.0/XKAPPA

      DO IJ=1,SIZE(USTAR)
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        U10 = USTAR(IJ)*XKAPPAD*(LOG10-LOG(Z0(IJ)))
        U10 = MAX(U10,WSPMIN)
        U10M1=1.0_rk/U10
        C2U10P1=C2*U10**P1
        U10P2=U10**P2
        C_D = (C1 + C2U10P1)*U10P2
        DC_DDU = (P2*C1+(P1+P2)*C2U10P1)*U10P2*U10M1
        SIG_CONV = 1.0_rk + 0.5_rk*U10/C_D*DC_DDU
        SIG_N(IJ) = MIN(0.1_rk, SIG_CONV * U10M1*(BG_GUST*USTAR(IJ)**3 + &
     &                  0.5_rk*XKAPPA*WSTAR(IJ)**3)**ONETHIRD )
      ENDDO

END SUBROUTINE WSIGSTAR

SUBROUTINE WSIGSTAR_OPENACC (USTAR, Z0, WSTAR, SIG_N)

! ----------------------------------------------------------------------

!**** *WSIGSTAR* - COMPUTATION OF THE RELATIVE STANDARD DEVIATION OF USTAR.

!*    PURPOSE.
!     ---------

!     COMPUTES THE STANDARD DEVIATION OF USTAR DUE TO SMALL SCALE GUSTINESS
!     RELATIVE TO USTAR

!     METHOD.
!     -------

!     USE PANOFSKY (1991) TO EXPRESS THE STANDARD DEVIATION OF U10 IN TERMS
!     USTAR AND (Zi/L) THE MONIN-OBOKHOV LENGTH (Zi THE INVERSION HEIGHT).
!     (but with the background gustiness set to 0.)
!     and USTAR=SQRT(Cd)*U10 to DERIVE THE STANDARD DEVIATION OF USTAR.
!     WITH CD=A+B*U10 (see below).

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     REFERENCE.
!     ----------

!     SEE SECTION 3.2.1 OF THE WAM DOCUMENTATION.
!     USE HERSBACH 2011 FOR CD(U10) (SEE ALSO EDSON et al. 2013)

! ---------------------------------------------------------------------------- !

!
!     INTERFACE VARIABLES.
!
!     --------------------
!

REAL(real_kind),    INTENT(IN)    :: USTAR(:) !! FRICTION VELOCITY
REAL(real_kind),    INTENT(IN)    :: Z0(:)    !! ROUGHNESS LENGTH.
REAL(real_kind),    INTENT(IN)    :: WSTAR(:) !! FREE CONVECTION VELOCITY SCALE.
REAL(real_kind),    INTENT(OUT)   :: SIG_N(:) !! RELATIVE STANDARD DEVIATION OF USTAR

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL(real_kind), PARAMETER :: WSPMIN = 1.0 ! MINIMUM WIND SPEED
REAL(real_kind), PARAMETER :: BG_GUST = 0.0 ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
REAL(real_kind), PARAMETER :: ONETHIRD = 1.0/3.0
REAL(real_kind), PARAMETER :: LOG10 = LOG(10.0)
REAL(real_kind), PARAMETER :: C1 = 1.03E-3
REAL(real_kind), PARAMETER :: C2 = 0.04E-3
REAL(real_kind), PARAMETER :: P1 = 1.48
REAL(real_kind), PARAMETER :: P2 = -0.21

INTEGER :: IJ

REAL(real_kind) :: U10, C_D, DC_DDU, SIG_CONV
REAL(real_kind) :: XKAPPAD
REAL(real_kind) :: U10M1, C2U10P1, U10P2

! ---------------------------------------------------------------------------- !

      XKAPPAD=1.0/XKAPPA

      DO IJ=1,SIZE(USTAR)
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        U10 = USTAR(IJ)*XKAPPAD*(LOG10-LOG(Z0(IJ)))
        U10 = MAX(U10,WSPMIN)
        U10M1=1.0/U10
        C2U10P1=C2*U10**P1
        U10P2=U10**P2
        C_D = (C1 + C2U10P1)*U10P2
        DC_DDU = (P2*C1+(P1+P2)*C2U10P1)*U10P2*U10M1
        SIG_CONV = 1.0 + 0.5*U10/C_D*DC_DDU
        SIG_N(IJ) = MIN(0.1_rk, SIG_CONV * U10M1*(BG_GUST*USTAR(IJ)**3 + &
     &                  0.5*XKAPPA*WSTAR(IJ)**3)**ONETHIRD )
      ENDDO

END SUBROUTINE WSIGSTAR_OPENACC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TOTAL_ENERGY_B (F3, EMEAN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TOTAL_ENERGY_B - COMPUTES TOTAL ENERGY (VECTOR VERSION).                   !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OVER DIRECTION AND FREQUENCY. A TAIL CORRECTION IS ADDED.  !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)            :: F3(:,:,:)    !! BLOCK OF SPECTRA.
REAL(real_kind),    INTENT(OUT)           :: EMEAN(:)     !! TOTAL ENERGY.
LOGICAL, INTENT(IN), OPTIONAL  :: MASK(:,:,:)  !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IJ, M, ML
REAL(real_kind)     :: TEMP(SIZE(F3,1),SIZE(F3,3))



ML = SIZE(F3,3)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTION (WITHOUT DELTH).                             !
!        -----------------------------------------                             !

IF (PRESENT(MASK)) THEN
   TEMP(:,:) = SUM(F3(:,:,:), DIM=2, MASK=MASK(:,:,:))
ELSE
   TEMP(:,:) = SUM(F3(:,:,:), DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !

EMEAN(:) = MO_TAIL*TEMP(:,ML)    !! TAIL ENERGY

DO M = 1,ML
   DO IJ = 1,SIZE(F3,1)
      EMEAN(IJ) = EMEAN(IJ) + TEMP(IJ,M)*DFIM(M)
   END DO
END DO

EMEAN(:) = MAX(EMEAN(:), EMIN)

END SUBROUTINE TOTAL_ENERGY_B

SUBROUTINE TOTAL_ENERGY_B_OPENACC (F3, EMEAN, MASK)
! ---------------------------------------------------------------------------- !
!                                                                              !
!   TOTAL_ENERGY_B - COMPUTES TOTAL ENERGY (VECTOR VERSION).                   !
!                                                                              !
!     S.D. HASSELMANN.                                                         !
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER                               !
!     H. GUENTHER     GKSS   DECEMBER 2001  FT90                               !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OVER DIRECTION AND FREQUENCY. A TAIL CORRECTION IS ADDED.  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)            :: F3(:,:,:)    !! BLOCK OF SPECTRA.
REAL(real_kind),    INTENT(OUT)           :: EMEAN(:)     !! TOTAL ENERGY.
LOGICAL, INTENT(IN), OPTIONAL  :: MASK(:,:,:)  !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IJ, M, K
REAL(real_kind)     :: TEMP(SIZE(F3,1),SIZE(F3,3))
REAL(real_kind)     :: TEMPS

IF (PRESENT(MASK)) THEN
!NOTE: abs(MASK) assumes .True. corresponding to 1 or -1. For some compilers
!.True. may correspond to -1. 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTION (WITHOUT DELTH).                             !
!        -----------------------------------------                             !
   DO M = 1,SIZE(F3,3)
      DO IJ = 1,SIZE(F3,1)
         TEMPS = 0.
         DO K = 1, SIZE(F3,2) 
            TEMPS = TEMPS + F3(IJ,K,M)*abs(MASK(IJ,K,M))
         END DO
         TEMP(IJ,M) = TEMPS
      END DO
   END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !
   DO IJ = 1,SIZE(F3,1)
      TEMPS = MO_TAIL*TEMP(IJ,ML)    !! TAIL ENERGY
      DO M = 1,ML
         TEMPS = TEMPS + TEMP(IJ,M)*DFIM(M)
      END DO
      EMEAN(IJ) = MAX(TEMPS,EMIN)
   END DO

ELSE

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTION (WITHOUT DELTH).                             !
!        -----------------------------------------                             !
   DO M = 1,SIZE(F3,3)
      DO IJ = 1,SIZE(F3,1)
         TEMPS = 0.
         DO K = 1,SIZE(F3,2) 
            TEMPS = TEMPS + F3(IJ,K,M)
         END DO
         TEMP(IJ,M) = TEMPS
      END DO
   END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !
   DO IJ = 1,SIZE(F3,1)
      TEMPS = MO_TAIL*TEMP(IJ,ML)    !! TAIL ENERGY
      DO M = 1,SIZE(F3,3)
         TEMPS = TEMPS + TEMP(IJ,M)*DFIM(M)
      END DO
      EMEAN(IJ) = MAX(TEMPS,EMIN)
   END DO
END IF

END SUBROUTINE TOTAL_ENERGY_B_OPENACC

SUBROUTINE FEMEAN (F, EMEAN, FM, MASK)

! ---------------------------------------------------------------------------- !
!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.                             !
REAL(real_kind),    INTENT(IN)            :: F(:,:,:)      !! BLOCK OF SPECTRA.
REAL(real_kind),    INTENT(IN)            :: EMEAN(:)      !! TOTAL ENERGY.
REAL(real_kind),    INTENT(OUT)           :: FM   (:)      !! MEAN FREQUENCY.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK (:,:,:)  !! INTERATION MASK.

REAL(real_kind)     :: TEMP2(SIZE(F,1),SIZE(F,3))
INTEGER :: IJ, M


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP2 = SUM(F,DIM=2, MASK=MASK)
ELSE
   TEMP2 = SUM(F,DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !

DO IJ = 1,SIZE(F,1)
   FM(IJ) = MM1_TAIL*TEMP2(IJ,ML)  !! TAIL ENERGY
END DO

DO M = 1,size(F,3)
   DO IJ = 1,SIZE(F,1)
      FM(IJ) = FM(IJ) + TEMP2(IJ,M)*DFIMOFR(M)
   END DO
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. NORMALIZE.                                                            !
!        ----------                                                            !

DO IJ = 1, SIZE(F,1)
   FM(IJ) = EMEAN(IJ)/MAX(FM(IJ),EMIN)
END DO

END SUBROUTINE FEMEAN

SUBROUTINE FEMEAN_OPENACC (F, EMEAN, FM, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FEMEAN - COMPUTATION OF MEAN FREQUENCY.                                    !
!                                                                              !
!     S.D. HASSELMANN                                                          !
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER                              !
!     H. GUNTHER     GKSS         DECEMBER 2001    FT90                        !
!                                                                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN)            :: F(:,:,:)      !! BLOCK OF SPECTRA.
REAL(real_kind),    INTENT(IN)            :: EMEAN(:)      !! TOTAL ENERGY.
REAL(real_kind),    INTENT(OUT)           :: FM   (:)      !! MEAN FREQUENCY.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK (:,:,:)  !! INTERATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL(real_kind)     :: TEMP2(SIZE(F,1),SIZE(F,3))
REAL(real_kind)     :: TEMPS
INTEGER :: IJ, M, K


IF (PRESENT(MASK)) THEN
!NOTE: abs(MASK) assumes .True. corresponding to 1 or -1. For some compilers
!.True. may correspond to -1. 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !
   DO M = 1,SIZE(F,3)
      DO IJ = 1,SIZE(F,1)
         TEMPS = 0.0_rk
         DO K = 1,SIZE(F,2) 
            TEMPS = TEMPS + F(IJ,K,M)*abs(MASK(IJ,K,M))
            !TEMPS = TEMPS + FTMP(IJ,K,M)
         END DO
         TEMP2(IJ,M) = TEMPS
      END DO
   END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !
   DO IJ = 1,SIZE(F,1)
      TEMPS = MM1_TAIL*TEMP2(IJ,ML)  !! TAIL ENERGY
      DO M = 1,SIZE(F,3)
         TEMPS = TEMPS + TEMP2(IJ,M)*DFIMOFR(M)
      END DO
! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. NORMALIZE.                                                            !
!        ----------                                                            !
      FM(IJ) = EMEAN(IJ)/MAX(TEMPS,EMIN)
END DO

ELSE


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !
   DO M = 1,SIZE(F,3)
      DO IJ = 1,SIZE(F,1)
         TEMPS = 0.0_rk
         DO K = 1,SIZE(F,2)  
            TEMPS = TEMPS + F(IJ,K,M)
         END DO
         TEMP2(IJ,M) = TEMPS
      END DO
   END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !
   DO IJ = 1,SIZE(F,1)
      TEMPS = MM1_TAIL*TEMP2(IJ,ML)  !! TAIL ENERGY
      DO M = 1,SIZE(F,3)
         TEMPS = TEMPS + TEMP2(IJ,M)*DFIMOFR(M)
      END DO
! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. NORMALIZE.                                                            !
!        ----------                                                            !
      FM(IJ) = EMEAN(IJ)/MAX(TEMPS,EMIN)
END DO

END IF

END SUBROUTINE FEMEAN_OPENACC

SUBROUTINE TM1_TM2_PERIODS_B (F, EMEAN, TM1, TM2, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TM1_TM2_PERIODS_B - COMPUTES TM1 AND/OR TM2 PERIODS (VECTOR VESION).       !
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TM1 AND TM2 PERIODS.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !

REAL(real_kind),    INTENT(IN )           :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL(real_kind),    INTENT(IN )           :: EMEAN(:)    !! TOTAL ENERGY [M*M].
REAL(real_kind),    INTENT(OUT), OPTIONAL :: TM1(:)      !! TM1 PERIOD [S].
REAL(real_kind),    INTENT(OUT), OPTIONAL :: TM2(:)      !! TM2 PERIOD [S].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, M
REAL(real_kind)    :: TEMP(SIZE(F,1),SIZE(F,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !
IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=2, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM1 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM1)) THEN
   TM1(:) = MP1_TAIL * TEMP(:,ML)     !! TAIL
   DO M = 1,ML
      DO IJ = 1,SIZE(TEMP,1)
         TM1(IJ) = TM1(IJ) + TEMP(IJ,M)*DFIM_FR(M)
      END DO
   END DO

   WHERE (EMEAN.GT.EMIN)                          !! NORMALIZE WITH ENERGY.
      TM1 = EMEAN/TM1
   ELSEWHERE
      TM1 = 1.
   END WHERE
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TM2 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM2)) THEN
   TM2(:) = MP2_TAIL * TEMP(:,ML)     !! TAIL
   DO M = 1,ML
      DO IJ = 1,SIZE(TEMP,1)
         TM2(IJ) = TM2(IJ) + TEMP(IJ,M)*DFIM_FR2(M)
      END DO
   END DO

   WHERE (EMEAN.GT.EMIN)                          !! NORMALIZE WITH ENERGY.
      TM2 = SQRT(EMEAN/TM2)
   ELSEWHERE
      TM2 = 1.
   END WHERE
END IF

END SUBROUTINE TM1_TM2_PERIODS_B

SUBROUTINE TM1_TM2_PERIODS_B_OPENACC (F, EMEAN, TM1, TM2, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TM1_TM2_PERIODS_B - COMPUTES TM1 AND/OR TM2 PERIODS (VECTOR VESION).       !
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TM1 AND TM2 PERIODS.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN )           :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL(real_kind),    INTENT(IN )           :: EMEAN(:)    !! TOTAL ENERGY [M*M].
REAL(real_kind),    INTENT(OUT), OPTIONAL :: TM1(:)      !! TM1 PERIOD [S].
REAL(real_kind),    INTENT(OUT), OPTIONAL :: TM2(:)      !! TM2 PERIOD [S].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, M, K
REAL(real_kind)    :: TEMP(SIZE(F,1),SIZE(F,3))
REAL(real_kind)    :: TEMPS


IF (PRESENT(MASK)) THEN
!NOTE: abs(MASK) assumes .True. corresponding to 1 or -1. For some compilers
!.True. may correspond to -1. 

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !
   DO M = 1,SIZE(F,3)
      DO IJ = 1,SIZE(F,1)
         TEMPS = 0.0_rk
         DO K = 1,SIZE(F,2) 
            TEMPS = TEMPS + F(IJ,K,M)*abs(MASK(IJ,K,M))
         END DO
         TEMP(IJ,M) = TEMPS
      END DO
   END DO

   IF (PRESENT(TM1)) THEN
! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM1 PERIOD.                                                           !
!        -----------                                                           !
   DO IJ = 1,SIZE(F,1)
      TEMPS = MP1_TAIL * TEMP(IJ,ML)     !! TAIL
      DO M = 1,SIZE(F,3)
         TEMPS = TEMPS + TEMP(IJ,M)*DFIM_FR(M)
      END DO
      IF (EMEAN(IJ).GT.EMIN) THEN
         TM1(IJ) = EMEAN(IJ)/TEMPS            !! NORMALIZE WITH ENERGY.
      ELSE
        TM1(IJ) = 1.0_rk
      END IF
   END DO
   END IF

   IF (PRESENT(TM2)) THEN
! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM2 PERIOD.                                                           !
!        -----------                                                           !
   DO IJ = 1,SIZE(F,1)
      TEMPS = MP2_TAIL * TEMP(IJ,ML)     !! TAIL
      DO M = 1,SIZE(F,3)
         TEMPS = TEMPS + TEMP(IJ,M)*DFIM_FR2(M)
      END DO
      IF (EMEAN(IJ).GT.EMIN) THEN
         TM2(IJ) = SQRT(EMEAN(IJ)/TEMPS)            !! NORMALIZE WITH ENERGY.
      ELSE
        TM2(IJ) = 1.0_rk
      END IF
   END DO
   END IF

ELSE        !! if MASK not present


! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !
   DO M = 1,SIZE(F,3)
      DO IJ = 1,SIZE(F,1)
         TEMPS = 0.0_rk
         DO K = 1,SIZE(F,2)
            TEMPS = TEMPS + F(IJ,K,M)
         END DO
         TEMP(IJ,M) = TEMPS
      END DO
   END DO

   IF (PRESENT(TM1)) THEN
! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM1 PERIOD.                                                           !
!        -----------                                                           !
   DO IJ = 1,SIZE(F,1)
      TEMPS = MP1_TAIL * TEMP(IJ,ML)     !! TAIL
      DO M = 1,ML
         TEMPS = TEMPS + TEMP(IJ,M)*DFIM_FR(M)
      END DO
      IF (EMEAN(IJ).GT.EMIN) THEN
         TM1(IJ) = EMEAN(IJ)/TEMPS            !! NORMALIZE WITH ENERGY.
      ELSE
        TM1(IJ) = 1.0_rk
      END IF
   END DO
   END IF

   IF (PRESENT(TM2)) THEN
! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM2 PERIOD.                                                           !
!        -----------                                                           !
   DO IJ = 1,SIZE(F,1)
      TEMPS = MP2_TAIL * TEMP(IJ,ML)     !! TAIL
      DO M = 1,SIZE(F,3)
         TEMPS = TEMPS + TEMP(IJ,M)*DFIM_FR2(M)
      END DO
      IF (EMEAN(IJ).GT.EMIN) THEN
         TM2(IJ) = SQRT(EMEAN(IJ)/TEMPS)            !! NORMALIZE WITH ENERGY.
      ELSE
        TM2(IJ) = 1.0_rk
      END IF
   END DO
   END IF

END IF

END SUBROUTINE TM1_TM2_PERIODS_B_OPENACC

SUBROUTINE WM1_WM2_WAVENUMBER_B_OPENACC(F, EMEAN, WM1, WM2, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WM1_WM2_WAVENUMBER_B - COMPUTES WM1 AND/OR WM2 WAVENUMBERS (VECTOR VESION).!
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE  WM1 AND/OR WM2 WAVENUMBERS                                    !
!          WM1 IS SQRT(1/K)*F INTGRATION                                       !
!          WM2 IS SQRT(K)*F INTGRATION                                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL(real_kind),    INTENT(IN )           :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL(real_kind),    INTENT(IN )           :: EMEAN(:)    !! TOTAL ENERGY [M*M].
REAL(real_kind),    INTENT(OUT), OPTIONAL :: WM1(:)      !! WM1 WAVENUMBER [M].
REAL(real_kind),    INTENT(OUT), OPTIONAL :: WM2(:)      !! WM2 WAVENUMBER [M].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER :: IJ, M, K, I
REAL(real_kind)    :: DEL2
REAL(real_kind)    :: TEMPS
REAL(real_kind)    :: TEMP(SIZE(F,1),SIZE(F,3)), TEMP2(SIZE(F,1),SIZE(F,3))

IF (PRESENT(MASK)) THEN
!NOTE: abs(MASK) assumes .True. corresponding to 1 or -1. For some compilers
!.True. may correspond to -1. 
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !
         DO M = 1,ML
            DO IJ = 1,SIZE(F,1)
               TEMPS = 0.
               DO K = 1,KL
                  TEMPS = TEMPS + F(IJ,K,M)*abs(MASK(IJ,K,M))
               END DO
               TEMP(IJ,M) = TEMPS
               I = cellIdxTab(IJ)
               TEMP2(IJ,M) =  SQRT(kwave(I,M))
            END DO
         END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM1.                                         !
!        -----------------------------                                         !
      IF (PRESENT(WM1)) THEN
         DO IJ = 1,SIZE(F,1)
            DEL2 = SQRT(G)/Pi2
            TEMPS = MM1_TAIL*DEL2*TEMP(IJ,ML)   !! TAIL.
            DO M = 1,ML
               TEMPS = TEMPS + TEMP(IJ,M) / TEMP2(IJ,M) * DFIM(M)
            END DO
            IF(EMEAN(IJ).GT.EMIN) THEN
               WM1(IJ) = (EMEAN(IJ)/TEMPS)**2        !! NORMALIZE WITH ENERGY.
            ELSE
               WM1(IJ) = 1.
            END IF
         END DO
      END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM2.                                         !
!        -----------------------------                                         !
      IF (PRESENT(WM2)) THEN
         DO IJ = 1,SIZE(TEMP,1)
            DEL2 = Pi2/SQRT(G)
            TEMPS = MP1_TAIL*DEL2*TEMP(IJ,ML)   !! ADD TAIL.
            DO M = 1,ML
               TEMPS = TEMPS + TEMP(IJ,M) * TEMP2(IJ,M) * DFIM(M)
            END DO
            IF(EMEAN(IJ).GT.EMIN) THEN
               WM2(IJ) = (TEMPS/EMEAN(IJ))**2        !! NORMALIZE WITH ENERGY.
            ELSE
               WM2(IJ) = 1.
            END IF
         END DO
      END IF


ELSE   ! if MASK not present       
! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !
         DO M = 1,ML
            DO IJ = 1,SIZE(F,1)
               TEMPS = 0.
               DO K = 1,KL
                  TEMPS = TEMPS + F(IJ,K,M)
               END DO
               TEMP(IJ,M) = TEMPS
               I = cellIdxTab(IJ)
               TEMP2(IJ,M) =  SQRT(kwave(I,M))
            END DO
         END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM1.                                         !
!        -----------------------------                                         !
      IF (PRESENT(WM1)) THEN
         DO IJ = 1,SIZE(F,1)
            DEL2 = SQRT(G)/Pi2
            TEMPS = MM1_TAIL*DEL2*TEMP(IJ,ML)   !! TAIL.
            DO M = 1,ML
               TEMPS = TEMPS + TEMP(IJ,M) / TEMP2(IJ,M) * DFIM(M)
            END DO
            IF(EMEAN(IJ).GT.EMIN) THEN
               WM1(IJ) = (EMEAN(IJ)/TEMPS)**2        !! NORMALIZE WITH ENERGY.
            ELSE
               WM1(IJ) = 1.
            END IF
         END DO
      END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM2.                                         !
!        -----------------------------                                         !
      IF (PRESENT(WM2)) THEN
         DO IJ = 1,SIZE(TEMP,1)
            DEL2 = Pi2/SQRT(G)
            TEMPS = MP1_TAIL*DEL2*TEMP(IJ,ML)   !! ADD TAIL.
            DO M = 1,ML
               TEMPS = TEMPS + TEMP(IJ,M) * TEMP2(IJ,M) * DFIM(M)
            END DO
            IF(EMEAN(IJ).GT.EMIN) THEN
               WM2(IJ) = (TEMPS/EMEAN(IJ))**2        !! NORMALIZE WITH ENERGY.
            ELSE
               WM2(IJ) = 1.
            END IF
         END DO
      END IF

END IF   ! ENDIF FOR PRESENT(MASK)

END SUBROUTINE WM1_WM2_WAVENUMBER_B_OPENACC


REAL(real_kind) FUNCTION TRANSF (XK, D)
    
    !    TRANSF   DETERMINE NARROW BAND LIMIT BENJAMIN-FEIR INDEX FOR
    
    !
    !           BF**2 = (2 S^2)/SIG_OM^2) . TRANSF2(XK,D)
    
    REAL(real_kind),    INTENT(IN)    :: XK   !! WAVE NUMBER
    REAL(real_kind),    INTENT(IN)    :: D    !! DEPTH
    REAL(real_kind) ::  EPS, X, T_0, OM, C_0, V_G, DV_G, XNL_1, XNL_2, XNL
    
    EPS = 0.0001
    !     1. DETERMINE TRANSFER FUNCTION.
    IF (D.LT.999. .AND. D.GT.0.) THEN
         X   = XK*D
         IF (X .GT. DKMAX) THEN
                TRANSF = 1.
         ELSE
            T_0 = TANH(X)
            OM  = SQRT(G*XK*T_0)
            C_0 = OM/XK
            IF (X .LT. EPS) THEN
                   V_G = 0.5*C_0
                   V_G = C_0
            ELSE
               V_G = 0.5*C_0*(1.+2.*X/SINH(2.*X))
            ENDIF
            DV_G = (T_0-X*(1.-T_0**2))**2+4.*X**2*T_0**2*(1.-T_0**2)
          
            XNL_1 = (9.*T_0**4-10.*T_0**2+9.)/(8.*T_0**3)
            XNL_2 = ((2.*V_G-0.5*C_0)**2/(G*D-V_G**2)+1.)/X
          
            XNL = XNL_1-XNL_2
            TRANSF = XNL**2/(DV_G*T_0**8)
         ENDIF
    ELSE
       TRANSF = 1.
    ENDIF
    
END FUNCTION TRANSF


INTEGER FUNCTION JAFU (CL, J, IAN)

    ! ------------------------------------------------------------------------- !
    !                                                                           !
    !   JAFU - FUNCTION TO COMPUTE THE INDEX ARRAY FOR THE ANGLES OF THE        !
    !          INTERACTING WAVENUMBERS.                                         !
    !                                                                           !
    !     S. HASSELMANN        MPIFM        01/12/1985.                         !
    !                                                                           !
    !     PURPOSE.                                                              !
    !     --------                                                              !
    !                                                                           !
    !       INDICES DEFINING BINS IN FREQUENCY AND DIRECTION PLANE INTO         !
    !       WHICH NONLINEAR ENERGY TRANSFER INCREMENTS ARE STORED. NEEDED       !
    !       FOR COMPUTATION OF THE NONLINEAR ENERGY TRANSFER.                   !
    !                                                                           !
    !     METHOD.                                                               !
    !     -------                                                               !
    !                                                                           !
    !       SEE REFERENCE.                                                      !
    !                                                                           !
    !     REFERENCE.                                                            !
    !     ----------                                                            !
    !                                                                           !
    !        S. HASSELMANN AND K. HASSELMANN,JPO, 1985 B.                       !
    !                                                                           !
    ! ------------------------------------------------------------------------- !
    
    REAL(real_kind),    INTENT(IN) :: CL    !! WEIGHTS.
    INTEGER, INTENT(IN) :: J     !! INDEX IN ANGULAR ARRAY.
    INTEGER, INTENT(IN) :: IAN   !! NUMBER OF ANGLES IN ARRAY.
    
    JAFU = J + INT(CL)
    IF (JAFU.LE.0)   JAFU = JAFU+IAN
    IF (JAFU.GT.IAN) JAFU = JAFU-IAN
    
END FUNCTION JAFU



END MODULE WAM_SOURCE_MODULE
