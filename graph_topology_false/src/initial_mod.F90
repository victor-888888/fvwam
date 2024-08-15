module initial_mod

  use params_mod
  use const_mod
  use mesh_mod
  use state_mod
  use mask_mod
  use forcing_mod
  use fre_dir_mod
  use restart_mod


  implicit none

  private

  public set_initial_condition
  public FETCH_LAW
  public JONSWAP
  
contains

subroutine set_initial_condition()
    implicit none
    integer :: iCell, K, M
    if (restart_input_option) then
        call restart_input(wave(1) )  
    else
        call cold_start_wam(forcing%u10spd, forcing%u10Dir, wave(1)%N) 
    end if

    DO M = 1,nFre
       DO K = 1,nDir
          DO iCell = 1, nCells
             wave(2 )%N(iCell,K,M) = wave(1)%N(iCell,K,M)
             wave(-1)%N(iCell,K,M) = wave(1)%N(iCell,K,M)
          END DO
       END DO
    END DO
end subroutine set_initial_condition

subroutine cold_start_wam(u10Spd, u10Dir, N0)
  implicit none
  real(real_kind), intent(in) :: u10Spd(:) ! m/s 
  real(real_kind), intent(in) :: u10Dir(:) ! radian
  real(real_kind), intent(out):: N0(0:,:,:) 

  REAL(real_kind), PARAMETER :: ZDP=2./PI
  REAL(real_kind), PARAMETER :: GAMMA = 3.000000
  REAL(real_kind), PARAMETER :: SA=7.000000E-02
  REAL(real_kind), PARAMETER :: SB=9.000000E-02
  
  REAL(real_kind) :: FM=0.2 !! PEAK FREQUENCY [HZ]
  
  INTEGER :: iCell, k, m
  real(real_kind) :: ST(nCells, nDir), ET(nCells, nFre)
  
  real(real_kind) :: FP(nCells), ALPHJ(nCells)
  
  !FM = sum(fre)/nFre
  !FM = maxval(fre)
  
  CALL FETCH_LAW(FETCH, FM, u10Spd, ALPHJ, FP)
  ! ---------------------------------------------------------------------------- !
  !                                                                              !
  !     1. COMPUTE JONSWAP SPECTRUM.                                             !
  !        -------------------------                                             !
  
  CALL JONSWAP(fre, ALPHJ, GAMMA, SA, SB, FP, ET)

  
  ! ---------------------------------------------------------------------------- !
  !                                                                              !
  !     2. COMPUTATION OF SPREADING FUNCTION.                                    !
  !        ----------------------------------                                    !
  
  DO iCell = 1, nCells
     DO K = 1, nDir
        ST(iCell,K) = MAX(0. ,COS(theta(K)-u10Dir(iCell)))
        ST(iCell,K) = ZDP*ST(iCell,K)**2
        IF (ST(iCell,K) .LT.0.1E-08) ST(iCell,K) = 0.0
     END DO
  END DO
  
  ! ---------------------------------------------------------------------------- !
  !                                                                              !
  !     3. COMPUTATION OF 2-D SPECTRUM.                                          !
  !        ----------------------------                                          !
  DO M = 1,nFre
     DO K = 1,nDir
        DO iCell = 1, nCells
           N0(iCell,K,M) = ET(iCell,M) * ST(iCell,K)
        END DO
     END DO
  END DO
  
  DO M = 1,nFre
     DO K = 1,nDir
        N0(0,K,M) = 0.0_rk
     END DO
  END DO
     
  !print*, "ST max min  sum", maxval(ST), minval(ST), sum(ST)
  !print*, "ET max min  sum", maxval(ET), minval(ET), sum(ET)
  !print*, "FL3 max min sum", maxval(N0), minval(N0), sum(N0)
   
END SUBROUTINE cold_start_wam


SUBROUTINE FETCH_LAW (FETCH, FPMAX, U10, ALPHJ, FP)
! ---------------------------------------------------------------------------- !
!   FETCH_LAW - COMPUTE JONSWAP PARAMETERS FROM FETCH LAW.                     !
!   COPY from WAM Cycle 6                                                      !
! ---------------------------------------------------------------------------- !

REAL(real_kind),    INTENT(IN)  :: FETCH         !! FETCH TO BE USED (METRES).
REAL(real_kind),    INTENT(IN)  :: FPMAX         !! MAXIMUM PEAK FREQUENCY (HERTZ).
REAL(real_kind),    INTENT(IN)  :: U10(:)        !! MODULUS OF WIND VELOCITY [M/S].
REAL(real_kind),    INTENT(OUT) :: ALPHJ(:)      !! JONSWAP ALPHA.
REAL(real_kind),    INTENT(OUT) :: FP(:)         !! JONSWAP PEAK FREQUENCY [HERTZ].

REAL(real_kind), PARAMETER :: A = 2.84,  D = -(3./10.) !! PEAKFREQUENCY FETCH LAW CONSTANTS
REAL(real_kind), PARAMETER :: B = 0.033, E = 2./3.     !! ALPHA-PEAKFREQUENCY LAW CONSTANTS

REAL(real_kind) :: UG(SIZE(U10))

! COMPUTE VALUES FROM FETCH LAWS.                                     

WHERE (U10 .GT. 0.1E-08)
   UG = g/U10
   FP = MAX(0.13, A*((g*FETCH)/(U10**2))**D)
   FP = MIN(FP, FPMAX/UG)
   ALPHJ = MAX(0.0081, B * FP**E)
   FP = FP*UG
ELSEWHERE
   ALPHJ = 0.0081
   FP = FPMAX
END WHERE

END SUBROUTINE FETCH_LAW


SUBROUTINE JONSWAP (FR, ALPHAJ, GAMMA, SA, SB, FP, ET)
! ---------------------------------------------------------------------------- !
!   JONSWAP - ROUTINE TO COMPUTE THE 1-D JONSWAP SPECTRUM.                     !
!   COPY from WAM Cycle 6                                                      !
! ---------------------------------------------------------------------------- !
REAL(real_kind), INTENT(IN)  :: FR(:)       !! FREQUENCIES.[Hz]
REAL(real_kind), INTENT(IN)  :: ALPHAJ(:)   !! OVERALL ENERGY LEVEL OF JONSWAP SPECTRA.
REAL(real_kind), INTENT(IN)  :: GAMMA       !! OVERSHOOT FACTOR.
REAL(real_kind), INTENT(IN)  :: SA          !! LEFT PEAK WIDTH.
REAL(real_kind), INTENT(IN)  :: SB          !! RIGHT PEAK WIDTH.
REAL(real_kind), INTENT(IN)  :: FP(:)       !! PEAK FREQUENCIES.
REAL(real_kind), INTENT(OUT) :: ET(:,:)     !! JONSWAP SPECTRA.

INTEGER :: I, M
REAL(real_kind)    :: FRH, LOG_GAMMA
REAL(real_kind)    :: G2ZPI4FRH5M(1:SIZE(FR))
REAL(real_kind)    :: ARG, SAA, SBB

G2ZPI4FRH5M(:) = G**2/pi2**4*FR(:)**(-5)
LOG_GAMMA = LOG(GAMMA)
SAA = 0.5/(SA*SA)
SBB = 0.5/(SB*SB)

DO M = 1,SIZE(FR)
   FRH = FR(M)
   DO I = 1,SIZE(FP)
      ARG = 1.25*(FP(I)/FRH)**4
      IF (ARG.LT.50.) THEN
         ET(I,M) = ALPHAJ(I)*G2ZPI4FRH5M(M)*EXP(-ARG)
      ELSE
         ET(I,M) = 0.
         CYCLE
      END IF
      IF (FRH.GT.FP(I)) THEN
         ARG = SBB*(FRH/FP(I)-1.)**2
      ELSE
         ARG = SAA*(FRH/FP(I)-1.)**2
      ENDIF
      IF (ARG.LT.99.) ET(I,M) = ET(I,M)*EXP(LOG_GAMMA*EXP(-ARG))
   END DO
END DO

END SUBROUTINE JONSWAP


end module initial_mod
