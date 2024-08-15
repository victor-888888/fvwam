! Purpose : Discretization of frequency and direction, and provide some parameters.
! Version  Date        Author        Description
! V1.0     2021-08-12  GAO Yuanyong  Fisrt Version

module fre_dir_mod
  
  use const_mod
  use params_mod

  implicit none
  
  ! Frequency

  !! set by namelist
  !integer          :: nFre     = 25     ! Size of frequency
  !real(real_kind)  :: freMin   = 0.04177248   ! Minimum frequency, Hz
  
  real(real_kind)  :: ratioFre = 1.1    ! Ratio of frequency, usually equal to 1.1
  real(real_kind), allocatable :: fre(:)      ! Frequency, Hz or 1/s
  real(real_kind), allocatable :: sigma(:)    ! Frequency, Radian, sigma = 2*pi/fre
  real(real_kind), allocatable :: delFre(:)   ! delta Fre, Hz
  real(real_kind), allocatable :: delSigma(:) ! delta sigma, Radian

  ! Direaction
  real(real_kind)              :: delTheta    ! Delta theta, Radian
  real(real_kind), allocatable :: theta(:)    ! Direction, Radian, True north is 0, clockwise, 0 to 2Pi
  real(real_kind), allocatable :: sinTheta(:) ! Sin of direction
  real(real_kind), allocatable :: cosTheta(:) ! Cos of direction

  integer, allocatable :: KPM(:,:)  ! Direction neightbour index, k+1(Plus), k, k-1(Minus)

  ! FREQUENCY DIRECTION INTERVALS AND AREAS
  REAL(real_kind),ALLOCATABLE :: DFIM(:)      !! FREQUENCY INTERVAL*DIRECTION INTER.
  REAL(real_kind),ALLOCATABLE :: DF_FR(:)   !! DF*FR [HZ*RAD].
  REAL(real_kind),ALLOCATABLE :: DF_FR2(:)   !! DF*FR*FR [HZ**HZ*RAD].
  REAL(real_kind),ALLOCATABLE :: DFIM_FR(:)   !! DFIM*FR [HZ*HZ*RAD].
  REAL(real_kind),ALLOCATABLE :: DFIM_FR2(:)  !! DFIM*FR**2  [HZ**3*RAD].
  REAL(real_kind),ALLOCATABLE :: FR5(:)       !! FR(M)**5
  REAL(real_kind),ALLOCATABLE :: FRM5(:)      !! 1./FR(M)**5
  REAL(real_kind),ALLOCATABLE :: RHOWG_DFIM(:)!! ROWATER*G*DELTH*LOG(CO)
  REAL(real_kind),ALLOCATABLE :: DFIMOFR(:)      

  ! TAIL FACTOR
  INTEGER, PARAMETER :: EX_TAIL = -5   !! TAIL FREQUENCY EXPONENT.
  REAL(real_kind)    :: MO_TAIL        !! MOMENT  O TAIL FACTOR.
  REAL(real_kind)     :: MM1_TAIL       !! MOMENT -1 TAIL FACTOR.
  REAL(real_kind)     :: MP1_TAIL       !! MOMENT +1 TAIL FACTOR.
  REAL(real_kind)     :: MP2_TAIL       !! MOMENT +2 TAIL FACTOR.

  contains

  Subroutine fre_dir_init()
      
      call generate_frequency() 
      
      call generate_direction()

      ! COMPUTATION FREQUENCY DIRECTION INTERVALS AND AREAS
      IF (.NOT.ALLOCATED (DFIM))       ALLOCATE (DFIM(1:nFre))       
      IF (.NOT.ALLOCATED (DFIMOFR))    ALLOCATE (DFIMOFR(1:nFre))    
      IF (.NOT.ALLOCATED (DF_FR))      ALLOCATE (DF_FR(1:nFre))      !! DF*FR.
      IF (.NOT.ALLOCATED (DF_FR2))     ALLOCATE (DF_FR2(1:nFre))     !! DF*FR*FR.
      IF (.NOT.ALLOCATED (FR5 ))       ALLOCATE (FR5 (nFre))       !! FR**5.
      IF (.NOT.ALLOCATED (FRM5))       ALLOCATE (FRM5(nFre))       !! 1. / FR**5.
      IF (.NOT.ALLOCATED (DFIM_FR) )   ALLOCATE (DFIM_FR(1:nFre))    !! DFIM*FR.
      IF (.NOT.ALLOCATED (DFIM_FR2))   ALLOCATE (DFIM_FR2(1:nFre))   !! DFIM*FR**2.
      IF (.NOT.ALLOCATED (RHOWG_DFIM)) ALLOCATE (RHOWG_DFIM(1:nFre))
      
      DFIM(:)     = delFre(:)*delTheta  !! MO  INTEGRATION WEIGHTS.
      DFIMOFR(:)  = DFIM(:)/Fre(:)      !! M-1 INTEGRATION WEIGHTS.
      
      DF_FR(:)  = delFre(:)*Fre(:)
      DF_FR2(:) = DF_FR(:)*Fre(:)
      
      DFIM_FR(:)  = DF_FR(:)*delTheta                       !! M+1 INTEGRATION WEIGHTS.
      DFIM_FR2(:) = DF_FR2(:)*delTheta                      !! M+2 INTEGRATION WEIGHTS.
      RHOWG_DFIM(:)=rho_sea*G*delTheta*LOG(ratioFre)*Fre(:) !! MOMENTUM AND ENERGY FLUX WEIGHTS.
      RHOWG_DFIM(1)=0.5*RHOWG_DFIM(1)                       !! TRAPEZOIDAL INTEGRATION
      RHOWG_DFIM(nFre)=0.5*RHOWG_DFIM(nFre)                   !! WITH CHANGE OF VARIABLE x=LOG(f)
                                                            !! HENCE df = f dx,
                                                            !! dx=LOG(FR(n+1))-LOG(FR(n)) = LOG(CO)
      
      ! COMPUTATION TAIL FACTOR
      MO_TAIL  = -delTheta/ REAL(EX_TAIL+1)*Fre(nFre)       !! MO  TAIL FACTOR.
      MM1_TAIL = -delTHeta/ REAL(EX_TAIL)                   !! M-1 TAIL FACTOR.
      MP1_TAIL = -delTheta/ REAL(EX_TAIL+2)*Fre(nFre)**2    !! M+1 TAIL FACTOR.
      MP2_TAIL = -delTheta/ REAL(EX_TAIL+3)*Fre(nFre)**3    !! M+2 TAIL FACTOR.

  end subroutine fre_dir_init

  Subroutine generate_frequency()
      implicit none
      integer m
      real(real_kind) :: CO1

      If (.not. allocated(fre     )) allocate(fre(nFre     )) 
      If (.not. allocated(sigma   )) allocate(sigma(nFre   )) 
      If (.not. allocated(delSigma)) allocate(delSigma(nFre)) 
      If (.not. allocated(delFre  )) allocate(delFre(nFre  )) 

      fre(1)   = freMin
      sigma(1) = 2.0*Pi*fre(1)

      do m = 2, nFre
         fre(m)   = ratioFre*fre(m-1)
         sigma(m) = 2.0*Pi*fre(m)
      end do
     
      CO1 = 0.5*(ratioFre-1.0)
      delFre(1)   = CO1*Fre(1)
      delSigma(1) = 2.0*Pi*delFre(1)
      do m = 2, nFre-1
         delFre(m)   = CO1*(Fre(m) + Fre(m-1))
         delSigma(m) = 2.0*Pi*delFre(m)
      end do
      delFre(nFre)   = CO1*Fre(nFre-1) 
      delSigma(nFre) = 2.0*Pi*delFre(nFre)

  end subroutine generate_frequency

  Subroutine generate_direction()
      implicit none
      integer k
      If (.not. allocated(theta   )) allocate(theta(nDir   ))
      If (.not. allocated(sinTheta)) allocate(sinTheta(nDir))
      If (.not. allocated(cosTheta)) allocate(cosTheta(nDir))
      If (.not. allocated(KPM)) allocate(KPM(nDir, -1:1))
      
      delTheta = 2.0*PI/REAL(nDir) ! delta theta
      
      do k = 1, nDir
         theta(k) = REAL(k-1)*delTheta + 0.5*delTheta
         sinTheta(k) = sin(theta(k))
         cosTheta(k) = cos(theta(k))
      end do

      do k = 1, nDir
         KPM(K,-1) = K-1
         if (KPM(k,-1) .lt. 1   ) KPM(k,-1) = nDir
  
         KPM(k,0) = k

         KPM(k,1) = k+1
         if (KPM(k,1)  .gt. nDir) KPM(k,1) = 1
      end do

  end subroutine generate_direction
  
 
end module

