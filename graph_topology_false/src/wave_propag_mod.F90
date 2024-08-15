module wave_propag_mod
  use const_mod
  use params_mod
  use mesh_mod
  use fre_dir_mod
  use mask_mod

  implicit none

  real(real_kind), allocatable :: Cg(:,:)       ! Group velocity, m/s  (nDepth, nFre)
  real(real_kind), allocatable :: cgEdgeX(:,:,:)  ! X-direction Group velocity on Edge, m/s  (nDepth, nDir, nFre)
  real(real_kind), allocatable :: cgEdgeY(:,:,:)  ! Y-Direction Group velocity, m/s  (nDepth, nDir, nFre)

  real(real_kind), allocatable :: dotTheta(:,:,:) ! Theta dot (nCells, nDir, nFre)
  
  real(real_kind), allocatable :: COEF_geog(:,:)  ! (maxEdges, nCells)
  real(real_kind), allocatable :: COEF_theta(:,:,:,:) ! (nCells, nFre, nDir, -1:1) k+1, k, k-1

  real(real_kind), allocatable :: kwave(:,:)         ! wave number(nDepth, nFre)
  
  integer, allocatable :: cellIdxTab(:)         ! cell index table
  integer, allocatable :: edgeIdxTab(:)         ! edge index table
  
  real(real_kind), allocatable :: depthBin(:)     
  
  real(real_kind) :: depthMax

  real(real_kind), allocatable :: dHdx(:)         ! Depth gradient, x-direction (nCells)
  real(real_kind), allocatable :: dHdy(:)         ! Depth gradient, y-direction (nCells)
  real(real_kind), allocatable :: dHdxSmooth(:)         ! Depth gradient, x-direction (nCells)
  real(real_kind), allocatable :: dHdySmooth(:)         ! Depth gradient, y-direction (nCells)
  real(real_kind), allocatable :: bottomDepthEdge(:) 
  integer, allocatable :: cellIdxTabSmooth(:)         ! cell index table
  real(real_kind), allocatable :: bottomDepthSmooth(:) 
 
  contains

  subroutine propag_init()
    implicit none
    allocate(Cg(nDepth, nFre))
    allocate(cgEdgeX(nDepth, nDir, nFre))
    allocate(cgEdgeY(nDepth, nDir, nFre))

    allocate(dotTheta(nCells, nDir, nFre))

    allocate(COEF_geog(maxEdges, nCells))
    allocate(COEF_theta(nCells, nDir, nFre, -1:1))

    allocate(kwave(nDepth, nFre))
    
    allocate(cellIdxTab(nCells))
    allocate(edgeIdxTab(nEdges))
    
    allocate(depthBin(nDepth))

    allocate(dHdx(nCells))
    allocate(dHdy(nCells))
    allocate(dHdxSmooth(mpi_num_ncells_halo))
    allocate(dHdySmooth(mpi_num_ncells_halo))
    allocate(bottomDepthEdge(nEdges))

    allocate(cellIdxTabSmooth(nCells))
    allocate(bottomDepthSmooth(mpi_num_ncells_halo))

    Cg        = 0.0_rk
    COEF_geog = 0.0_rk

    call prepare_propag()

  end subroutine propag_init

  subroutine prepare_propag()
      implicit none
      integer        :: iCell, iEdge, k, m
      real(real_kind) :: kd, kd2, temp
     
      integer nCOC, coc, j, i, eoc
      real(real_kind) :: bt

      if (smooth_depth_num .lt. 1) then
         ! only smooth for cellIdxTabSmooth
         bottomDepthSmooth(:) = bottomDepth(:)
         if(mpi_procs>1) then
           call mpi_data_exchange_oncell(bottomDepthSmooth(1:mpi_num_ncells_halo))
         end if

         do i =1, 2
             do iCell = 1, nCells
                 bt = bottomDepthSmooth(iCell)
                 nCOC = 1
                 do j=1, nEdgesOnCell(iCell)
                    coc = cellsOnCell(j,iCell)
                    if (coc .ne. 0) then 
                        nCOC = nCOC + 1
                        bt = bt + bottomDepthSmooth(coc)
                    end if
                 end do
                 bottomDepthSmooth(iCell) = bt/nCOC
             end do
             !bottomDepth updated in the it, so we need to exchange data
             if(mpi_procs>1) then
               call mpi_data_exchange_oncell(bottomDepthSmooth(1:mpi_num_ncells_halo))
             end if
         end do

      ! smooth bottomDepth 
      else
         if (mpi_rank == 0) then
            print*, 'Smoothing the bottom depth, smooth_depth_num = ', smooth_depth_num
         end if

         if(mpi_procs>1) then
           call mpi_data_exchange_oncell(bottomDepth(1:mpi_num_ncells_halo))
         end if

         do i =1, smooth_depth_num 
             do iCell = 1, nCells
                 bt = bottomDepth(iCell)
                 nCOC = 1
                 do j=1, nEdgesOnCell(iCell)
                    coc = cellsOnCell(j,iCell)
                    if (coc .ne. 0) then 
                        nCOC = nCOC + 1
                        bt = bt + bottomDepth(coc)
                    end if
                 end do
                 bottomDepth(iCell) = bt/nCOC
             end do
             !bottomDepth updated in the it, so we need to exchange data
             if(mpi_procs>1) then
               call mpi_data_exchange_oncell(bottomDepth(1:mpi_num_ncells_halo))
             end if
         end do
         bottomDepthSmooth(:) = bottomDepth(:)
      end if

      ! depth bins
      do i = 1, nDepth
         depthBin(i) = depthMin*depthRatio**(i-1)
      end do
      depthMax = depthBin(nDepth)
      
      do iCell = 1, nCells
         if (bottomDepth(iCell) .gt. depthMax) bottomDepth(iCell) = depthMax
      end do


      ! index table
      do iCell = 1, nCells
         cellIdxTab(iCell) = nint(log(bottomDepth(iCell)/depthMin)/log(depthRatio)+1.0)
         cellIdxTab(iCell) = max(cellIdxTab(iCell), 1)
         cellIdxTab(iCell) = min(cellIdxTab(iCell), nDepth)
         
         cellIdxTabSmooth(iCell) = nint(log(bottomDepthSmooth(iCell)/depthMin)/log(depthRatio)+1.0)
         cellIdxTabSmooth(iCell) = max(cellIdxTabSmooth(iCell), 1)
         cellIdxTabSmooth(iCell) = min(cellIdxTabSmooth(iCell), nDepth)
      end do

      if (mpi_procs>1) then
         call mpi_data_exchange_oncell(bottomDepth(1:mpi_num_ncells_halo))
      end if
      call scalar_c2e_interp_operator(bottomDepth(0:mpi_num_ncells_halo), bottomDepthEdge)
      do iEdge = 1, nEdges
         edgeIdxTab(iEdge) = nint(log(bottomDepthEdge(iEdge)/depthMin)/log(depthRatio)+1.0)
         edgeIdxTab(iEdge) = max(edgeIdxTab(iEdge), 1)
         edgeIdxTab(iEdge) = min(edgeIdxTab(iEdge), nDepth)
      end do
      deallocate(bottomDepthEdge)
  
      
      ! Compute wave number
      do m = 1, nFre
         do i = 1, nDepth
               kwave(i, m) =  wave_numb(sigma(m), depthBin(i))
         end do
      end do

      ! Compute group velocity
      do m = 1, nFre
         do i = 1, nDepth
            kd = kwave(i,m)*depthBin(i)
            if(kd > 10.) then ! Deep water 
              Cg(i, m) = g/(2.0*sigma(m))
            else              ! Shallow water 
              !Cg(i, m) = sigma(m)/kwave(i,m)*(0.5+kd/sinh(2.0*kd))
              Cg(i, m) = 0.5*sqrt(g*tanh(kd)/kwave(i,m))*(1.0+2.0*kd/sinh(2.0*kd))
            end if
         end do
      end do

      do i = 1, nDepth
         do k = 1, nDir
            do m = 1, nFre
               cgEdgeX(i, k, m) = Cg(i, m)*sinTheta(k)  ! Theta, North = 0, clockwise
               cgEdgeY(i, k, m) = Cg(i, m)*cosTheta(k)
            end do
         end do
      end do

      ! Compute depth gradient
      call calc_depth_gradient(bottomDepth, dHdx, dHdy)

      ! Compute depth part of dot theta
      do m = 1, nFre
         do k = 1, nDir
            do iCell = 1, nCells
               !i = cellIdxTab(iCell)
               i = cellIdxTabSmooth(iCell)
               kd = kwave(i,m)*depthBin(i)
               if (kd .le. 10.0) then
                  temp = sinTheta(k)*dHdy(iCell) - cosTheta(k)*dHdx(iCell)
                  kd2 = 2.0*kd
                  dotTheta(iCell, k, m) = sigma(m)/sinh(kd2) * temp
               else
                  dotTheta(iCell, k, m) = 0.0
               end if
            end do
         end do
      end do

      ! Compute great circle part of dot theta
      do m = 1, nFre
         do k = 1, nDir
            do iCell = 1, nCells
               i = cellIdxTab(iCell)
               dotTheta(iCell, k, m) = dotTheta(iCell, k, m) - Cg(i, m)*tanLat(iCell)*cosTheta(k)/radius
            end do
         end do
      end do

      call prepare_coef_geog()   ! COEF_geog
      call prepare_coef_theta()  ! COEF_theta
  
  end subroutine prepare_propag
  
  subroutine calc_tend_geog(N_old, tend)
    implicit none
    real(real_kind), intent(in)  :: N_old(0:,:,:)   ! wave actioon, old time
    real(real_kind), intent(out) :: tend(:,:,:)   

    integer         :: iCell, k, m, j, iE, nEOC, upwindCell, eoc
    real(real_kind) :: temp, cgEdgeNorm
    
    ! geog tend
    do m = 1, nFre
       do k = 1, nDir
          do iCell = 1, nCells
             temp = 0.0_rk
             do j = 1, nEdgesOnCell(iCell)
                  eoc = edgesOnCell(j, iCell)
                  iE = edgeIdxTab(eoc)
                  cgEdgeNorm = cgEdgeX(iE, k, m)*angleEdgeCos(eoc) + cgEdgeY(iE, k, m)*angleEdgeSin(eoc)
                  if ( nSignEdge(j, iCell)*cgEdgeNorm .gt. 0.0 ) then ! outflow, upwind source is iCell, sink is coc
                     upwindCell  = iCell
                  else
                     upwindCell  = cellsOnCell(j, iCell) ! coc
                  end if
                  temp = temp + N_old(upwindCell, k, m)*COEF_geog(j, iCell)*cgEdgeNorm
             end do
             tend(iCell, k, m) = temp*maskBdy(iCell)
          end do
       end do
    end do

  end subroutine calc_tend_geog
  
  subroutine calc_tend_theta(N_old, tend)
    implicit none
    real(real_kind), intent(in)  :: N_old(0:,:,:)   ! wave actioon, old time
    real(real_kind), intent(out) :: tend(:,:,:)   

    integer         :: iCell, k, m, kp, km, i

    ! theta
    do m = 1, nFre
       do k = 1, nDir
          do iCell = 1, nCells
             kp = KPM(k,  1) ! k+1, plus index
             km = KPM(k, -1) ! k-1, minus index
             i = cellIdxTab(iCell)
             tend(iCell, k, m)  =   COEF_theta(i, k, m,  1)*N_old(iCell, kp, m)  &
                                  + COEF_theta(i, k, m, -1)*N_old(iCell, km, m)
             tend(iCell, k, m)  = tend(iCell, k, m)*maskBdy(iCell)
          end do
       end do
    end do

  end subroutine calc_tend_theta
  
  subroutine prepare_coef_geog()
      implicit none
      integer         :: iCell, j, nEOC, eoc
      real(real_kind) :: dAreaCell
     
      do iCell = 1, nCells
         dAreaCell = 1.0/areaCell(iCell)
         nEOC = nEdgesOnCell(iCell)
         do j = 1, nEOC
            eoc = edgesOnCell(j, iCell)
            COEF_geog(j, iCell) = -nSignEdge(j, iCell)*dvEdge(eoc)*dAreaCell
         end do
      end do

  end subroutine prepare_coef_geog

  subroutine prepare_coef_theta()
      implicit none
      integer         :: iCell, k, m, kp, km
      real(real_kind) :: arg 

      arg = 0.5_rk/delTheta

      do m = 1, nFre
         do k = 1, nDir 
            do iCell = 1, nCells
               kp = KPM(k, 1) ! k + 1
               km = KPM(k,-1) ! k - 1
               COEF_theta(iCell, k, m,  1) = -arg*dotTheta(iCell, kp, m) ! k+1
               COEF_theta(iCell, k, m, -1) =  arg*dotTheta(iCell, km, m) ! k-1
               COEF_theta(iCell, k, m,  0) =  0.0_rk                     ! k, no use
            end do
         end do
      end do
  end subroutine prepare_coef_theta
     
  subroutine calc_depth_gradient(depth, dHdx, dHdy)

      ! Calculate depth gradient.
      ! Establish a polar coordinate system in each cell, and convert the
      ! depth gradient in the rectangular coordinate system to the polar 
      ! coordinate system.
      ! Author : GAO Yuanyong  2021-08-20
  
      implicit none
      real(real_kind), intent(in)  :: depth(:)  ! Water depth, meter
      real(real_kind), intent(out) :: dHdx(:)   ! Depth gradient, x-direction 
      real(real_kind), intent(out) :: dHdy(:)   ! Depth gradient, y-direction

      integer         :: iCell, i, j, nCOC, coc
      real(real_kind) :: temp1, temp2

      dHdx(:) = dHdx_read(:)*maskBdy
      dHdy(:) = dHdy_read(:)*maskBdy
      
      if (smooth_depth_num .gt. 0) then

         dHdxSmooth(1:nCells) = dHdx(:)
         dHdySmooth(1:nCells) = dHdy(:)
         if (mpi_procs>1) then
            call mpi_data_exchange_oncell(dHdxSmooth(1:mpi_num_ncells_halo))
            call mpi_data_exchange_oncell(dHdySmooth(1:mpi_num_ncells_halo))
         end if

         do i =1, smooth_depth_num
             do iCell = 1, nCells
                 temp1= dHdxSmooth(iCell)
                 temp2= dHdySmooth(iCell)
                 nCOC = 1
                 do j=1, nEdgesOnCell(iCell)
                    coc = cellsOnCell(j,iCell)
                    if (coc .ne. 0) then 
                        nCOC = nCOC + 1
                        temp1 = temp1 + dHdxSmooth(coc)
                        temp2 = temp2 + dHdySmooth(coc)
                    end if
                 end do
                 dHdxSmooth(iCell) = temp1/nCOC
                 dHdySmooth(iCell) = temp2/nCOC
             end do
             !bottomDepth updated in the it, so we need to exchange data
             if (mpi_procs>1) then
                call mpi_data_exchange_oncell(dHdxSmooth(1:mpi_num_ncells_halo))
                call mpi_data_exchange_oncell(dHdySmooth(1:mpi_num_ncells_halo))
             end if
         end do
         dHdx = dHdxSmooth(1:nCells)*maskBdy
         dHdy = dHdySmooth(1:nCells)*maskBdy
     end if
  
  end subroutine calc_depth_gradient

  subroutine scalar_c2e_interp_operator(f_cell, f_edge)
     
     real(real_kind), intent(in ) :: f_cell(0:mpi_num_ncells_halo) 
     real(real_kind), intent(out) :: f_edge(:)
     integer iEdge
     do iEdge = 1, nEdges
           f_edge(iEdge) = 0.5_rk * (f_cell(cellsOnEdge(1,iEdge)) +  f_cell(cellsOnEdge(2,iEdge)))
     end do
      
  end subroutine scalar_c2e_interp_operator

  subroutine update_state(dt, tend, N_old, N_new)
    implicit none 
    real(real_kind), intent(in)  :: dt
    real(real_kind), intent(in)  :: tend(:,:,:)
    real(real_kind), intent(in)  :: N_old(0:,:,:)
    real(real_kind), intent(out) :: N_new(0:,:,:)
    integer :: iCell, k, m
     
    do k = 1, nDir
       do m = 1, nFre
          do iCell = 1, nCells
             N_new(iCell, k, m) = N_old(iCell, k, m) + dt*tend(iCell, k, m)
          end do
       end do
    end do

  end subroutine update_state

!  real(real_kind) function wave_numb(sigma ,depth)
!      ! calculate wave number, copy from WW3
!      ! GAO YY, 2021-08-12
!      implicit none
!      real(real_kind), intent(in)  ::  sigma  ! radian frequency (rad)
!      real(real_kind), intent(in)  ::  depth  ! water depth (m)
!      real(real_kind) :: x,xx,y,omega
!
!      if(depth<=0.0 .or. sigma<= 0.) then
!         wave_numb = -10.
!      else
!        omega = sigma**2/g   ! g, gravitational acceleration (m/s^2)
!        y     = omega*depth
!        xx    = y*(y+1./(1.+y*(0.66667+y*(0.35550+y*(0.16084+y*(0.06320+y* &
!                  (0.02174+y*(0.00654+y*(0.00171+y*(0.00039+y*0.00011))))))))))
!        x     = sqrt(xx)
!        wave_numb = x/depth
!      end if
!      return
!  end function

  REAL(real_kind) FUNCTION wave_numb (OM, BETA)
 
  !  FUNCTION TO COMPUTE WAVE NUMBER(Copy from WAM model)
  !  NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW WATER.    
 
 
  REAL, INTENT(IN) :: OM    !! CIRCULAR FREQUENCY.
  REAL, INTENT(IN) :: BETA  !! WATER DEPTH.
 
  !     LOCAL VARIABLES.                                                         !
  REAL, PARAMETER :: EBS = 0.0001  !! RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
  REAL, PARAMETER :: DKMAX = 40.        !! MAXIMUM VALUE OF DEPTH*WAVENUMBER. 
 
  REAL :: AKP, BO, TH, STH
 
  !     1. START WITH MAXIMUM FROM DEEP AND EXTREM SHALLOW WATER WAVE NUMBER.    !
  wave_numb   = MAX ( OM**2/(4.*G), OM/(2.*SQRT(G*BETA)) )
  !                                                                              !
  !     2. ITERATION LOOP.                                                       !
 
  AKP = 10000.
  DO WHILE (ABS(AKP-wave_numb).GT.EBS*wave_numb)
     BO = BETA*wave_numb
     IF (BO.GT.DKMAX) THEN
        wave_numb = OM**2/G
        EXIT
     ELSE
        AKP = wave_numb
        TH = G*wave_numb*TANH(BO)
        STH = SQRT(TH)
        wave_numb = wave_numb+(OM-STH)*STH*2./(TH/wave_numb+G*BO/COSH(BO)**2)
     END IF
  END DO

 END FUNCTION

end module


