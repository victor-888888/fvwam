module forcing_mod

  use const_mod
  use mesh_mod
  use log_mod
  use read_forcing_mod
  use time_mod
  use params_mod
  use string
  use mod_mpi_interfaces

  implicit none
  private

  public forcing_init
  public update_atm_forcing
  public update_atm_forcing_openacc
  public update_forcing_NWP
  public update_forcing_TC
  public interp_wind_to_steptime
  public dt_src

  public forcing_type
  public atm_forcing_type
  public TC_forcing_type
  public forcing
  public atm_forcing
  public TC_forcing

  type atm_forcing_type
    real(real_kind), allocatable :: u10(:)    ! zonal wind, m/s
    real(real_kind), allocatable :: v10(:)    ! meridian wind, m/s
  end type atm_forcing_type

  type TC_forcing_type
    real(real_kind), allocatable :: distCell   (:) ! Distance from TC center on cell, meter
    real(real_kind), allocatable :: uc         (:) ! TC x direction combination wind on cell, c2*u + c1*ut, m/s 
    real(real_kind), allocatable :: vc         (:) ! TC y direction combination wind on cell, c2*v + c1*vt, m/s 
    real(real_kind), allocatable :: latTC      (:) ! Latitude and longitude of TC center, radian
    real(real_kind), allocatable :: lonTC      (:) ! Latitude and longitude of TC center, radian
    real(real_kind), allocatable :: utTC       (:) ! Tansfer velocity of TC center, m/s
    real(real_kind), allocatable :: vtTC       (:) ! Tansfer velocity of TC center, m/s
    real(real_kind), allocatable :: presTC     (:) ! Pressure On TC center, Pa
    real(real_kind), allocatable :: rmwTC      (:) ! Maximum wind speed radius of typhoon, meter
  end type TC_forcing_type

  type forcing_type
    real(real_kind), allocatable :: u10Spd(:) ! total wind speed, m/s
    real(real_kind), allocatable :: u10Dir(:) ! wind direction, radian
    real(real_kind), allocatable :: u10(:)    ! zonal wind, m/s
    real(real_kind), allocatable :: v10(:)    ! meridian wind, m/s
    integer                      :: atm_time_size                  ! Time size of atm input file 
    real(real_kind)              :: atm_interval_seconds           ! Atmpshpere forcing interval
    integer                      :: atm_time_steps_from_update = 0 ! Integration steps after updating the atmospheric field
    real(real_kind)              :: atm_time_rate                  ! Time rate used for time interpolation
  end type forcing_type

  type(forcing_type)                  :: forcing
  type(atm_forcing_type), allocatable :: atm_forcing(:)
  type(TC_forcing_type)               :: TC_forcing

  integer :: atm_forcing_time_step = 0 
  real(real_kind) :: dt_src

contains

  subroutine forcing_init()
    implicit none
    character(10) time_value, time_units

    dt_src = dt * dt_src_ratio 

    if (atm_forcing_option .gt. 0) then
       time_value = split_string(atm_forcing_interval, ' ', 1)
       time_units = split_string(atm_forcing_interval, ' ', 2)
       read(time_value, *) forcing%atm_interval_seconds
       select case (time_units)
       case ('days')
         forcing%atm_interval_seconds = forcing%atm_interval_seconds * 86400
       case ('hours')
         forcing%atm_interval_seconds = forcing%atm_interval_seconds * 3600
       case ('minutes')
         forcing%atm_interval_seconds = forcing%atm_interval_seconds * 60
       case ('seconds')
         forcing%atm_interval_seconds = forcing%atm_interval_seconds
       case default
         if(mpi_rank==0) then
         call log_error('Invalid forcing interval ' // trim(atm_forcing_interval) // '!')
         end if
       end select

       call time_add_alert('atm_forcing_update', seconds=forcing%atm_interval_seconds)
       forcing%atm_time_size = int(run_duration_time%total_seconds()/forcing%atm_interval_seconds) + 1

       allocate(forcing%u10Spd(nCells))
       allocate(forcing%u10Dir(nCells))
       allocate(forcing%u10(nCells))
       allocate(forcing%v10(nCells))

    end if
   
    if (.not. allocated(atm_forcing)) then
       allocate(atm_forcing(1:2))
       call allocate_atm_forcing(atm_forcing(1))
       call allocate_atm_forcing(atm_forcing(2))
    end if



    if (atm_forcing_option .eq. 2)then
    
       if (.not. allocated(TC_forcing%distCell   ))  allocate(TC_forcing%distCell   (nCells))
       if (.not. allocated(TC_forcing%uc     ))  allocate(TC_forcing%uc     (nCells))
       if (.not. allocated(TC_forcing%vc     ))  allocate(TC_forcing%vc     (nCells))
       if (.not. allocated(TC_forcing%lonTC      ))  allocate(TC_forcing%lonTC      (forcing%atm_time_size))
       if (.not. allocated(TC_forcing%latTC      ))  allocate(TC_forcing%latTC      (forcing%atm_time_size))
       if (.not. allocated(TC_forcing%utTC       ))  allocate(TC_forcing%utTC       (forcing%atm_time_size))
       if (.not. allocated(TC_forcing%vtTC       ))  allocate(TC_forcing%vtTC       (forcing%atm_time_size))
       if (.not. allocated(TC_forcing%presTC     ))  allocate(TC_forcing%presTC     (forcing%atm_time_size))
       if (.not. allocated(TC_forcing%rmwTC      ))  allocate(TC_forcing%rmwTC      (forcing%atm_time_size))
       if (mpi_rank==0) then
       call log_notice('Using TC empirical model to force the model!')
       call read_tc_key_parameters(forcing%atm_time_size,  &   
                                   TC_forcing%lonTC,       &
                                   TC_forcing%latTC,       &
                                   TC_forcing%rmwTC,       &
                                   TC_forcing%presTC,      &
                                   TC_forcing%utTC,        &
                                   TC_forcing%vtTC         )

      end if


      if(mpi_procs>1) then

        call mpi_bcast_key_parameters( TC_forcing%lonTC,       &
                                       TC_forcing%latTC,       &
                                       TC_forcing%rmwTC,       &
                                       TC_forcing%presTC,      &
                                       TC_forcing%utTC,        &
                                       TC_forcing%vtTC         )
      end if
    end if

    call time_add_alert('atm_forcing_update', seconds=forcing%atm_interval_seconds)
    call time_add_alert('integrate_source',   seconds=dt_src)
    call update_atm_forcing()

  end subroutine forcing_init

  subroutine interp_wind_to_steptime()
    ! Interpolate meteorological data to each integration step.
    
    integer i
    forcing%atm_time_rate = forcing%atm_time_steps_from_update*dt_src/forcing%atm_interval_seconds

#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("atm_forcing1_u10.nc"),atm_forcing(1)%u10(1:nCells))
    call mpi_test_output(TRIM("atm_forcing1_v10.nc"),atm_forcing(1)%v10(1:nCells))
    call mpi_test_output(TRIM("atm_forcing2_u10.nc"),atm_forcing(2)%u10(1:nCells))
    call mpi_test_output(TRIM("atm_forcing2_v10.nc"),atm_forcing(2)%v10(1:nCells))
    end if
#endif
    do i = 1, nCells
       forcing%u10(i) = atm_forcing(1)%u10(i)*(1.0_rk - forcing%atm_time_rate) + atm_forcing(2)%u10(i)*forcing%atm_time_rate
       forcing%v10(i) = atm_forcing(1)%v10(i)*(1.0_rk - forcing%atm_time_rate) + atm_forcing(2)%v10(i)*forcing%atm_time_rate
    end do

#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("_forcing_u10.nc"),forcing%u10(1:nCells))
    call mpi_test_output(TRIM("_forcing_v10.nc"),forcing%v10(1:nCells))
    end if
#endif
    call transfer_uv_to_spd_dir_openacc(forcing%u10, forcing%v10, forcing%u10Spd, forcing%u10Dir)
#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("_forcing_u10Spd.nc"),forcing%u10Spd(1:nCells))
    call mpi_test_output(TRIM("_forcing_u10Dir.nc"),forcing%u10Dir(1:nCells))
    end if
#endif

  end subroutine interp_wind_to_steptime
  
  subroutine update_atm_forcing( )
    implicit none
    real(real_kind) Cd ! wind drag coefficient
    integer i

    if(mpi_rank==0) then
    call log_notice('Update forcing field' )
    end if
    atm_forcing_time_step = atm_forcing_time_step + 1
    forcing%atm_time_steps_from_update = 0             ! reset to 0

    select case (atm_forcing_option)
    case (1) ! Numerical Weather Prediction(NWP) model, such as ERA5, WRF
      call update_forcing_NWP(atm_forcing_time_step)
    
    case (2) ! TC wind model, such as Holland, Fujita model
      if (atm_forcing_time_step .gt. 1 .and. atm_forcing_time_step .lt. forcing%atm_time_size) then
         do i = 1, nCells
           atm_forcing(1)%u10(i) = atm_forcing(2)%u10(i)
           atm_forcing(1)%v10(i) = atm_forcing(2)%v10(i)
         end do
         call update_forcing_TC(atm_forcing_time_step+1 ,atm_forcing(2)%u10, atm_forcing(2)%v10)
      else if (atm_forcing_time_step .eq. 1) then
         call update_forcing_TC(atm_forcing_time_step   ,atm_forcing(1)%u10 ,atm_forcing(1)%u10)
         call update_forcing_TC(atm_forcing_time_step+1 ,atm_forcing(2)%v10 ,atm_forcing(2)%v10)
      else
         do i = 1, nCells
           atm_forcing(1)%u10(i) = atm_forcing(2)%u10(i)
           atm_forcing(1)%v10(i) = atm_forcing(2)%v10(i)
         end do
      end if
    case (0) ! Turn off atmosphere forcing
      continue
    case default
      if(mpi_rank==0) then
      call log_error('Invalid atm_forcing_option, 0(turn off), 1(NWP), 2(TC) !')
      end if
    end select

!    !!!! dont interpolate to time step
    do i = 1, nCells
       forcing%u10(i) = atm_forcing(1)%u10(i)
       forcing%v10(i) = atm_forcing(1)%v10(i)
    end do

    call transfer_uv_to_spd_dir(forcing%u10, forcing%v10, forcing%u10Spd, forcing%u10Dir)
  end subroutine update_atm_forcing


  subroutine update_forcing_NWP(time_step)
    integer,         intent(in)  :: time_step
    integer :: i
    allocate(mpi_real_buf_nCells(mpi_total_nCells))
    allocate(mpi_real_buf_nCells_2(mpi_total_nCells))
    if (time_step .gt. 1 .and. time_step .lt. forcing%atm_time_size) then
       do i=1,nCells
         atm_forcing(1)%u10(i) = atm_forcing(2)%u10(i)
         atm_forcing(1)%v10(i) = atm_forcing(2)%v10(i)
       end do
        
       if(mpi_procs>1) then
         if(mpi_rank==0) then
           call read_met_forcing_var(mpi_total_nCells, time_step+1, mpi_real_buf_nCells, mpi_real_buf_nCells_2)
           call mpi_input_exchange_root(mpi_real_buf_nCells,   atm_forcing(2)%u10, nCells, 'oncell')
           call mpi_input_exchange_root(mpi_real_buf_nCells_2, atm_forcing(2)%v10, nCells, 'oncell')
         else
           call mpi_input_exchange_nonroot(mpi_real_buf_nCells,   atm_forcing(2)%u10, nCells, 'oncell')
           call mpi_input_exchange_nonroot(mpi_real_buf_nCells_2, atm_forcing(2)%v10, nCells, 'oncell')
         end if
       else
         call read_met_forcing_var(nCells, time_step+1,atm_forcing(2)%u10,atm_forcing(2)%v10)
       end if
    else if (time_step .eq. 1) then
     
      if(mpi_procs>1) then
        if(mpi_rank==0) then
          call read_met_forcing_var(mpi_total_nCells, time_step, mpi_real_buf_nCells, mpi_real_buf_nCells_2)
          call mpi_input_exchange_root(mpi_real_buf_nCells,   atm_forcing(1)%u10, nCells, 'oncell')
          call mpi_input_exchange_root(mpi_real_buf_nCells_2, atm_forcing(1)%v10, nCells, 'oncell')
        else
           call mpi_input_exchange_nonroot(mpi_real_buf_nCells,   atm_forcing(1)%u10, nCells, 'oncell')
           call mpi_input_exchange_nonroot(mpi_real_buf_nCells_2, atm_forcing(1)%v10, nCells, 'oncell')
        end if
    
        if(mpi_rank==0) then
           call read_met_forcing_var(mpi_total_nCells, time_step+1, mpi_real_buf_nCells, mpi_real_buf_nCells_2)
           call mpi_input_exchange_root(mpi_real_buf_nCells,   atm_forcing(2)%u10, nCells, 'oncell')
           call mpi_input_exchange_root(mpi_real_buf_nCells_2, atm_forcing(2)%v10, nCells, 'oncell')
        else
           call mpi_input_exchange_nonroot(mpi_real_buf_nCells,   atm_forcing(2)%u10, nCells, 'oncell')
           call mpi_input_exchange_nonroot(mpi_real_buf_nCells_2, atm_forcing(2)%v10, nCells, 'oncell')
        end if
      else
         call read_met_forcing_var(nCells, time_step,atm_forcing(1)%u10,atm_forcing(1)%v10)
         call read_met_forcing_var(nCells, time_step+1,atm_forcing(2)%u10,atm_forcing(2)%v10)

      end if

    else
       do i=1,nCells
         atm_forcing(1)%u10(i) = atm_forcing(2)%u10(i)
         atm_forcing(1)%v10(i) = atm_forcing(2)%v10(i)
       end do
    end if

    deallocate(mpi_real_buf_nCells)
    deallocate(mpi_real_buf_nCells_2)
  end subroutine update_forcing_NWP

  subroutine update_atm_forcing_openacc()
    real(real_kind):: Cd ! wind drag coefficient
    integer :: i

    if(mpi_rank==0) then
      call log_notice('Update forcing field' )
    end if

    atm_forcing_time_step = atm_forcing_time_step + 1
    forcing%atm_time_steps_from_update = 0             ! reset to 0

    select case (atm_forcing_option)
    case (1) ! Numerical Weather Prediction(NWP) model, such as ERA5, WRF
      call update_forcing_NWP_openacc(atm_forcing_time_step)
    
    case (2) ! TC wind model, such as Holland, Fujita model
      if (atm_forcing_time_step .gt. 1 .and. atm_forcing_time_step .lt. forcing%atm_time_size) then
         do i = 1, nCells
           atm_forcing(1)%u10(i) = atm_forcing(2)%u10(i)
           atm_forcing(1)%v10(i) = atm_forcing(2)%v10(i)
         end do
         call update_forcing_TC_openacc(atm_forcing_time_step+1 ,atm_forcing(2)%u10, atm_forcing(2)%v10)
      else if (atm_forcing_time_step .eq. 1) then
         call update_forcing_TC_openacc(atm_forcing_time_step   ,atm_forcing(1)%u10 ,atm_forcing(1)%u10)
         call update_forcing_TC_openacc(atm_forcing_time_step+1 ,atm_forcing(2)%v10 ,atm_forcing(2)%v10)
      else
         do i = 1, nCells
           atm_forcing(1)%u10(i) = atm_forcing(2)%u10(i)
           atm_forcing(1)%v10(i) = atm_forcing(2)%v10(i)
         end do
      end if
    case default
      if(mpi_rank==0) then
      call log_error('Invalid atm_forcing_option, 1(NWP), 2(TC) !')
      end if
    end select

!!!! dont interpolate to time step
    do i = 1, nCells
       forcing%u10(i) = atm_forcing(1)%u10(i)
       forcing%v10(i) = atm_forcing(1)%v10(i)
    end do
    
    call transfer_uv_to_spd_dir_openacc(forcing%u10, forcing%v10, forcing%u10Spd, forcing%u10Dir)

  end subroutine update_atm_forcing_openacc

  subroutine update_forcing_NWP_openacc(time_step)
    integer,         intent(in)  :: time_step
    integer :: i
    allocate(mpi_real_buf_nCells(mpi_total_nCells))
    allocate(mpi_real_buf_nCells_2(mpi_total_nCells))

    if (time_step .gt. 1 .and. time_step .lt. forcing%atm_time_size) then
       do i=1,nCells
         atm_forcing(1)%u10(i) = atm_forcing(2)%u10(i)
         atm_forcing(1)%v10(i) = atm_forcing(2)%v10(i)
       end do
        
       if(mpi_procs>1) then
         if(mpi_rank==0) then
           call read_met_forcing_var(mpi_total_nCells, time_step+1, mpi_real_buf_nCells, mpi_real_buf_nCells_2)
           call mpi_input_exchange_root(mpi_real_buf_nCells,atm_forcing(2)%u10,nCells,'oncell')
           call mpi_input_exchange_root(mpi_real_buf_nCells_2,atm_forcing(2)%v10,nCells,'oncell')
         else
           call mpi_input_exchange_nonroot(mpi_real_buf_nCells,atm_forcing(2)%u10,nCells,'oncell')
           call mpi_input_exchange_nonroot(mpi_real_buf_nCells_2,atm_forcing(2)%v10,nCells,'oncell')
         end if
       else
         call read_met_forcing_var(nCells, time_step+1,atm_forcing(2)%u10,atm_forcing(2)%v10)
       end if
    else if (time_step .eq. 1) then
     
     if(mpi_procs>1) then
       if(mpi_rank==0) then
          call read_met_forcing_var(mpi_total_nCells, time_step, mpi_real_buf_nCells, mpi_real_buf_nCells_2)
          call mpi_input_exchange_root(mpi_real_buf_nCells,atm_forcing(1)%u10,nCells,'oncell')
          call mpi_input_exchange_root(mpi_real_buf_nCells_2,atm_forcing(1)%v10,nCells,'oncell')
       else
          call mpi_input_exchange_nonroot(mpi_real_buf_nCells,atm_forcing(1)%u10,nCells,'oncell')
          call mpi_input_exchange_nonroot(mpi_real_buf_nCells_2,atm_forcing(1)%v10,nCells,'oncell')
       end if
     else
       call read_met_forcing_var(nCells, time_step,atm_forcing(1)%u10,atm_forcing(1)%v10)
     end if

     if(mpi_procs>1) then
       if(mpi_rank==0) then
          call read_met_forcing_var(mpi_total_nCells, time_step+1, mpi_real_buf_nCells, mpi_real_buf_nCells_2)
          call mpi_input_exchange_root(mpi_real_buf_nCells,atm_forcing(2)%u10,nCells,'oncell')
          call mpi_input_exchange_root(mpi_real_buf_nCells_2,atm_forcing(2)%v10,nCells,'oncell')
       else
          call mpi_input_exchange_nonroot(mpi_real_buf_nCells,atm_forcing(2)%u10,nCells,'oncell')
          call mpi_input_exchange_nonroot(mpi_real_buf_nCells_2,atm_forcing(2)%v10,nCells,'oncell')
       end if
     else
       call read_met_forcing_var(nCells,time_step+1,atm_forcing(2)%u10,atm_forcing(2)%v10)
     end if

    else
       do i=1,nCells
         atm_forcing(1)%u10(i) = atm_forcing(2)%u10(i)
         atm_forcing(1)%v10(i) = atm_forcing(2)%v10(i)
       end do
    end if
    deallocate(mpi_real_buf_nCells)
    deallocate(mpi_real_buf_nCells_2)
  end subroutine update_forcing_NWP_openacc

  subroutine update_forcing_TC(time_step, u10, v10)
    integer, intent(in)  :: time_step
    real(real_kind), intent(out):: u10 (:) ! m/s
    real(real_kind), intent(out):: v10 (:) ! m/s

    call calc_distance_from_TC(TC_forcing%lonTC(time_step), TC_forcing%latTC(time_step), TC_forcing%distCell)

    !call calc_TC_pres_field(TC_forcing%rmwTC(time_step),   &
    !                        TC_forcing%presTC(time_step),  &
    !                        TC_forcing%distCell,           &
    !                        pres)
    
    call calc_TC_wind_field(TC_forcing%lonTC(time_step),   &
                            TC_forcing%latTC(time_step),   &
                            TC_forcing%rmwTC(time_step),   &
                            TC_forcing%presTC(time_step),  &
                            TC_forcing%utTC(time_step),    &
                            TC_forcing%vtTC(time_step),    &
                            TC_forcing%distCell,           &
                            u10, v10)


  end subroutine update_forcing_TC

  subroutine update_forcing_TC_openacc(time_step, u10, v10)
    integer, intent(in)  :: time_step
    real(real_kind), intent(out):: u10 (:) ! m/s
    real(real_kind), intent(out):: v10 (:) ! m/s

    call calc_distance_from_TC_openacc(TC_forcing%lonTC(time_step), TC_forcing%latTC(time_step), TC_forcing%distCell)

    !call calc_TC_pres_field(TC_forcing%rmwTC(time_step),   &
    !                        TC_forcing%presTC(time_step),  &
    !                        TC_forcing%distCell,           &
    !                        pres)
    
    call calc_TC_wind_field_openacc(TC_forcing%lonTC(time_step),   &
                            TC_forcing%latTC(time_step),   &
                            TC_forcing%rmwTC(time_step),   &
                            TC_forcing%presTC(time_step),  &
                            TC_forcing%utTC(time_step),    &
                            TC_forcing%vtTC(time_step),    &
                            TC_forcing%distCell,           &
                            u10, v10)

  end subroutine update_forcing_TC_openacc

  subroutine calc_TC_wind_field(lonTC, latTC, rmwTC, presTC, utTC, vtTC, distTC, uc, vc)
      implicit none 
      real(real_kind), intent(in) :: lonTC      ! radian
      real(real_kind), intent(in) :: latTC   
      real(real_kind), intent(in) :: rmwTC      ! m
      real(real_kind), intent(in) :: presTC     ! Pa
      real(real_kind), intent(in) :: utTC       ! m/s
      real(real_kind), intent(in) :: vtTC       ! m/s
      real(real_kind), intent(in) :: distTC(:)  ! m
      real(real_kind), intent(out):: uc(:)  ! Combination wind field, ugEdge + utEdge, m/s
      real(real_kind), intent(out):: vc(:)

      integer iCell
      real(real_kind) :: deltaP
      real(real_kind) :: B         ! Holland B parameter
      real(real_kind) :: rmwTC2
      real(real_kind) :: rmwTC8    ! distTC > 8*rmwTC, give up, use NWP model
      real(real_kind) :: coeff

      real(real_kind) :: X_rad  (lbound(uc,1) : ubound(uc,1)) ! Cartesian coordinate, origin point is TC center
      real(real_kind) :: Y_rad  (lbound(uc,1) : ubound(uc,1))
      real(real_kind) :: OP_rad (lbound(uc,1) : ubound(uc,1)) ! distance betwwn a given Point and TC center, radian 
      
      real(real_kind) :: VgCell (lbound(uc,1) : ubound(uc,1)) ! Gradient wind field
      real(real_kind) :: theta  (lbound(uc,1) : ubound(uc,1))
 

      deltaP = Pe - presTC
      rmwTC2 = 2.0_rk*rmwTC
      do iCell = lbound(lonCell, 1), ubound(lonCell, 1)
         X_rad(iCell)  = lonCell(iCell) - lonTC
         Y_rad(iCell)  = latCell(iCell) - latTC
         OP_rad(iCell) = sqrt(X_rad(iCell)**2 + Y_rad(iCell)**2)
      end do

      do iCell = lbound(theta,1), ubound(theta,1)
        theta(iCell) = asin(Y_rad(iCell)/OP_rad(iCell))
        if(X_rad(iCell) .lt. 0.0_rk) theta(iCell) = pi  - theta(iCell)
        if(theta(iCell) .lt. 0.0_rk) theta(iCell) = pi2 + theta(iCell)
        theta(iCell) = theta(iCell) + inflowA
      end do
      
      select case (trim(adjustL(tc_model_option)))
      case('Holland')
          B = 1.5_rk + ((98000.0_rk - presTC)/100.0_rk)/120.0_rk
         do iCell = lbound(VgCell,1), ubound(VgCell,1) 
            VgCell(iCell) = sqrt(  deltaP*B/rho_air*((rmwTC/distTC(iCell))**B)*exp(-1.0*(rmwTC/distTC(iCell))**B) &
                                  + (distTC(iCell)*fCell(iCell)*0.5)**2  )                                        &
                                  - distTC(iCell)*ABS(fCell(iCell))*0.5

            if (distTC(iCell) .gt. 800000.0) VgCell(iCell) = 0.0 
         end do
     case default
         if (mpi_rank==0) then
         call log_error('Unknow tc_model_option in calc_TC_wind_field, please choose from (Holland)')
         end if
     end select

     do iCell = lbound(uc,1), ubound(uc,1)
          coeff = exp(-0.25_rk*pi*abs(distTC(iCell)-rmwTC)/rmwTC) 
          uc(iCell) = coeff*utTC - 0.8*VgCell(iCell)*sin(theta(iCell))
          vc(iCell) = coeff*vtTC + 0.8*VgCell(iCell)*cos(theta(iCell))
     end do
     
   end subroutine calc_TC_wind_field

  subroutine calc_TC_wind_field_openacc(lonTC, latTC, rmwTC, presTC, utTC, vtTC, distTC, uc, vc)
      implicit none 
      real(real_kind), intent(in) :: lonTC      ! radian
      real(real_kind), intent(in) :: latTC   
      real(real_kind), intent(in) :: rmwTC      ! m
      real(real_kind), intent(in) :: presTC     ! Pa
      real(real_kind), intent(in) :: utTC       ! m/s
      real(real_kind), intent(in) :: vtTC       ! m/s
      real(real_kind), intent(in) :: distTC(:)  ! m
      real(real_kind), intent(out):: uc(:)  ! Combination wind field, ugEdge + utEdge, m/s
      real(real_kind), intent(out):: vc(:)

      integer iCell
      real(real_kind) :: deltaP
      real(real_kind) :: B         ! Holland B parameter
      real(real_kind) :: rmwTC2
      real(real_kind) :: rmwTC8    ! distTC > 8*rmwTC, give up, use NWP model
      real(real_kind) :: coeff

      real(real_kind) :: X_rad  (lbound(uc,1) : ubound(uc,1)) ! Cartesian coordinate, origin point is TC center
      real(real_kind) :: Y_rad  (lbound(uc,1) : ubound(uc,1))
      real(real_kind) :: OP_rad (lbound(uc,1) : ubound(uc,1)) ! distance betwwn a given Point and TC center, radian 
      
      real(real_kind) :: VgCell (lbound(uc,1) : ubound(uc,1)) ! Gradient wind field
      real(real_kind) :: theta  (lbound(uc,1) : ubound(uc,1))
 

      deltaP = Pe - presTC
      rmwTC2 = 2.0_rk*rmwTC
      do iCell = lbound(lonCell, 1), ubound(lonCell, 1)
         X_rad(iCell)  = lonCell(iCell) - lonTC
         Y_rad(iCell)  = latCell(iCell) - latTC
         OP_rad(iCell) = sqrt(X_rad(iCell)**2 + Y_rad(iCell)**2)
      end do

      do iCell = lbound(theta,1), ubound(theta,1)
        theta(iCell) = asin(Y_rad(iCell)/OP_rad(iCell))
        if(X_rad(iCell) .lt. 0.0_rk) theta(iCell) = pi  - theta(iCell)
        if(theta(iCell) .lt. 0.0_rk) theta(iCell) = pi2 + theta(iCell)
        theta(iCell) = theta(iCell) + inflowA
      end do
      
      select case (trim(adjustL(tc_model_option)))
      case('Holland')
          B = 1.5_rk + ((98000.0_rk - presTC)/100.0_rk)/120.0_rk
         do iCell = lbound(VgCell,1), ubound(VgCell,1) 
            VgCell(iCell) = sqrt(  deltaP*B/rho_air*((rmwTC/distTC(iCell))**B)*exp(-1.0*(rmwTC/distTC(iCell))**B) &
                                  + (distTC(iCell)*fCell(iCell)*0.5)**2  )                                        &
                                  - distTC(iCell)*ABS(fCell(iCell))*0.5

            if (distTC(iCell) .gt. 800000.0) VgCell(iCell) = 0.0 
         end do
     case default
         if (mpi_rank==0) then
         call log_error('Unknow tc_model_option in calc_TC_wind_field, please choose from (Holland)')
         end if
     end select

     do iCell = lbound(uc,1), ubound(uc,1)
          coeff = exp(-0.25_rk*pi*abs(distTC(iCell)-rmwTC)/rmwTC) 
          uc(iCell) = coeff*utTC - 0.8*VgCell(iCell)*sin(theta(iCell))
          vc(iCell) = coeff*vtTC + 0.8*VgCell(iCell)*cos(theta(iCell))
     end do
     
   end subroutine calc_TC_wind_field_openacc

   subroutine calc_TC_pres_field(rmwTC, presTC, distTC, pres)
      implicit none 
      real(real_kind), intent(in) :: rmwTC      ! m
      real(real_kind), intent(in) :: presTC     ! Pa
      real(real_kind), intent(in) :: distTC(:)  ! m
      real(real_kind), intent(out):: pres(:)
      
      integer i      

      real(real_kind) :: deltaP
      real(real_kind) :: B         ! Holland B parameter
      real(real_kind) :: rmwTC2
      real(real_kind) :: rmwTC8    ! distTC > 8*rmwTC, give up, use NWP model
 
      deltaP = Pe - presTC
      rmwTC2 = 2.0_rk*rmwTC     

      select case (trim(adjustL(tc_model_option)))
      case('Holland')
          B = 1.5_rk + ((98000.0_rk - presTC)/100.0_rk)/120.0_rk
          do i = lbound(pres,1), ubound(pres,1)
             pres(i) = presTC + deltaP*exp(-1.0_rk*(rmwTC/distTC(i))**B)
          end do
      case('Fujita_Takahashi')
          do i = lbound(pres,1), ubound(pres,1)
             if (distTC(i) .gt. rmwTC2) then ! Takahashi
               pres(i) = Pe - deltaP/sqrt(1.0_rk + 2.0_rk*(distTC(i)/rmwTC)**2)
             else ! Fujita
               pres(i) = Pe - deltaP/sqrt(1.0_rk + 2.0_rk*(distTC(i)/rmwTC)**2)
             end if
          end do
      case('Fujita')
          do i = lbound(pres,1), ubound(pres,1)
             pres(i) = Pe - deltaP/sqrt(1.0_rk + 2.0_rk*(distTC(i)/rmwTC)**2)
          end do
      case('Takahashi')
          do i = lbound(pres,1), ubound(pres,1)
             pres(i) = Pe - deltaP/(1.0_rk + distTC(i)/rmwTC) 
          end do
      case default
          call log_error('Unknow tc_model_option in calc_TC_pres_field, please choose from (Holland, Fujita_Takahashi, Fujita, Takahashi)')
      end select
      
   end subroutine calc_TC_pres_field

   subroutine calc_distance_from_TC(lonTC, latTC, distTC)
      implicit none
      real(real_kind), intent(in) :: lonTC
      real(real_kind), intent(in) :: latTC
      real(real_kind), intent(out):: distTC(:)

      integer i
      do i = lbound(distTC,1), ubound(distTC,1)
        distTC(i) = radius * acos(min(1.0_rk, max(-1.0_rk, sin(latTC) * sin(latCell(i)) + cos(latTC) * cos(latCell(i)) * cos(lonTC - lonCell(i)))))
        distTC(i) = max(distTC(i), 10.0_rk)
      end do
   end subroutine calc_distance_from_TC     

   subroutine calc_distance_from_TC_openacc(lonTC, latTC, distTC)
      implicit none
      real(real_kind), intent(in) :: lonTC
      real(real_kind), intent(in) :: latTC
      real(real_kind), intent(out):: distTC(:)

      integer i
      do i = lbound(distTC,1), ubound(distTC,1)
        distTC(i) = radius * acos(min(1.0_rk, max(-1.0_rk, sin(latTC) * sin(latCell(i)) + cos(latTC) * cos(latCell(i)) * cos(lonTC - lonCell(i)))))
        distTC(i) = max(distTC(i), 10.0_rk)
      end do
   end subroutine calc_distance_from_TC_openacc     

  subroutine allocate_atm_forcing(atm_forcing)
    type(atm_forcing_type), intent(inout) :: atm_forcing
    if (.not. allocated(atm_forcing%u10    ))  allocate(atm_forcing%u10(nCells))
    if (.not. allocated(atm_forcing%v10    ))  allocate(atm_forcing%v10(nCells))
  end subroutine allocate_atm_forcing
   

subroutine transfer_uv_to_spd_dir(u10, v10, u10Spd, u10Dir)
    !!! Copy from WAM model(wam_wind_module.f90)
    real(real_kind), intent(in) :: u10(:)    ! u wind, m/s
    real(real_kind), intent(in) :: v10(:)    ! v wind, m/s
    real(real_kind), intent(out):: u10Spd(:) ! speed, m/s 
    real(real_kind), intent(out):: u10Dir(:) ! direction, north=0, clockwise, radian
    
    real(real_kind) :: wspd        ! total wind speed, m/s
    integer i
    do i = lbound(u10Spd,1), ubound(u10Spd,1)
       u10Spd(i) = sqrt(u10(i)**2 + v10(i)**2)
       if ( u10Spd(i) .ne. 0.0) then
           u10Dir(i) = atan2(u10(i), v10(i)) 
       else
           u10Dir(i) = 0.0
       end if
       if (u10Dir(i) .lt. 0.0) u10Dir(i) = u10Dir(i) + Pi2
    end do
end subroutine transfer_uv_to_spd_dir

subroutine transfer_uv_to_spd_dir_openacc(u10, v10, u10Spd, u10Dir)
    !!! Copy from WAM model(wam_wind_module.f90)
    real(real_kind), intent(in) :: u10(:)    ! u wind, m/s
    real(real_kind), intent(in) :: v10(:)    ! v wind, m/s
    real(real_kind), intent(out):: u10Spd(:) ! speed, m/s 
    real(real_kind), intent(out):: u10Dir(:) ! direction, north=0, clockwise, radian
    
    real(real_kind) :: wspd        ! total wind speed, m/s
    integer i

    do i = lbound(u10Spd,1), ubound(u10Spd,1)
       u10Spd(i) = sqrt(u10(i)**2 + v10(i)**2)
       if ( u10Spd(i) .ne. 0.0) then
           u10Dir(i) = atan2(u10(i), v10(i)) 
       else
           u10Dir(i) = 0.0
       end if
       if (u10Dir(i) .lt. 0.0) u10Dir(i) = u10Dir(i) + Pi2
    end do

end subroutine transfer_uv_to_spd_dir_openacc

end module forcing_mod
