module params_mod

  use const_mod
  use mod_mpi_variables

  implicit none

  integer run_days
  integer run_hours
  integer run_minutes
  integer run_seconds
  integer :: start_time(5) = [0, 0, 0, 0, 0]
  integer :: end_time  (5) = [0, 0, 0, 0, 0]
  real(real_kind) :: dt
  integer dt_src_ratio  ! dt_src = dt * dt_src_ratio

  namelist /time_control/ &
    run_days,             &
    run_hours,            &
    run_minutes,          &
    run_seconds,          &
    start_time,           &
    end_time,             &
    dt,                   &
    dt_src_ratio                        

  character(max_name_len)      :: history_interval(1)
  character(max_name_len)      :: time_units         = 'days'
  character(max_file_path_len) :: output_file_prefix = 'N/A'
  character(max_name_len)      :: frames_per_file    = '1 days'
  character(max_file_path_len) :: mesh_file_path
  logical                      :: restart_output_option = .False.     ! output last time step wave%N state to nc file
  logical                      :: restart_input_option  = .False.     ! read wave%N state from nc file at first time step
  logical                      :: station_output_option  = .False.     
  character(max_file_path_len) :: station_info_path      = 'station_info.txt'
  logical                      :: spectrum_output_option = .False.     


  namelist /io_control/          &
    history_interval           , &
    time_units                 , &
    output_file_prefix         , &
    frames_per_file            , &
    mesh_file_path             , &
    restart_output_option      , &
    restart_input_option       , & 
    station_output_option      , &   
    station_info_path          , &
    spectrum_output_option         

  character(max_name_len)      :: case_name = 'wave' 
  character(max_name_len)      :: atm_forcing_interval
  character(max_file_path_len) :: nwp_forcing_file_path
  logical                      :: multiple_nwp_forcing_file = .False.        ! if True, NWP forcing file read by nwp_forcing_file_prefix 
  character(max_file_path_len) :: nwp_forcing_file_prefix   = 'wind_forcing' ! prefix of NWP forcing files
  integer                      :: nwp_forcing_file_num      = 100            ! number of NWP forcing files
  character(max_file_path_len) :: tc_forcing_file_path
  character(max_name_len)      :: tc_model_option         = 'Holland'  ! 
  integer                      :: atm_forcing_option      = 2          ! 1(NWP), 2(TC model)
  logical                      :: interp_wind_steptime    = .True.     ! interpolate wind data to time of integration step
  
  integer                      :: nDir   = 24         ! Size of direction
  integer                      :: nFre   = 25         ! Size of frequency
  real(real_kind)              :: freMin = 0.04_rk    ! Minimum frequency, Hz
  real(real_kind)              :: FETCH  = 30000.0_rk ! fetch length, meter

  real(real_kind)              :: depthMin     = 1.0_rk ! meter
  real(real_kind)              :: depthRatio   = 1.1_rk     
  integer                      :: nDepth       = 69 
  integer                      :: smooth_depth_num = 0 ! 0(no smooth),>=1(smooth number)

  namelist /wave_setting/       &
    case_name                 , &
    atm_forcing_interval      , &
    nwp_forcing_file_path     , &
    multiple_nwp_forcing_file , &
    nwp_forcing_file_prefix   , &
    nwp_forcing_file_num      , &
    tc_forcing_file_path      , &
    tc_model_option           , &
    atm_forcing_option        , &
    interp_wind_steptime      , &
    freMin                    , & 
    nDir                      , &
    nFre                      , &
    FETCH                     , &
    depthMin                  , & 
    depthRatio                , & 
    nDepth                    , & 
    smooth_depth_num                   

contains

  subroutine params_parse_namelist(file_path)

    character(*), intent(in) :: file_path

    if(mpi_rank==0) then
    open(10, file=file_path)
    read(10, nml=time_control)
    read(10, nml=io_control)
    read(10, nml=wave_setting)
    close(10)
    end if
    
  end subroutine params_parse_namelist

end module params_mod
