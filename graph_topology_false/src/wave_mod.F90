module wave_mod

  use params_mod
  use log_mod
  use mesh_mod
  use mask_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use state_mod
  use history_mod
  use time_scheme_mod
  use forcing_mod
  use initial_mod
  use fre_dir_mod
  use wam_source_module
  use mod_mpi_interfaces
  use restart_mod

#ifdef MPI_DEBUG
  use mod_mpi_test
#endif

#ifdef MPI_GDB
  use mod_mpi_test
#endif

#ifdef OPENACC
  use mod_openacc_interfaces
#endif
  implicit none

  public wave_mod_init
  public wave_mod_run

contains

  subroutine wave_mod_init(namelist_file_path)

    character(*), intent(in) :: namelist_file_path

    call mpi_init_basic_env
#ifdef MPI_GDB
    call mpi_gdb_test
#endif
#ifdef MPI_DEBUG
  call mpi_test_debug_init
#endif
#ifdef MPI_TEST                                                                     
   CALL mpi_test_time_start                                                         
#endif
    call params_parse_namelist(namelist_file_path)
    if(mpi_procs>1) then
      call mpi_send_recv_namelist
    end if
    if(mpi_rank==0) then
      call log_init()
    end if
    call time_init()
    call fiona_init(time_units, start_time_format)
    call mesh_init()
    call mask_init()
    call state_init()
    call fre_dir_init()
    call propag_init()
    call wam_source_init() 
    call forcing_init()
    call history_init()
    call set_initial_condition()
  end subroutine wave_mod_init

  subroutine wave_mod_run()
    implicit none
    real(real_kind) :: max_Hs
    real(real_kind) :: max_wind
    real(real_kind) :: reduce_val

#ifdef OPENACC
    call openacc_init
#endif

    call calc_output(wave(old)%N, SWH, MWP, MWD)
    call history_write()
    if (station_output_option) call station_output_write(wave(old)%N)

    if(mpi_procs>1) then
      reduce_val=maxval(SWH)
      CALL MPI_REDUCE(reduce_val,max_Hs,1,mpi_real_kind,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)
      reduce_val=maxval(forcing%u10Spd)
      CALL MPI_REDUCE(reduce_val,max_wind,1,mpi_real_kind,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)
    else
      max_Hs=maxval(SWH)
      max_wind=maxval(forcing%u10Spd)
    end if
    
    if(mpi_rank==0) then
      call log_add_diag('Hs_max',max_Hs)
      call log_add_diag('wind_max',max_wind)
      call log_print_diag(curr_time_format)
    end if
    
    do while (.not. time_is_finished())
#ifdef MPI_DEBUG
    myIter=myIter+1
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("_bottomDepth_before.nc"),bottomDepth(1:nCells))
    call mpi_test_output(TRIM("_Hs_before.nc"),Hs(1:nCells))
    call mpi_test_output(TRIM("_u10Spd_before.nc"),forcing%u10Spd(1:nCells))
    call mpi_test_output(TRIM("_u10Dir_before.nc"),forcing%u10Dir(1:nCells))
    call mpi_test_output(TRIM("_wave_before.nc"),wave(old)%N(1:nCells,1:nDir,1:nFre),nDir,nFre)
    call mpi_test_output(TRIM("_cgCell_before.nc"),cgCell(1:nCells,1:nDir),nDir)
    end if
#endif

      call time_integrate_advection()

      if (time_is_alerted('atm_forcing_update')) then
          call update_atm_forcing_openacc
      end if
      
      if (time_is_alerted('integrate_source')) then
         call time_integrate_source(wave(new)%N(1:nCells,:,:))

         if (interp_wind_steptime) then 
             forcing%atm_time_steps_from_update = forcing%atm_time_steps_from_update + 1
             call interp_wind_to_steptime()
         end if
      end if
      
      call time_advance()
      
      if (time_is_alerted('history_new_file')) then
          call history_create_new_file()
      end if

      if (time_is_alerted('history_write')) then
         call calc_output(wave(old)%N, SWH, MWP, MWD)

         if(mpi_procs>1) then
           reduce_val=maxval(SWH)
           CALL MPI_REDUCE(reduce_val,max_Hs,1,mpi_real_kind,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)
           reduce_val=maxval(forcing%u10Spd)
           CALL MPI_REDUCE(reduce_val,max_wind,1,mpi_real_kind,MPI_MAX,0,MPI_COMM_WORLD,mpi_err)
         else
           max_Hs=maxval(SWH)
           max_wind=maxval(forcing%u10Spd)
         end if

         if(mpi_rank==0) then
           call log_add_diag('Hs_max', max_Hs)
           call log_add_diag('wind_max',max_wind)
           call log_print_diag(curr_time_format)
         end if
        
       call history_write()
       if (station_output_option) call station_output_write(wave(old)%N)
         
      end if
#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("_bottomDepth_after.nc"),bottomDepth(1:nCells))
    call mpi_test_output(TRIM("_Hs_after.nc"),Hs(1:nCells))
    call mpi_test_output(TRIM("_u10Spd_after.nc"),forcing%u10Spd(1:nCells))
    call mpi_test_output(TRIM("_wave_after.nc"),wave(old)%N(1:nCells,1:nDir,1:nFre),nDir,nFre)
    call mpi_test_output(TRIM("_cgCell_after.nc"),cgCell(1:nCells,1:nDir),nDir)
    end if
#endif

    end do

    if (restart_output_option) then
        call restart_output(wave(old))
    end if
  end subroutine wave_mod_run

end module
