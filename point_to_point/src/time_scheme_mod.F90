module time_scheme_mod

  use params_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use log_mod
  use state_mod
  use wave_propag_mod
  use wam_source_module
  use forcing_mod

#ifdef MPI_DEBUG
  use mod_mpi_test
#endif

  implicit none

  public time_integrate_advection
  public time_integrate_source

contains

  subroutine time_integrate_advection()

    integer, parameter :: one   = -1

    integer iCell, k, m
    if(mpi_procs>1) then
!      call mpi_data_exchange_oncell(maskCell(1:mpi_num_ncells_halo))
      call mpi_data_exchange_oncell(wave(old)%N(1:mpi_num_ncells_halo,1:nDir,1:nFre))
    end if
#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("_wave_old_before_tend.nc"),wave(old)%N(1:nCells,1:nDir,1:nFre),nDir,nFre)
    end if
#endif MPI_DEBUG
    call calc_tend_geog(wave(old)%N(0:,:,:), tend_geog)
#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("_tend_geog.nc"),tend_geog(1:nCells,1:nDir,1:nFre),nDir,nFre)
    end if
#endif MPI_DEBUG
    call update_state(dt, tend_geog, wave(old)%N(0:,:,:), wave(one)%N(0:,:,:))
#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("_wave_one.nc"),wave(one)%N(1:nCells,1:nDir,1:nFre),nDir,nFre)
    end if
#endif MPI_DEBUG
    
    call calc_tend_theta(wave(old)%N(0:,:,:), tend_theta)
#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("_tend_theta.nc"),tend_theta(1:nCells,1:nDir,1:nFre),nDir,nFre)
    end if
#endif MPI_DEBUG
    call update_state(dt, tend_theta,wave(one)%N(0:,:,:), wave(new)%N(0:,:,:))
    
#ifdef MPI_DEBUG
    if(mod(myIter,interval)==0) then
    call mpi_test_output(TRIM("wave_new_N_after_implsch.nc"),wave(new)%N(1:nCells,1:nDir,1:nFre),nDir,nFre)
    call mpi_test_output(TRIM("TAUW_after_implsch.nc"),TAUW(1:nCells))
    call mpi_test_output(TRIM("Z0_after_implsch.nc"),Z0(1:nCells))
    call mpi_test_output(TRIM("USTAR_after_implsch.nc"),USTAR(1:nCells))
    end if
#endif MPI_DEBUG
  end subroutine time_integrate_advection

subroutine time_integrate_source(N)

  real(real_kind), intent(inout)  :: N(1:,:,:)   ! wave actioon
  call IMPLSCH(N(1:nCells,:,:), dt_src, forcing%u10Spd, forcing%u10Dir, TAUW, USTAR, Z0, ROAIRN, WSTAR, bottomDepth(1:nCells))

end subroutine time_integrate_source


end module

