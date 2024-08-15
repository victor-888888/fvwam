program wave_main

  use params_mod
  use log_mod
  use wave_mod
  use mod_mpi_test

  implicit none

  integer(kind=8) :: start_clock_main, end_clock_main,end_clock_init,clock_comp
  integer(kind=8),allocatable,dimension(:) :: clock_comp_all
  real(8) :: run_time_main, count_rate

  call system_clock(count=start_clock_main, count_rate=count_rate)

  call wave_mod_init('./namelist.wave')

  call system_clock(count=end_clock_init)
  mpi_clock_neigh_commu=0

#ifdef MPI_TEST                                                                  
      call mpi_test_init_time_end                                                
      call mpi_test_compute_time_start
#endif 

  call wave_mod_run()

#ifdef OPENACC
  call openacc_finalize()
#endif

#ifdef MPI_DEBUG
  call mpi_test_debug_end
#endif
#ifdef MPI_TEST
call mpi_test_compute_time_end
call mpi_test_time_end
#endif
call system_clock(count=end_clock_main)
clock_comp=end_clock_main-end_clock_init

if(mpi_rank==0) then
  write(*,*) "begin to output"
  allocate(clock_comp_all(mpi_procs))
  call mpi_gather(clock_comp,1,mpi_integer8,clock_comp_all,1,mpi_integer8,0,MPI_COMM_WORLD,mpi_err)
  allocate(mpi_clock_neigh_commu_all(mpi_procs))
  call mpi_gather(mpi_clock_neigh_commu,1,mpi_integer8,mpi_clock_neigh_commu_all,1,mpi_integer8,0,MPI_COMM_WORLD,mpi_err)
  write(*,*) 'It took ',(end_clock_main-start_clock_main)/count_rate,' seconds for running whole program in root process'
  write(*,*) 'It took ',(end_clock_init-start_clock_main)/count_rate,' seconds for initializing in root process'
  write(*,*) 'It took ',mpi_clock_neigh_commu/count_rate,' seconds for 3d N communication in root process'
  write(*,*) 'It took ',maxval(mpi_clock_neigh_commu_all)/count_rate,' seconds for maximum 3d N communication'
  write(*,*) 'It took ',(sum(mpi_clock_neigh_commu_all)/mpi_procs)/count_rate,' seconds for average 3d N communication'
  write(*,*) 'It took ',minval(mpi_clock_neigh_commu_all)/count_rate,' seconds for minimum 3d N communication'
  write(*,*) 'It took ',clock_comp/count_rate,' seconds for computing in root process'
  write(*,*) 'It took ',maxval(clock_comp_all)/count_rate,' seconds for maximum computing'
  write(*,*) 'It took ',(sum(clock_comp_all)/mpi_procs)/count_rate,' seconds for average computing'
  write(*,*) 'It took ',minval(clock_comp_all)/count_rate,' seconds for minimum computing'
else
  call mpi_gather(clock_comp,1,mpi_integer8,clock_comp_all,1,mpi_integer8,0,MPI_COMM_WORLD,mpi_err)
  call mpi_gather(mpi_clock_neigh_commu,1,mpi_integer8,mpi_clock_neigh_commu_all,1,mpi_integer8,0,MPI_COMM_WORLD,mpi_err)
end if

  CALL MPI_FINALIZE(mpi_err)

end program wave_main
