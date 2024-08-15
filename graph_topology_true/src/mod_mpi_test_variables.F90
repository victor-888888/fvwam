MODULE mod_mpi_test_variables 
  use const_mod
  IMPLICIT NONE
   integer, public, parameter  :: lc = 256
   integer, public, parameter  :: i1 = selected_int_kind(2)    ! integer byte       
   integer, public, parameter  :: i2 = selected_int_kind(4)    ! integer short   
   integer, public, parameter  :: i4 = selected_int_kind(9)    ! integer long       
   integer, public, parameter  :: i8 = selected_int_kind(14)    ! integer larger
   type, public :: nc_attr                                                           
      character(lc)  :: var_units = 'none'                                           
      character(lc)  :: var_longname = 'none'                                        
   end type nc_attr 
   ! for mpi test function
   INTEGER, PUBLIC, PARAMETER :: MPI_TIME_UNIT = 93   ! unit for outputing performance test results for computing processes
   INTEGER, PUBLIC, PARAMETER :: MPI_TEST_UNIT = 94   ! unit for outputing performance test results for computing processes
   INTEGER, PUBLIC, PARAMETER :: MPI_LOG_UNIT = 95   ! unit for outputing log information
   INTEGER :: myIter=0  ! unit for outputing log information
   integer :: interval=1
   ! only for test
   integer,allocatable,dimension(:,:) :: mpi_cellsOnCell
   integer,allocatable,dimension(:,:) :: mpi_edgesOnCell
   real(kind=real_kind),allocatable,dimension(:) :: mpi_real_buf_nCells_debug
   REAL(kind=real_kind), ALLOCATABLE, DIMENSION(:, :) :: mpi_real_test_2d ! test for size of real in fortran
   REAL(kind=real_kind), ALLOCATABLE, DIMENSION(:, :) :: mpi_real_test_sec_2d ! test for size of real in fortran
   CHARACTER(LEN = 10), PUBLIC :: mpi_time ! use time to name test file
   CHARACTER(LEN = lc), PUBLIC :: mpi_test_filename ! file name of test file
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_start  ! used for counting program runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_end  ! used for counting program runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_rate  ! used for counting program runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_count_max  ! used for counting program runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_compute_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_compute_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_init_start  ! used for counting initializing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_init_end  ! used for counting initializing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_communicate_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_communicate_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_neigh_commu_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_neigh_commu_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_neigh_commu  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_prep_recv  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_prep_recv_start  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), PUBLIC :: mpi_clock_prep_recv_end  ! used for counting computing runtime with SYSTEM_CLOCK function
   INTEGER(KIND = 8), public,allocatable,dimension(:) :: mpi_clock_neigh_commu_all  ! used for counting computing runtime with SYSTEM_CLOCK function
!   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_read_start  ! used for counting computing runtime with SYSTEM_CLOCK function
!   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_read_end  ! used for counting computing runtime with SYSTEM_CLOCK function
!   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_write_start  ! used for counting computing runtime with SYSTEM_CLOCK function
!   INTEGER(KIND = 8), PUBLIC :: mpi_clock_io_write_end  ! used for counting computing runtime with SYSTEM_CLOCK function
!   
!   ! for mpi debug function
!   INTEGER, PUBLIC :: mpi_debug_size  ! test for size of integer and real in mpi
!   INTEGER, PUBLIC :: mpi_int_type ! test for size of integer in fortran
!   REAL(KIND = dp) :: mpi_real_type ! test for size of real in fortran
!   INTEGER, PUBLIC, PARAMETER :: MPI_DEBUG_UNIT = 92   ! unit for outputing to debug file
!   CHARACTER(LEN = lc), PUBLIC :: mpi_debug_filename ! file name of debug file
!   CHARACTER(LEN = lc), PUBLIC :: mpi_var_filename ! file name of debug file
!
!   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :, :) :: mpi_variable_3d_nk_2 ! test for size of real in fortran
!   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_variable_2d_nk ! test for size of real in fortran
!
!   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_variable_1d ! test for size of real in fortran
!   REAL(dp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_dp_1d ! test for size of real in fortran
!   REAL(dp), PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: mpi_dp_2d_nk ! test for size of real in fortran

END MODULE

