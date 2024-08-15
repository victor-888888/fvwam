MODULE mod_mpi_test
   use mod_mpi_test_variables
   USE mod_mpi_variables
   use mod_mpi_interfaces
   use mod_mpi_io_netcdf
!   use ifport
   IMPLICIT NONE
   INTERFACE mpi_test_output
     MODULE PROCEDURE mpi_data_1d_output_real
     MODULE PROCEDURE mpi_data_2d_output_real
     MODULE PROCEDURE mpi_data_3d_output_real
     MODULE PROCEDURE mpi_data_2d_output_int
   END INTERFACE
  

CONTAINS
   ! time test                                                                      
   SUBROUTINE mpi_test_time_start                                                   
! test in computing processes                                                       
   IF (mpi_rank == 0) THEN                                                          
      CALL DATE_AND_TIME(TIME = mpi_time)                                           
      WRITE (mpi_test_filename, '(I4.4,A1,A,A)'), mpi_procs, '_', TRIM(mpi_time), "_comp_test.txt"
      OPEN (UNIT = MPI_TIME_UNIT, FILE = TRIM(mpi_test_filename), STATUS ='REPLACE', ACTION ='WRITE')
      WRITE (MPI_TIME_UNIT, *) "Number of whole parallel processors is:", mpi_procs
      mpi_clock_io=0
      CALL SYSTEM_CLOCK(mpi_clock_start,mpi_clock_rate,mpi_clock_count_max)         
      CALL SYSTEM_CLOCK(mpi_clock_start)                                            
      CALL SYSTEM_CLOCK(mpi_clock_init_start)                                       
   END IF                                                                           
! test in io processes                                                              
!   IF (mpi_rank == mpi_comp_procs) THEN                                            
!      CALL DATE_AND_TIME(TIME = mpi_time)                                          
!      WRITE (mpi_test_filename, '(I3.3,A1,A,A)'), mpi_procs, '_', TRIM(mpi_time), "_io_test.txt"
!      OPEN (UNIT = MPI_TEST_IO_UNIT, FILE = TRIM(mpi_test_filename), STATUS ='REPLACE', ACTION ='WRITE')
!      WRITE (MPI_TEST_IO_UNIT, *) "Number of whole parallel processors is:", mpi_procs
!      WRITE (MPI_TEST_IO_UNIT, *) "Number of whole IO processors is:", mpi_io_procs
!      mpi_time_real_io_write = 0                                                   
!   END IF                                                                          
   END SUBROUTINE 

   SUBROUTINE mpi_test_compute_time_end                                             
      IF (mpi_rank == 0) THEN                                                       
        CALL SYSTEM_CLOCK(mpi_clock_compute_end)                                    
      END IF                                                                        
   END SUBROUTINE

   SUBROUTINE mpi_test_compute_time_start                                        
      IF (mpi_rank == 0) THEN                                                       
        CALL SYSTEM_CLOCK(mpi_clock_compute_start)                                  
      END IF                                                                        
   END SUBROUTINE 

   SUBROUTINE mpi_test_io_time_start                                        
      IF (mpi_rank == 0) THEN                                                       
        CALL SYSTEM_CLOCK(mpi_clock_io_start)                                  
      END IF                                                                        
   END SUBROUTINE 

   SUBROUTINE mpi_test_io_time_end                                        
      IF (mpi_rank == 0) THEN                                                       
        CALL SYSTEM_CLOCK(mpi_clock_io_end)
        mpi_clock_io=mpi_clock_io+mpi_clock_io_end-mpi_clock_io_start
      END IF                                                                        
   END SUBROUTINE 

    ! test initializing time                                                         
   SUBROUTINE mpi_test_init_time_end                                                
      IF (mpi_rank == 0) THEN                                                       
        CALL SYSTEM_CLOCK(mpi_clock_init_end)                                       
      END IF                                                                        
   END SUBROUTINE 
 
   SUBROUTINE mpi_test_time_end                                                     
   IF (mpi_rank == 0) THEN                                                          
     CALL SYSTEM_CLOCK(mpi_clock_end)                                               
      WRITE (MPI_TIME_UNIT, *) "the whole execution time (seconds) for proc 0 is:"  
      WRITE (MPI_TIME_UNIT, *) dble(mpi_clock_end - mpi_clock_start)/ dble(mpi_clock_rate)
      WRITE (MPI_TIME_UNIT, *) "the whole computing time (seconds) for proc 0 is:"  
      WRITE (MPI_TIME_UNIT, *) dble(mpi_clock_compute_end - mpi_clock_compute_start-mpi_clock_io)/ dble(mpi_clock_rate)
      WRITE (MPI_TIME_UNIT, *) "the whole initializing time (seconds) for proc 0 is:"
      WRITE (MPI_TIME_UNIT, *) dble(mpi_clock_init_end - mpi_clock_init_start)/ dble(mpi_clock_rate)
      WRITE (MPI_TIME_UNIT, *) "the whole IO time (seconds) for proc 0 is:"
      WRITE (MPI_TIME_UNIT, *) dble(mpi_clock_io)/dble(mpi_clock_rate)                                   
      CLOSE (MPI_TIME_UNIT)                                                         
   END IF                                                                           
!   IF (mpi_rank == mpi_comp_procs) THEN                                            
!      WRITE (MPI_TEST_IO_UNIT, *) "the real whole IO - writing time (seconds) for proc ", mpi_comp_procs, " is:", mpi_time_real_io_write
!      CLOSE (MPI_TEST_IO_UNIT)                                                  
!   END IF                                                                       
   END SUBROUTINE 
 
  ! for gdb debuging
!  SUBROUTINE mpi_gdb_test
!    integer :: i,hostname_len,mpi_err,pid
!    character(MPI_MAX_PROCESSOR_NAME) :: hostname
!    call MPI_Get_processor_name(hostname, hostname_len, mpi_err)
!!    write(*,*) "PID= ", pid,", my_rank=",mpi_rank,", hostname= ",  trim(hostname), ", ready for attach"
!    if(mpi_rank==0) then
!      pid=getpid()
!      write(*,*) "my_rank=",mpi_rank,", pid=",", hostname= ",  trim(hostname), ", ready for attach"
!      i=1
!      do 
!        call sleep(1)
!        if (i==0) exit
!      end do
!    end if
!    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
!  END SUBROUTINE

   ! initializing operation for mpi debug and test
   SUBROUTINE mpi_test_debug_init
     IF (mpi_rank == 0) THEN
        OPEN (UNIT = MPI_LOG_UNIT, FILE = TRIM("mpi_test_log.txt"), STATUS ='REPLACE', ACTION ='WRITE')
     END IF
!     ALLOCATE (mpi_variable_3d_nk_2(nlpbz, nk, 2))
!     ALLOCATE (mpi_variable_2d_nk(nlpbz, nk))
!     ALLOCATE (mpi_variable_1d(nlpbz))
!     ALLOCATE (mpi_dp_1d(nlpb))
!     ALLOCATE (mpi_dp_2d_nk(nlpb, nk))
   END SUBROUTINE
   
   ! finalizing operation for mpi debug and test
   SUBROUTINE mpi_test_debug_end
     IF (mpi_rank == 0) THEN
       close(MPI_LOG_UNIT)
     END IF
   END SUBROUTINE

   ! output integer data of 1 dimension
   SUBROUTINE mpi_gather_output_integer(mpi_input, mpi_filename)
      INTEGER, INTENT(IN) :: mpi_input
      CHARACTER(lc), INTENT(IN) :: mpi_filename
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_int_array
      INTEGER :: i
      IF (mpi_rank == 0) THEN
         ALLOCATE (mpi_int_array(mpi_procs))
      END IF

      CALL MPI_GATHER(mpi_input, 1, MPI_INTEGER, mpi_int_array, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      IF (mpi_rank == 0) THEN
         OPEN (UNIT = MPI_TEST_UNIT, FILE = TRIM(mpi_filename), STATUS ='REPLACE', ACTION ='WRITE')
         DO i = 1, mpi_procs
            WRITE (MPI_TEST_UNIT, *) mpi_int_array(i)
         END DO
         CLOSE (MPI_TEST_UNIT)
      END IF
   END SUBROUTINE
  

   SUBROUTINE mpi_data_1d_output_real(file_name, send_var)
      character(*), intent(in) :: file_name
      real(real_kind), intent(in), dimension(:) :: send_var
      character(lc) :: file_name_update
      character(2) :: real_type
      real(real_kind),allocatable,dimension(:) :: recv_buf
      real(real_kind), ALLOCATABLE, DIMENSION(:) :: adjust_buf
      integer :: i


      IF (mpi_rank == 0) THEN
        allocate(recv_buf(mpi_total_nCells))
        ALLOCATE (adjust_buf(mpi_total_nCells))
        CALL MPI_GATHERV(send_var, mpi_nCells, mpi_real_kind, recv_buf, mpi_cell_indexes_all_ordered_counts, &
                       mpi_cell_indexes_all_ordered_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
        CALL mpi_1d_data_adjust(recv_buf, adjust_buf)
        WRITE (file_name_update, '(I4.4,A)') myIter, TRIM(file_name)
!         OPEN (UNIT = 98, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
!!     WRITE (98, *) adjust_buf
!         DO i = 1, mpi_total_nCells
!            WRITE (98, *) adjust_buf(i)
!         END DO
!         CLOSE (98)
        if (real_kind==real_kind_4) then
          real_type='sp'
        else
          real_type='dp'
        end if
        CALL netcdf_write(TRIM(file_name_update),"var_name",adjust_buf,real_type,d1 = mpi_total_nCells,ts = 1)
        deallocate(recv_buf)
        deallocate(adjust_buf)
      ELSE
        CALL MPI_GATHERV(send_var, mpi_nCells, mpi_real_kind, recv_buf, mpi_cell_indexes_all_ordered_counts, &
                       mpi_cell_indexes_all_ordered_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
        
      END IF
   END SUBROUTINE

!==============================================================================  
   SUBROUTINE mpi_1d_data_adjust(buf_original, buf_adjust)                
!==============================================================================  
      REAL(real_kind), INTENT(IN), DIMENSION(:) :: buf_original                            
      REAL(real_kind), INTENT(INOUT), DIMENSION(:) :: buf_adjust                           
      INTEGER ::  i
                                                                     
      DO i = 1, mpi_total_nCells                                                          
        buf_adjust(mpi_cell_indexes_all_ordered(i)) = buf_original(i)                                           
      END DO                                                                        
                                                                                    
   END SUBROUTINE 

!==============================================================================  
   SUBROUTINE mpi_2d_data_adjust_real(buf_original, buf_adjust, d1)                
!==============================================================================  
      REAL(real_kind), INTENT(IN), DIMENSION(:) :: buf_original                            
      REAL(real_kind), INTENT(INOUT), DIMENSION(:,:) :: buf_adjust                           
      INTEGER, INTENT(IN) :: d1                           
      INTEGER ::  i,pos_1d,pos_2d,k,j

      pos_2d=0
      DO k=1,mpi_procs
        DO j=1,d1
          DO i = 1, mpi_cell_indexes_all_ordered_counts(k)    
            pos_2d=pos_2d+1
            pos_1d=mpi_cell_indexes_all_ordered_displs(k)+i
            buf_adjust(mpi_cell_indexes_all_ordered(pos_1d),j) = buf_original(pos_2d)                                           
          END DO     
        END DO
      END DO
                                                                                    
   END SUBROUTINE 

!==============================================================================  
   SUBROUTINE mpi_2d_data_adjust_int(buf_original, buf_adjust, d1)                
!==============================================================================  
      INTEGER, INTENT(IN), DIMENSION(:) :: buf_original                            
      INTEGER, INTENT(INOUT), DIMENSION(:,:) :: buf_adjust                           
      INTEGER, INTENT(IN) :: d1                           
      INTEGER ::  i,pos_1d,pos_2d,k,j

      pos_2d=0
      DO k=1,mpi_procs
        DO j=1,d1
          DO i = 1, mpi_cell_indexes_all_ordered_counts(k)    
            pos_2d=pos_2d+1
            pos_1d=mpi_cell_indexes_all_ordered_displs(k)+i
            buf_adjust(mpi_cell_indexes_all_ordered(pos_1d),j) = buf_original(pos_2d)                                           
          END DO     
        END DO
      END DO
                                                                                    
   END SUBROUTINE 

   SUBROUTINE mpi_data_2d_output_int(file_name, send_var, d1)
      character(*), intent(in) :: file_name
      integer, intent(in) :: d1
      integer, intent(in), dimension(:, :) :: send_var
      character(lc) :: file_name_update
      character(2) :: real_type
      integer :: send_buf(mpi_nCells, d1)
      integer, allocatable,dimension(:) :: recv_buf
      integer, ALLOCATABLE, DIMENSION(:, :) :: adjust_buf
      integer :: i, j, send_count
      integer :: counts_all_2d(mpi_procs)
      integer :: displs_all_2d(mpi_procs)

      send_buf(1:mpi_nCells, 1:d1) = send_var(1:mpi_nCells, 1:d1)
      send_count = mpi_nCells * d1

      IF (mpi_rank == 0) THEN
         ALLOCATE (adjust_buf(mpi_total_nCells, d1))
         ALLOCATE (recv_buf(mpi_total_nCells*d1))
         
         DO i=1,mpi_procs
           counts_all_2d(i)=mpi_cell_indexes_all_ordered_counts(i)*d1
           displs_all_2d(i)=mpi_cell_indexes_all_ordered_displs(i)*d1
         END DO
        
         CALL MPI_GATHERV(send_buf, send_count, MPI_INTEGER, recv_buf, counts_all_2d, &
                       displs_all_2d, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         CALL mpi_2d_data_adjust_int(recv_buf, adjust_buf,d1)
         WRITE (file_name_update, '(I4.4,A)') myIter, TRIM(file_name)

!         OPEN (UNIT = 98, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
!!     WRITE (98, *) adjust_buf
!         DO j = 1, d1
!            DO i = 1, mpi_total_nCells
!               WRITE (98, *) adjust_buf(i, j)
!            END DO
!         END DO
!         CLOSE (98)
        CALL netcdf_write(TRIM(file_name_update),"var_name",adjust_buf,'i4',d1 = mpi_total_nCells, d2 = d1, ts = 1)
        DEALLOCATE(adjust_buf)
        DEALLOCATE(recv_buf)
      ELSE
         CALL MPI_GATHERV(send_buf, send_count, MPI_INTEGER, recv_buf, counts_all_2d, &
                       displs_all_2d, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

      END IF
   END SUBROUTINE

   SUBROUTINE mpi_data_2d_output_real(file_name, send_var, d1)
      character(*), intent(in) :: file_name
      integer, intent(in) :: d1
      real(real_kind), intent(in), dimension(:, :) :: send_var
      character(lc) :: file_name_update
      character(2) :: real_type
      real(real_kind) :: send_buf(mpi_nCells, d1)
      real(real_kind), allocatable,dimension(:) :: recv_buf
      real(real_kind), ALLOCATABLE, DIMENSION(:, :) :: adjust_buf
      integer :: i, j, send_count
      integer :: counts_all_2d(mpi_procs)
      integer :: displs_all_2d(mpi_procs)

      send_buf(1:mpi_nCells, 1:d1) = send_var(1:mpi_nCells, 1:d1)
      send_count = mpi_nCells * d1

      IF (mpi_rank == 0) THEN
         ALLOCATE (adjust_buf(mpi_total_nCells, d1))
         ALLOCATE (recv_buf(mpi_total_nCells*d1))
         
         DO i=1,mpi_procs
           counts_all_2d(i)=mpi_cell_indexes_all_ordered_counts(i)*d1
           displs_all_2d(i)=mpi_cell_indexes_all_ordered_displs(i)*d1
         END DO
        
         CALL MPI_GATHERV(send_buf, send_count, mpi_real_kind, recv_buf, counts_all_2d, &
                       displs_all_2d, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
         CALL mpi_2d_data_adjust_real(recv_buf, adjust_buf,d1)
         WRITE (file_name_update, '(I4.4,A)') myIter, TRIM(file_name)

!         OPEN (UNIT = 98, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
!!     WRITE (98, *) adjust_buf
!         DO j = 1, d1
!            DO i = 1, mpi_total_nCells
!               WRITE (98, *) adjust_buf(i, j)
!            END DO
!         END DO
!         CLOSE (98)
        if (real_kind==real_kind_4) then
          real_type='sp'
        else
          real_type='dp'
        end if
        CALL netcdf_write(TRIM(file_name_update),"var_name",adjust_buf,real_type,d1 = mpi_total_nCells, d2 = d1, ts = 1)
        DEALLOCATE(adjust_buf)
        DEALLOCATE(recv_buf)
      ELSE
         CALL MPI_GATHERV(send_buf, send_count, mpi_real_kind, recv_buf, counts_all_2d, &
                       displs_all_2d, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)

      END IF
   END SUBROUTINE

   SUBROUTINE mpi_data_3d_output_real(file_name, send_var, d1,d2)
      character(*), intent(in) :: file_name
      integer, intent(in) :: d1
      integer, intent(in) :: d2
      real(real_kind), intent(in), dimension(:, :,:) :: send_var
      character(lc) :: file_name_update
      character(2) :: real_type
      real(real_kind) :: send_buf(mpi_nCells, d1,d2)
      real(real_kind), allocatable,dimension(:) :: recv_buf
      real(real_kind), ALLOCATABLE, DIMENSION(:, :,:) :: adjust_buf
      integer :: i, j, send_count
      integer :: counts_all_3d(mpi_procs)
      integer :: displs_all_3d(mpi_procs)

      send_buf(1:mpi_nCells, 1:d1,1:d2) = send_var(1:mpi_nCells, 1:d1,1:d2)
      send_count = mpi_nCells * d1*d2

      IF (mpi_rank == 0) THEN
         ALLOCATE (adjust_buf(mpi_total_nCells, d1,d2))
         ALLOCATE (recv_buf(mpi_total_nCells*d1*d2))
         
         DO i=1,mpi_procs
           counts_all_3d(i)=mpi_cell_indexes_all_ordered_counts(i)*d1*d2
           displs_all_3d(i)=mpi_cell_indexes_all_ordered_displs(i)*d1*d2
         END DO
        
         CALL MPI_GATHERV(send_buf, send_count, mpi_real_kind, recv_buf, counts_all_3d, &
                       displs_all_3d, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
         CALL mpi_3d_data_adjust(recv_buf, adjust_buf,d1,d2)
         WRITE (file_name_update, '(I4.4,A)') myIter, TRIM(file_name)

!         OPEN (UNIT = 98, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
!!     WRITE (98, *) adjust_buf
!         DO j = 1, d1
!            DO i = 1, mpi_total_nCells
!               WRITE (98, *) adjust_buf(i, j)
!            END DO
!         END DO
!         CLOSE (98)
        if (real_kind==real_kind_4) then
          real_type='sp'
        else
          real_type='dp'
        end if
        CALL netcdf_write(TRIM(file_name_update),"var_name",adjust_buf,real_type,d1 = mpi_total_nCells, d2 = d1, d3=d2,ts=1)
        DEALLOCATE(adjust_buf)
        DEALLOCATE(recv_buf)
      ELSE
         CALL MPI_GATHERV(send_buf, send_count, mpi_real_kind, recv_buf, counts_all_3d, &
                       displs_all_3d, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)

      END IF
   END SUBROUTINE

!==============================================================================  
   SUBROUTINE mpi_3d_data_adjust(buf_original, buf_adjust, d1,d2)                
!==============================================================================  
      REAL(real_kind), INTENT(IN), DIMENSION(:) :: buf_original                            
      REAL(real_kind), INTENT(INOUT), DIMENSION(:,:,:) :: buf_adjust                           
      INTEGER, INTENT(IN) :: d1                           
      INTEGER, INTENT(IN) :: d2                           
      INTEGER ::  i,pos_1d,pos_2d,k,j,m

      pos_2d=0
      DO k=1,mpi_procs
        DO m=1,d2
          DO j=1,d1
            DO i = 1, mpi_cell_indexes_all_ordered_counts(k)    
              pos_2d=pos_2d+1
              pos_1d=mpi_cell_indexes_all_ordered_displs(k)+i
              buf_adjust(mpi_cell_indexes_all_ordered(pos_1d),j,m) = buf_original(pos_2d)                                           
            END DO     
          END DO
        END DO
      END DO
                                                                                    
   END SUBROUTINE 

   ! write data of 1 dimension to file
   SUBROUTINE mpi_write_file_1d_dp(file_name, var, var_size)
      character(lc), intent(in) :: file_name
      real(real_kind_8), intent(in), dimension(:) ::var
      integer, intent(in) :: var_size
      integer :: i
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size
         WRITE (99, *) var(i)
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_1d_dp

   ! write data of 2 dimension to file
   SUBROUTINE mpi_write_file_2d_dp(file_name, var, var_size1, var_size2)
      character(lc), intent(in) :: file_name
      real(real_kind_8), intent(in), dimension(:, :) ::var
      integer, intent(in) :: var_size1, var_size2
      integer :: i, j
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size2
         DO j = 1, var_size1
            WRITE (99, *) var(j, i)
         END DO
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_2d_dp

   ! write data of 2 dimension to file
   SUBROUTINE mpi_write_file_2d_i2(file_name, var, var_size1, var_size2)
      character(lc), intent(in) :: file_name
      integer(i2), intent(in), dimension(:, :) ::var
      integer, intent(in) :: var_size1, var_size2
      integer :: i, j
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size2
         DO j = 1, var_size1
            WRITE (99, *) var(j, i)
         END DO
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_2d_i2

   ! write data of 3 dimension to file
   SUBROUTINE mpi_write_file_3d_dp(file_name, var, var_size1, var_size2, var_size3)
      character(lc), intent(in) :: file_name
      real(real_kind_8), intent(in), dimension(:, :,:) ::var
      integer, intent(in) :: var_size1, var_size2, var_size3
      integer :: i, j, k
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO k = 1, var_size3
        DO i = 1, var_size2
           DO j = 1, var_size1
              WRITE (99, *) var(j, i, k)
           END DO
        END DO
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_3d_dp

   ! write integer data of 1 dimension to file
   SUBROUTINE mpi_write_file_1d_integer(file_name, var, var_size)
      character(lc), intent(in) :: file_name
      integer, intent(in), dimension(:) ::var
      integer, intent(in) :: var_size
      integer :: i
      OPEN (UNIT = 99, FILE = TRIM(file_name), STATUS ='REPLACE', ACTION ='WRITE')
      DO i = 1, var_size
         WRITE (99, *) var(i)
      END DO
      CLOSE (99)
   END SUBROUTINE mpi_write_file_1d_integer

   
END MODULE
