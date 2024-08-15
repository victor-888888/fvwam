module mod_mpi_interfaces
  use mod_mpi_variables
  use params_mod
  use mod_mpi_reallocate
  use fiona_mod
  use mod_mpi_test_variables
#ifdef OPENACC
  use mod_openacc_variable
  use openacc
#endif

  implicit none

   ! exchanging variables in computing processes in blocking communication or non - blocking communication based on cell
   INTERFACE mpi_data_exchange_oncell
      MODULE PROCEDURE mpi_int_1d_block_exchange_oncell
      MODULE PROCEDURE mpi_real_1d_block_exchange_oncell
      MODULE PROCEDURE mpi_real_2d_block_exchange_oncell
      MODULE PROCEDURE mpi_real_3d_block_exchange_oncell
!      MODULE PROCEDURE mpi_real_3d_block_exchange_oncell_onebyone
   END INTERFACE

   ! exchanging variables in computing processes in blocking communication or non - blocking communication based on edge
   INTERFACE mpi_data_exchange_onedge
      MODULE PROCEDURE mpi_int_1d_block_exchange_onedge
   END INTERFACE


  ! distribute input data from process 0
  interface mpi_input_exchange_root
    module procedure mpi_input_exchange_root_1d_integer
    module procedure mpi_input_exchange_root_1d_real
    module procedure mpi_input_exchange_root_2d_integer
    module procedure mpi_input_exchange_root_2d_real

  end interface
  
  ! receive input data from process 0
  interface mpi_input_exchange_nonroot
    module procedure mpi_input_exchange_nonroot_1d_integer
    module procedure mpi_input_exchange_nonroot_1d_real
    module procedure mpi_input_exchange_nonroot_2d_integer
    module procedure mpi_input_exchange_nonroot_2d_real

  end interface

  !read data from netcdf files and scatter to other processes
  !it is different from mpi_input_exchange which index of cell is at the last
  !in this interface, index of cell is at the beginning
  interface mpi_input_oncell
    module procedure mpi_input_oncell_1d_real
    module procedure mpi_input_oncell_2d_real
    module procedure mpi_input_oncell_3d_real
  end interface

  !gather data from all process and output based on cell into netcdf files
  interface mpi_output_oncell
    module procedure mpi_output_oncell_1d_real
    module procedure mpi_output_oncell_1d_integer
    module procedure mpi_output_oncell_2d_real
    module procedure mpi_output_oncell_3d_real
  end interface

  !gather data from all process and output based on edge index into netcdf files
  interface mpi_output_onedge
    module procedure mpi_output_onedge_1d_real
  end interface

  ! sort elements in arrays
  interface mpi_sort
    module procedure mpi_sort_one
    module procedure mpi_sort_two

  end interface
  
  !adjust the buffer to the order based on input station 
  interface mpi_adjust_sta
    module procedure mpi_adjust_sta_1d_real
  end interface
  contains

    ! initialize mpi basic environment
    subroutine mpi_init_basic_env

#ifdef OPENACC                                                           
      !Multiple GPUs                                                             
      INTEGER :: InNodeComm,InNodeRank,dev_rank,ierr                                    
      CHARACTER (len=MPI_MAX_PROCESSOR_NAME) :: hostname                         
      INTEGER :: namelength    
      INTEGER :: GPU_AVAILABLE_NUM
      INTEGER, DIMENSION(:), ALLOCATABLE  :: GPU_ARRAY
      NAMELIST /GPU_INFO1/ GPU_AVAILABLE_NUM
      NAMELIST /GPU_INFO2/ GPU_ARRAY
      INTEGER :: IOS
      INTEGER, PARAMETER  :: IU_GPU_RESOURCE = 484 !! FILE UNIT FOR NAMELIST.NC
      LOGICAL :: FILE_EXISTS
#endif  

      call MPI_INIT(mpi_err)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_procs, mpi_err)
      call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_err)
      
      ! set default real type for mpi_real_kind
      if (real_kind == real_kind_4) then
        mpi_real_kind = MPI_REAL
      else
        if(real_kind == real_kind_8) then
          mpi_real_kind = MPI_DOUBLE_PRECISION
        else
          write(*,*) "Error:there is no match for mpi_real_kind"
          stop
        end if
      end if

#ifdef OPENACC                                                               
!Multiple GPUs                                                                       
!FIND GPU DEVICES IN EACH NODE AND BIND GPUS WITH LOCAL RANKS (IN EACH NODE)         
!THIS IS DONE BY USING MPI_COMM_SPLIT_TYPE WITH ATTRIBUTE MPI_COMM_TYPE_SHARED   
! -----------                                                                        
    ! ngpus per node                                                                 
    ngpus = ACC_GET_NUM_DEVICES(ACC_DEVICE_NVIDIA)                                   
    call MPI_GET_PROCESSOR_NAME(hostname, namelength, mpi_err)                       
    if (ngpus .le. 0) then                                                           
        if (mpi_rank .eq. 0) then                                                    
            write(*,'(a)'), '***GPU*** No NVIDIA GPUs available in node: ',  &    
                             hostname(1:namelength),  '. STOP!'                      
            call MPI_Abort(MPI_COMM_WORLD, 1, mpi_err)                               
        endif                                                                        
    else                                                                             
        call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,   &    
             MPI_INFO_NULL, InNodeComm, mpi_err)                                     
        call MPI_COMM_RANK(InNodeComm, InNodeRank, mpi_err)                          
        INQUIRE(FILE="GPU_RESOURCE_SINGLENODE.NML", EXIST=FILE_EXISTS)
        IF (FILE_EXISTS) THEN
            IF (mpi_rank .eq. 0) then
              WRITE(*,*)'READING GPU_ARRAY...'
              WRITE(*,*)'ONLY SUPPORT SPECIFYING GPUS IN SINGLE-NODE OPERATION'
              WRITE(*,*)'IN NAMELIST AVAILABLE GPU NUM LIST FROM SMALL TO BIG.'
            end if
            OPEN(UNIT=IU_GPU_RESOURCE,FILE="GPU_RESOURCE_SINGLENODE.NML",STATUS='OLD')
            READ(IU_GPU_RESOURCE,NML=GPU_INFO1,IOSTAT=IOS)
            IF (IOS.NE.0) THEN
                WRITE(*, *) 'GPU_RESOURCE NAMELIST READ FAILED'
                WRITE(*, *) 'CANNOT SPECIFY GPU DEVICES FOR EACH RANK'
                CALL MPI_ABORT(MPI_COMM_WORLD,1, IERR)
            END IF
            ALLOCATE(GPU_ARRAY(GPU_AVAILABLE_NUM) )
            READ(IU_GPU_RESOURCE,NML=GPU_INFO2,IOSTAT=IOS)
            call ACC_SET_DEVICE_NUM(GPU_ARRAY(InNodeRank+1), ACC_DEVICE_NVIDIA)
            write(*,"(A18,I4,A13,I2,A8,A20)") "***GPU*** RANK:",mpi_rank, &
                    " USING GPU: ",GPU_ARRAY(InNodeRank+1), " AT: ", HostName(1:namelength)
            DEALLOCATE(GPU_ARRAY)
        ELSE
            dev_rank=mod(InNodeRank,ngpus)
            call ACC_SET_DEVICE_NUM(dev_rank, ACC_DEVICE_NVIDIA)               
            write(*,"(A18,I4,A13,I2,A8,A20)") "***GPU*** RANK: ",mpi_rank," USING GPU: ",dev_rank, " AT: ", hostname(1:namelength)
        END IF

    call MPI_COMM_FREE(InNodeComm, mpi_err) 
    end if  
#endif 

    end subroutine
    
    ! send and receive data in namelist
    subroutine mpi_send_recv_namelist
      if(mpi_rank == 0) then
        call mpi_namelist_pack_send
      else
        call mpi_namelist_pack_recv
      end if
      
    end subroutine


    ! send data in namelist from process 0
    subroutine mpi_namelist_pack_send
      integer,parameter :: mpi_buf_size = 10000 ! 1300
      integer :: mpi_position = 0 
      integer :: mpi_size(1) 
      character, allocatable, dimension(:) :: mpi_buf
      allocate(mpi_buf(mpi_buf_size))
      ! namelist in time_control
      call MPI_PACK(run_days, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(run_hours, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(run_minutes, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(run_seconds, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(start_time, 5, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(end_time, 5, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(dt, 1, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(dt_src_ratio, 1, mpi_integer, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
                                                                            
      ! namelist in io_control                                                   
      call MPI_PACK(history_interval, max_name_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(time_units, max_name_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(output_file_prefix, max_file_path_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(frames_per_file, max_name_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(mesh_file_path, max_file_path_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(restart_output_option, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(restart_input_option, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(station_output_option, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(station_info_path, max_file_path_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(spectrum_output_option, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      
      ! namelist in wave_setting                                                 
      call MPI_PACK(case_name, max_name_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(atm_forcing_interval, max_name_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(nwp_forcing_file_path, max_file_path_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(tc_forcing_file_path, max_file_path_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(tc_model_option, max_name_len, MPI_CHARACTER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(atm_forcing_option, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(interp_wind_steptime, 1, MPI_LOGICAL, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(nDir, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(nFre, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(freMin, 1, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(FETCH, 1, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(depthMin, 1, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(depthRatio, 1, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(nDepth, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(smooth_depth_num, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)

      ! send data in namelist
      mpi_size(1)=mpi_position
      call MPI_BCAST(mpi_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      call MPI_BCAST(mpi_buf, mpi_position, MPI_PACKED, 0, MPI_COMM_WORLD, mpi_err)
      ! mpi_position is the real size of buffer in mpi_root_comp                     
      DEALLOCATE (mpi_buf) 
    end subroutine
  
    ! receive data in namelist from process 0
    subroutine mpi_namelist_pack_recv
      integer :: mpi_buf_size 
      integer :: mpi_size(1) 
      integer :: mpi_position = 0 
      character, allocatable, dimension(:) :: mpi_buf

      ! receive data in namelist
      call MPI_BCAST(mpi_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      mpi_buf_size=mpi_size(1)
      allocate(mpi_buf(mpi_buf_size))
      call MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, 0, MPI_COMM_WORLD, mpi_err)

      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, run_days, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, run_hours, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, run_minutes, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, run_seconds, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, start_time, 5, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, end_time, 5, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dt, 1, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, dt_src_ratio, 1, mpi_integer, MPI_COMM_WORLD, mpi_err)

      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, history_interval, max_name_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, time_units, max_name_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, output_file_prefix, max_file_path_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, frames_per_file, max_name_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, mesh_file_path, max_file_path_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, restart_output_option, 1, MPI_LOGICAL, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, restart_input_option, 1, MPI_LOGICAL, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, station_output_option, 1, MPI_LOGICAL, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, station_info_path, max_file_path_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, spectrum_output_option, 1, MPI_LOGICAL, MPI_COMM_WORLD, mpi_err)

      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, case_name, max_name_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, atm_forcing_interval, max_name_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nwp_forcing_file_path, max_file_path_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, tc_forcing_file_path, max_file_path_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, tc_model_option, max_name_len, MPI_CHARACTER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, atm_forcing_option, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, interp_wind_steptime, 1, MPI_LOGICAL, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nDir, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nFre, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, freMin, 1, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, FETCH, 1, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, depthMin, 1, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, depthRatio, 1, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, nDepth, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, smooth_depth_num, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
    
      deallocate(mpi_buf) 
      
    end subroutine
    
    ! pack information from input nc file and broadcast these information
    subroutine mpi_pack_bcast (vertexDegree, maxEdges, maxEdges2)
      integer, parameter :: mpi_buf_size = 50
      integer :: mpi_position = 0 
      character, allocatable, dimension(:) :: mpi_buf
      integer, intent(INOUT) :: vertexDegree
      integer, intent(INOUT) :: maxEdges
      integer, intent(INOUT) :: maxEdges2

      allocate(mpi_buf(mpi_buf_size))

      if(mpi_rank == 0) then
        ! namelist in time_control
        call MPI_PACK(mpi_total_nCells, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
        call MPI_PACK(mpi_total_nEdges, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
        call MPI_PACK(mpi_total_nVertices, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
        call MPI_PACK(vertexDegree, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
        call MPI_PACK(maxEdges, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
        call MPI_PACK(maxEdges2, 1, MPI_INTEGER, mpi_buf, mpi_buf_size, mpi_position, MPI_COMM_WORLD, mpi_err)
        
        call MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, 0, MPI_COMM_WORLD, mpi_err)
      else
        call MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, 0, MPI_COMM_WORLD, mpi_err)

        call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, mpi_total_nCells, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
        call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, mpi_total_nEdges, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
        call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, mpi_total_nVertices, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
        call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, vertexDegree, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
        call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, maxEdges, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)
        call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_position, maxEdges2, 1, MPI_INTEGER, MPI_COMM_WORLD, mpi_err)

      end if
        
      deallocate(mpi_buf)
      
    end subroutine


    ! read cell partition from process 0, and broadcast to other processes
    subroutine mpi_bcast_process_partition(nCells,maxEdges)
      integer, intent(inout) :: nCells
      integer, intent(in) :: maxEdges
      integer :: mpi_io_stat,i
      integer,allocatable,dimension(:) :: mpi_cell_indexes_tmp
      character(LEN = MPI_MIDIUM_LEN) :: mpi_io_msg
      allocate(mpi_cell_indexes_all(mpi_total_nCells))
      allocate(mpi_cell_indexes_tmp(mpi_total_nCells))

      if(mpi_rank == 0) then
        open (UNIT = MPI_INPUT_UNIT, FILE ='tmppartition', STATUS ='OLD', IOSTAT = mpi_io_stat, ACTION ='READ', IOMSG = mpi_io_msg)
                                                                                    
        if (mpi_io_stat /= 0) then                                                  
           write (*, *) 'opening tmppartion is failed, error message = ', mpi_io_msg
           stop                                                                     
        else                                                                        
           write (*, *) 'open tmppartion file successfully'              
        end if    


        read (UNIT = MPI_INPUT_UNIT, FMT =*, IOSTAT = mpi_io_stat, IOMSG = mpi_io_msg) mpi_cell_indexes_all(1:mpi_total_nCells)
                                                                                
        if (mpi_io_stat /= 0) then                                              
           write (*, *) 'reading tmppartion is failed, error message = ', mpi_io_msg
           stop                                                                  
        end if                                                                  
                                                                                
        close (UNIT = MPI_INPUT_UNIT)       

        if(maxval(mpi_cell_indexes_all(1:mpi_total_nCells))+1 .ne. mpi_procs) then
          write(*,*) "process number is not same between tmppartition file and argument specified by mpirun, please check them"
          stop
        end if
                                                                                
        call MPI_BCAST(mpi_cell_indexes_all, mpi_total_nCells, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
 
      else
        call MPI_BCAST(mpi_cell_indexes_all, mpi_total_nCells, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      end if

      mpi_nCells = 0

      do i = 1,mpi_total_nCells
        if(mpi_cell_indexes_all(i)== mpi_rank) then
          mpi_nCells = mpi_nCells + 1
          mpi_cell_indexes_tmp(mpi_nCells)= i
        end if
      end do

      nCells = mpi_nCells

      allocate(mpi_cell_indexes(nCells))
      mpi_cell_indexes(1:nCells) = mpi_cell_indexes_tmp(1:nCells)

      if (mpi_rank == 0) then
        allocate(mpi_cell_indexes_all_ordered(mpi_total_nCells))
        allocate(mpi_cell_indexes_all_ordered_counts(mpi_procs))
        allocate(mpi_cell_indexes_all_ordered_displs(mpi_procs))
        call MPI_GATHER(nCells, 1, MPI_INTEGER, mpi_cell_indexes_all_ordered_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

!        do i = 1,mpi_procs
!          mpi_indexes_all_ordered_counts_2d_maxedges(i)= mpi_cell_indexes_all_ordered_counts(i)* maxEdges
!          mpi_indexes_all_ordered_counts_2d_3size(i)= mpi_cell_indexes_all_ordered_counts(i)* 3
!          mpi_cell_indexes_all_ordered_counts_3d_3size_maxedges(i)= mpi_cell_indexes_all_ordered_counts(i)* 3 * maxEdges
!          mpi_cell_indexes_all_ordered_counts_3d_3size_2size(i)= mpi_cell_indexes_all_ordered_counts(i)* 3 * 2
!        end do

        mpi_cell_indexes_all_ordered_displs(1) = 0                                        
!        mpi_indexes_all_ordered_displs_2d_maxedges(1)= 0
!        mpi_indexes_all_ordered_displs_2d_3size(1)= 0
!        mpi_cell_indexes_all_ordered_displs_3d_3size_maxedges(1)= 0
!        mpi_cell_indexes_all_ordered_displs_3d_3size_2size(1)= 0

        DO i = 2, mpi_procs                                             
          mpi_cell_indexes_all_ordered_displs(i) = mpi_cell_indexes_all_ordered_displs(i - 1) + &
                                                   mpi_cell_indexes_all_ordered_counts(i - 1) 
!          mpi_indexes_all_ordered_displs_2d_maxedges(i)= mpi_indexes_all_ordered_displs_2d_maxedges(i - 1) + &
!                                                              mpi_indexes_all_ordered_counts_2d_maxedges(i - 1)
!          mpi_indexes_all_ordered_displs_2d_3size(i)= mpi_indexes_all_ordered_displs_2d_3size(i - 1) + &
!                                                           mpi_indexes_all_ordered_counts_2d_3size(i - 1)
!          mpi_cell_indexes_all_ordered_displs_3d_3size_maxedges(i)= mpi_cell_indexes_all_ordered_displs_3d_3size_maxedges(i - 1) + &
!                                                                    mpi_cell_indexes_all_ordered_counts_3d_3size_maxedges(i - 1)
!          mpi_cell_indexes_all_ordered_displs_3d_3size_2size(i)= mpi_cell_indexes_all_ordered_displs_3d_3size_2size(i - 1) + &
!                                                                 mpi_cell_indexes_all_ordered_counts_3d_3size_2size(i - 1)
        END DO 

        CALL MPI_GATHERV(mpi_cell_indexes, nCells, MPI_INTEGER, mpi_cell_indexes_all_ordered, mpi_cell_indexes_all_ordered_counts, &
                         mpi_cell_indexes_all_ordered_displs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

      else
        call MPI_GATHER(nCells, 1, MPI_INTEGER, mpi_cell_indexes_all_ordered_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

        CALL MPI_GATHERV(mpi_cell_indexes, nCells, MPI_INTEGER, mpi_cell_indexes_all_ordered, mpi_cell_indexes_all_ordered_counts, &
                         mpi_cell_indexes_all_ordered_displs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      end if

      deallocate(mpi_cell_indexes_tmp)
      
    end subroutine

    ! distribute input data from process 0
    subroutine mpi_input_exchange_root_1d_integer(input_buf,output_buf,output_size,mode)
      integer, intent(in), dimension(:) :: input_buf
      integer, intent(inout), dimension(:) :: output_buf
      integer, intent(in) :: output_size
      character(len = 6),intent(in) :: mode
      integer, allocatable, dimension(:) :: send_buf
      integer :: i
      select case (mode)
      case ("oncell")
        allocate(send_buf(mpi_total_nCells))
        do i = 1,mpi_total_nCells
          send_buf(i)= input_buf(mpi_cell_indexes_all_ordered(i))
        end do
        call MPI_SCATTERV(send_buf,mpi_cell_indexes_all_ordered_counts,mpi_cell_indexes_all_ordered_displs,MPI_INTEGER, &
                          output_buf,output_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
      case ("onedge")
        allocate(send_buf(mpi_edge_total_ordered_counts))
        do i = 1,mpi_edge_total_ordered_counts
          send_buf(i)= input_buf(mpi_edge_indexes_all_ordered(i))
        end do
        call MPI_SCATTERV(send_buf,mpi_edge_indexes_all_ordered_counts,mpi_edge_indexes_all_ordered_displs,MPI_INTEGER, &
                          output_buf,output_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
      case ("onvert")
        allocate(send_buf(mpi_vertex_total_ordered_counts))
        do i = 1,mpi_vertex_total_ordered_counts
          send_buf(i)= input_buf(mpi_vertex_indexes_all_ordered(i))
        end do
        call MPI_SCATTERV(send_buf,mpi_vertex_indexes_all_ordered_counts,mpi_vertex_indexes_all_ordered_displs,MPI_INTEGER, &
                          output_buf,output_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
      case DEFAULT
        write(*,*) "error: wrong mode, only oncell,onedge,onvert are allowed"
        stop
      end select
    end subroutine

    ! distribute input data from process 0
    subroutine mpi_input_exchange_root_2d_integer(input_buf,output_buf,dimen1,dimen2,mode)
      integer, intent(in), dimension(:,:) :: input_buf
      integer, intent(inout), dimension(:,:) :: output_buf
      integer, intent(in) :: dimen1
      integer, intent(in) :: dimen2
      character(len = 6),intent(in) :: mode
      integer, allocatable, dimension(:) :: send_buf
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_counts_2d
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_displs_2d
      integer :: i,j,buf_postion
      buf_postion = 1
      select case (mode)
      case ("oncell")
        allocate(send_buf(mpi_total_nCells * dimen1))
        allocate(mpi_indexes_all_ordered_counts_2d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_2d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_2d(i)= mpi_cell_indexes_all_ordered_counts(i)* dimen1
        end do
        mpi_indexes_all_ordered_displs_2d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_2d(i)= mpi_indexes_all_ordered_displs_2d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_2d(i - 1)
        end do
        do i = 1,mpi_total_nCells
          send_buf(buf_postion:buf_postion + dimen1 - 1) = input_buf(1:dimen1,mpi_cell_indexes_all_ordered(i))
          buf_postion = buf_postion + dimen1
        end do
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,MPI_INTEGER, &
                          output_buf,dimen1 * dimen2,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
        deallocate(mpi_indexes_all_ordered_counts_2d)
        deallocate(mpi_indexes_all_ordered_displs_2d)
      case ("onedge")
        allocate(send_buf(mpi_edge_total_ordered_counts * dimen1))
        allocate(mpi_indexes_all_ordered_counts_2d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_2d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_2d(i)= mpi_edge_indexes_all_ordered_counts(i)* dimen1
        end do
        mpi_indexes_all_ordered_displs_2d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_2d(i)= mpi_indexes_all_ordered_displs_2d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_2d(i - 1)
        end do
        do i = 1,mpi_edge_total_ordered_counts
          send_buf(buf_postion:buf_postion + dimen1 - 1) = input_buf(1:dimen1,mpi_edge_indexes_all_ordered(i))
          buf_postion = buf_postion + dimen1
        end do
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,MPI_INTEGER, &
                          output_buf,dimen1 * dimen2,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
        deallocate(mpi_indexes_all_ordered_counts_2d)
        deallocate(mpi_indexes_all_ordered_displs_2d)
      case ("onvert")
        allocate(send_buf(mpi_vertex_total_ordered_counts * dimen1))
        allocate(mpi_indexes_all_ordered_counts_2d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_2d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_2d(i)= mpi_vertex_indexes_all_ordered_counts(i)* dimen1
        end do
        mpi_indexes_all_ordered_displs_2d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_2d(i)= mpi_indexes_all_ordered_displs_2d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_2d(i - 1)
        end do
        do i = 1,mpi_vertex_total_ordered_counts
          send_buf(buf_postion:buf_postion + dimen1 - 1) = input_buf(1:dimen1,mpi_vertex_indexes_all_ordered(i))
          buf_postion = buf_postion + dimen1
        end do
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,MPI_INTEGER, &
                          output_buf,dimen1 * dimen2,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
        deallocate(mpi_indexes_all_ordered_counts_2d)
        deallocate(mpi_indexes_all_ordered_displs_2d)
      case DEFAULT
        write(*,*) "error: wrong mode, only oncell,onedge,onvert are allowed"
        stop
      end select
    end subroutine

    ! distribute input data from process 0
    subroutine mpi_input_exchange_root_2d_real(input_buf,output_buf,dimen1,dimen2,mode)
      real(kind = real_kind), intent(in), dimension(:,:) :: input_buf
      real(kind = real_kind), intent(inout), dimension(:,:) :: output_buf
      integer, intent(in) :: dimen1
      integer, intent(in) :: dimen2
      character(len = 6),intent(in) :: mode
      real(kind = real_kind), allocatable, dimension(:) :: send_buf
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_counts_2d
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_displs_2d
      integer :: i,j,buf_postion
      buf_postion = 1
      select case (mode)
      case ("oncell")
        allocate(send_buf(mpi_total_nCells * dimen1))
        allocate(mpi_indexes_all_ordered_counts_2d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_2d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_2d(i)= mpi_cell_indexes_all_ordered_counts(i)* dimen1
        end do
        mpi_indexes_all_ordered_displs_2d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_2d(i)= mpi_indexes_all_ordered_displs_2d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_2d(i - 1)
        end do
        do i = 1,mpi_total_nCells
          send_buf(buf_postion:buf_postion + dimen1 - 1) = input_buf(1:dimen1,mpi_cell_indexes_all_ordered(i))
          buf_postion = buf_postion + dimen1
        end do
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,MPI_INTEGER, &
                          output_buf,dimen1 * dimen2,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
        deallocate(mpi_indexes_all_ordered_counts_2d)
        deallocate(mpi_indexes_all_ordered_displs_2d)
      case ("onedge")
        allocate(send_buf(mpi_edge_total_ordered_counts * dimen1))
        allocate(mpi_indexes_all_ordered_counts_2d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_2d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_2d(i)= mpi_edge_indexes_all_ordered_counts(i)* dimen1
        end do
        mpi_indexes_all_ordered_displs_2d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_2d(i)= mpi_indexes_all_ordered_displs_2d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_2d(i - 1)
        end do
        do i = 1,mpi_edge_total_ordered_counts
          send_buf(buf_postion:buf_postion + dimen1 - 1) = input_buf(1:dimen1,mpi_edge_indexes_all_ordered(i))
          buf_postion = buf_postion + dimen1
        end do
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,mpi_real_kind, &
                          output_buf,dimen1 * dimen2,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
        deallocate(mpi_indexes_all_ordered_counts_2d)
        deallocate(mpi_indexes_all_ordered_displs_2d)
      case ("onvert")
        allocate(send_buf(mpi_vertex_total_ordered_counts * dimen1))
        allocate(mpi_indexes_all_ordered_counts_2d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_2d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_2d(i)= mpi_vertex_indexes_all_ordered_counts(i)* dimen1
        end do
        mpi_indexes_all_ordered_displs_2d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_2d(i)= mpi_indexes_all_ordered_displs_2d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_2d(i - 1)
        end do
        do i = 1,mpi_vertex_total_ordered_counts
          send_buf(buf_postion:buf_postion + dimen1 - 1) = input_buf(1:dimen1,mpi_vertex_indexes_all_ordered(i))
          buf_postion = buf_postion + dimen1
        end do
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,mpi_real_kind, &
                          output_buf,dimen1 * dimen2,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
        deallocate(mpi_indexes_all_ordered_counts_2d)
        deallocate(mpi_indexes_all_ordered_displs_2d)

      case DEFAULT
        write(*,*) "error: wrong mode, only oncell,onedge,onvert are allowed"
        stop
      end select
    end subroutine

    ! distribute input data from process 0
    subroutine mpi_input_exchange_root_1d_real(input_buf,output_buf,output_size,mode)
      real(kind = real_kind), intent(in), dimension(:) :: input_buf
      real(kind = real_kind), intent(inout), dimension(:) :: output_buf
      integer, intent(in) :: output_size
      character(len = 6),intent(in) :: mode
      real(kind = real_kind), allocatable, dimension(:) :: send_buf
      integer :: i
      select case (mode)
      case ("oncell")
        allocate(send_buf(mpi_total_nCells))
        do i = 1,mpi_total_nCells
          send_buf(i)= input_buf(mpi_cell_indexes_all_ordered(i))
        end do
        call MPI_SCATTERV(send_buf,mpi_cell_indexes_all_ordered_counts,mpi_cell_indexes_all_ordered_displs,mpi_real_kind, &
                          output_buf,output_size,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)

        deallocate(send_buf)
      case ("onedge")
        allocate(send_buf(mpi_edge_total_ordered_counts))
        do i = 1,mpi_edge_total_ordered_counts
          send_buf(i)= input_buf(mpi_edge_indexes_all_ordered(i))
        end do
        call MPI_SCATTERV(send_buf,mpi_edge_indexes_all_ordered_counts,mpi_edge_indexes_all_ordered_displs,mpi_real_kind, &
                          output_buf,output_size,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
      case ("onvert")
        allocate(send_buf(mpi_vertex_total_ordered_counts))
        do i = 1,mpi_vertex_total_ordered_counts
          send_buf(i)= input_buf(mpi_vertex_indexes_all_ordered(i))
        end do
        call MPI_SCATTERV(send_buf,mpi_vertex_indexes_all_ordered_counts,mpi_vertex_indexes_all_ordered_displs,mpi_real_kind, &
                          output_buf,output_size,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
      case DEFAULT
        write(*,*) "error: wrong mode, only oncell,onedge,onvert are allowed"
        stop
      end select
    end subroutine

    ! receive input data from process 0
    subroutine mpi_input_exchange_nonroot_1d_integer(input_buf,output_buf,output_size,mode)
      integer, intent(in), dimension(:) :: input_buf
      integer, intent(inout), dimension(:) :: output_buf
      integer, intent(in) :: output_size
      character(len = 6),intent(in) :: mode
      integer, allocatable, dimension(:) :: send_buf
      integer :: i
      select case (mode)
      case ("oncell")
        call MPI_SCATTERV(send_buf,mpi_cell_indexes_all_ordered_counts,mpi_cell_indexes_all_ordered_displs,MPI_INTEGER, &
                          output_buf,output_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      case ("onedge")
        call MPI_SCATTERV(send_buf,mpi_edge_indexes_all_ordered_counts,mpi_edge_indexes_all_ordered_displs,MPI_INTEGER, &
                          output_buf,output_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      case ("onvert")
        call MPI_SCATTERV(send_buf,mpi_vertex_indexes_all_ordered_counts,mpi_vertex_indexes_all_ordered_displs,MPI_INTEGER, &
                          output_buf,output_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      case DEFAULT
        write(*,*) "error: wrong mode, only oncell,onedge,onvert are allowed"
        stop
      end select
    end subroutine

    ! receive input data from process 0
    subroutine mpi_input_exchange_nonroot_1d_real(input_buf,output_buf,output_size,mode)
      real(kind = real_kind), intent(in), dimension(:) :: input_buf
      real(kind = real_kind), intent(inout), dimension(:) :: output_buf
      integer, intent(in) :: output_size
      character(len = 6),intent(in) :: mode
      real(kind = real_kind), allocatable, dimension(:) :: send_buf
      integer :: i
      select case (mode)
      case ("oncell")
        call MPI_SCATTERV(send_buf,mpi_cell_indexes_all_ordered_counts,mpi_cell_indexes_all_ordered_displs,mpi_real_kind, &
                          output_buf,output_size,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
      case ("onedge")
        call MPI_SCATTERV(send_buf,mpi_edge_indexes_all_ordered_counts,mpi_edge_indexes_all_ordered_displs,mpi_real_kind, &
                          output_buf,output_size,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
      case ("onvert")
        call MPI_SCATTERV(send_buf,mpi_vertex_indexes_all_ordered_counts,mpi_vertex_indexes_all_ordered_displs,mpi_real_kind, &
                          output_buf,output_size,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
      case DEFAULT
        write(*,*) "error: wrong mode, only oncell,onedge,onvert are allowed"
        stop
      end select
    end subroutine

    ! distribute input data from process 0
    subroutine mpi_input_exchange_nonroot_2d_integer(input_buf,output_buf,dimen1,dimen2,mode)
      integer, intent(in), dimension(:,:) :: input_buf
      integer, intent(inout), dimension(:,:) :: output_buf
      integer, intent(in) :: dimen1
      integer, intent(in) :: dimen2
      character(len = 6),intent(in) :: mode
      integer, allocatable, dimension(:) :: send_buf
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_counts_2d
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_displs_2d
      integer :: i,j
      select case (mode)
      case ("oncell")
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,MPI_INTEGER, &
                          output_buf,dimen1 * dimen2,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      case ("onedge")
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,MPI_INTEGER, &
                          output_buf,dimen1 * dimen2,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      case ("onvert")
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,MPI_INTEGER, &
                          output_buf,dimen1 * dimen2,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      case DEFAULT
        write(*,*) "error: wrong mode, only oncell,onedge,onvert are allowed"
        stop
      end select
    end subroutine

    ! distribute input data from process 0
    subroutine mpi_input_exchange_nonroot_2d_real(input_buf,output_buf,dimen1,dimen2,mode)
      real(kind = real_kind), intent(in), dimension(:,:) :: input_buf
      real(kind = real_kind), intent(inout), dimension(:,:) :: output_buf
      integer, intent(in) :: dimen1
      integer, intent(in) :: dimen2
      character(len = 6),intent(in) :: mode
      real(kind = real_kind), allocatable, dimension(:) :: send_buf
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_counts_2d
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_displs_2d
      integer :: i,j
      select case (mode)
      case ("oncell")
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,MPI_INTEGER, &
                          output_buf,dimen1 * dimen2,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
      case ("onedge")
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,mpi_real_kind, &
                          output_buf,dimen1 * dimen2,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
      case ("onvert")
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,mpi_real_kind, &
                          output_buf,dimen1 * dimen2,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
      case DEFAULT
        write(*,*) "error: wrong mode, only oncell,onedge,onvert are allowed"
        stop
      end select
    end subroutine

    ! gather 1 dimension of data and output into netcdf files
    subroutine mpi_output_onedge_1d_real(dataset_name,var_name,input)
      character(*), intent(IN) :: dataset_name
      character(*), intent(IN) :: var_name
      real(kind=real_kind),intent(IN),dimension(:) :: input
      real(kind=real_kind),allocatable,dimension(:) :: mpi_buf_recv
      real(kind=real_kind),allocatable,dimension(:) :: mpi_buf_output

      if(mpi_rank==0) then
        !pay attendtion mpi_edge_total_ordered_counts is bigger than mpi_total_nEdges, due to edge could be duplicated
        allocate(mpi_buf_recv(mpi_edge_total_ordered_counts))
        allocate(mpi_buf_output(mpi_total_nEdges))
        call MPI_GATHERV(input,mpi_nEdges,mpi_real_kind,mpi_buf_recv,mpi_edge_indexes_all_ordered_counts,&
                          mpi_edge_indexes_all_ordered_displs, mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        call mpi_onedge_output_real_adjust_1d(mpi_buf_recv,mpi_buf_output)
        call fiona_output(dataset_name, var_name, mpi_buf_output)
        deallocate(mpi_buf_recv)
        deallocate(mpi_buf_output)
      else
        call MPI_GATHERV(input,mpi_nEdges,mpi_real_kind,mpi_buf_recv,mpi_edge_indexes_all_ordered_counts,&
                          mpi_edge_indexes_all_ordered_displs, mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)

      end if

      
    end subroutine

    ! read netcdf files and scatter to other processes
    subroutine mpi_input_oncell_1d_real(input_buf,output_buf)
      real(kind = real_kind), intent(in), dimension(:) :: input_buf
      real(kind = real_kind), intent(inout), dimension(:) :: output_buf
      real(kind = real_kind), allocatable, dimension(:) :: send_buf
      integer :: i
      if(mpi_rank==0) then
        allocate(send_buf(mpi_total_nCells))
        do i = 1,mpi_total_nCells
          send_buf(i)= input_buf(mpi_cell_indexes_all_ordered(i))
        end do
        call MPI_SCATTERV(send_buf,mpi_cell_indexes_all_ordered_counts,mpi_cell_indexes_all_ordered_displs,mpi_real_kind, &
                          output_buf,mpi_nCells,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)

        deallocate(send_buf)
      else

        call MPI_SCATTERV(send_buf,mpi_cell_indexes_all_ordered_counts,mpi_cell_indexes_all_ordered_displs,mpi_real_kind, &
                          output_buf,mpi_nCells,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
      end if

      
    end subroutine

    ! read netcdf files and scatter to other processes
    subroutine mpi_input_oncell_2d_real(input_buf,output_buf)
      real(kind = real_kind), intent(in), dimension(:,:) :: input_buf
      real(kind = real_kind), intent(inout), dimension(:,:) :: output_buf
      real(kind = real_kind), allocatable, dimension(:) :: send_buf
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_counts_2d
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_displs_2d
      integer :: i,j,dimen1
      dimen1=SIZE(output_buf,2)
      if(mpi_rank==0) then
        allocate(send_buf(mpi_total_nCells * dimen1))
        allocate(mpi_indexes_all_ordered_counts_2d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_2d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_2d(i)= mpi_cell_indexes_all_ordered_counts(i)* dimen1
        end do
        mpi_indexes_all_ordered_displs_2d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_2d(i)= mpi_indexes_all_ordered_displs_2d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_2d(i - 1)
        end do
        call mpi_oncell_input_real_adjust_2d(input_buf,send_buf)
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,mpi_real_kind, &
                          output_buf,mpi_nCells*dimen1,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
        deallocate(mpi_indexes_all_ordered_counts_2d)
        deallocate(mpi_indexes_all_ordered_displs_2d)
      else
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_2d,mpi_indexes_all_ordered_displs_2d,mpi_real_kind, &
                          output_buf,mpi_nCells*dimen1,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)

      end if

      
    end subroutine

    ! read netcdf files and scatter to other processes
    subroutine mpi_input_oncell_3d_real(input_buf,output_buf)
      real(kind = real_kind), intent(in), dimension(:,:,:) :: input_buf
      real(kind = real_kind), intent(inout), dimension(:,:,:) :: output_buf
      real(kind = real_kind), allocatable, dimension(:) :: send_buf
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_counts_3d
      integer, allocatable, dimension(:) :: mpi_indexes_all_ordered_displs_3d
      integer :: i,j,dimen1,dimen2
      dimen1=SIZE(output_buf,2)
      dimen2=SIZE(output_buf,3)
      if(mpi_rank==0) then
        allocate(send_buf(mpi_total_nCells * dimen1*dimen2))
        allocate(mpi_indexes_all_ordered_counts_3d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_3d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_3d(i)= mpi_cell_indexes_all_ordered_counts(i)* dimen1*dimen2
        end do
        mpi_indexes_all_ordered_displs_3d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_3d(i)= mpi_indexes_all_ordered_displs_3d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_3d(i - 1)
        end do
        call mpi_oncell_input_real_adjust_3d(input_buf,send_buf)
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_3d,mpi_indexes_all_ordered_displs_3d,mpi_real_kind, &
                          output_buf,mpi_nCells*dimen1*dimen2,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        deallocate(send_buf)
        deallocate(mpi_indexes_all_ordered_counts_3d)
        deallocate(mpi_indexes_all_ordered_displs_3d)
      else
        call MPI_SCATTERV(send_buf,mpi_indexes_all_ordered_counts_3d,mpi_indexes_all_ordered_displs_3d,mpi_real_kind, &
                          output_buf,mpi_nCells*dimen1*dimen2,mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)

      end if

      
    end subroutine

    ! gather 1 dimension of data and output into netcdf files
    subroutine mpi_output_oncell_1d_real(dataset_name,var_name,input)
      character(*), intent(IN) :: dataset_name
      character(*), intent(IN) :: var_name
      real(kind=real_kind),intent(IN),dimension(:) :: input
      real(kind=real_kind),allocatable,dimension(:) :: mpi_buf_recv
      real(kind=real_kind),allocatable,dimension(:) :: mpi_buf_output

      if(mpi_rank==0) then
        allocate(mpi_buf_recv(mpi_total_nCells))
        allocate(mpi_buf_output(mpi_total_nCells))
        call MPI_GATHERV(input,mpi_nCells,mpi_real_kind,mpi_buf_recv,mpi_cell_indexes_all_ordered_counts,&
                          mpi_cell_indexes_all_ordered_displs, mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        call mpi_oncell_output_real_adjust_1d(mpi_buf_recv,mpi_buf_output)
        call fiona_output(dataset_name, var_name, mpi_buf_output)
        deallocate(mpi_buf_recv)
        deallocate(mpi_buf_output)
      else
        call MPI_GATHERV(input,mpi_nCells,mpi_real_kind,mpi_buf_recv,mpi_cell_indexes_all_ordered_counts,&
                          mpi_cell_indexes_all_ordered_displs, mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)

      end if

      
    end subroutine

    ! gather 2 dimension of data and output into netcdf files
    subroutine mpi_output_oncell_2d_real(dataset_name,var_name,input)
      character(*), intent(IN) :: dataset_name
      character(*), intent(IN) :: var_name
      real(kind=real_kind),intent(IN),dimension(:,:) :: input
      integer :: dimen1,i
      real(kind=real_kind),allocatable,dimension(:) :: mpi_buf_recv
      real(kind=real_kind),allocatable,dimension(:,:) :: mpi_buf_output
      integer,allocatable,dimension(:) :: mpi_indexes_all_ordered_counts_2d
      integer,allocatable,dimension(:) :: mpi_indexes_all_ordered_displs_2d
      
      dimen1=SIZE(input,2)

      if(mpi_rank==0) then
        allocate(mpi_buf_recv(mpi_total_nCells*dimen1))
        allocate(mpi_buf_output(mpi_total_nCells,dimen1))
        allocate(mpi_indexes_all_ordered_counts_2d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_2d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_2d(i)= mpi_cell_indexes_all_ordered_counts(i)* dimen1
        end do
        mpi_indexes_all_ordered_displs_2d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_2d(i)= mpi_indexes_all_ordered_displs_2d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_2d(i - 1)
        end do
        call MPI_GATHERV(input,mpi_nCells*dimen1,mpi_real_kind,mpi_buf_recv,mpi_indexes_all_ordered_counts_2d,&
                          mpi_indexes_all_ordered_displs_2d, mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        call mpi_oncell_output_real_adjust_2d(mpi_buf_recv,mpi_buf_output)
        call fiona_output(dataset_name, var_name, mpi_buf_output)
        deallocate(mpi_buf_recv)
        deallocate(mpi_buf_output)
        deallocate(mpi_indexes_all_ordered_counts_2d)
        deallocate(mpi_indexes_all_ordered_displs_2d)
      else
        call MPI_GATHERV(input,mpi_nCells*dimen1,mpi_real_kind,mpi_buf_recv,mpi_indexes_all_ordered_counts_2d,&
                          mpi_indexes_all_ordered_displs_2d, mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)

      end if
    end subroutine

    ! gather 3 dimension of data and output into netcdf files
    subroutine mpi_output_oncell_3d_real(dataset_name,var_name,input)
      character(*), intent(IN) :: dataset_name
      character(*), intent(IN) :: var_name
      real(kind=real_kind),intent(IN),dimension(:,:,:) :: input
      integer :: dimen1,dimen2,i
      real(kind=real_kind),allocatable,dimension(:) :: mpi_buf_recv
      real(kind=real_kind),allocatable,dimension(:,:,:) :: mpi_buf_output
      integer,allocatable,dimension(:) :: mpi_indexes_all_ordered_counts_3d
      integer,allocatable,dimension(:) :: mpi_indexes_all_ordered_displs_3d
      
      dimen1=SIZE(input,2)
      dimen2=SIZE(input,3)

      if(mpi_rank==0) then
        allocate(mpi_buf_recv(mpi_total_nCells*dimen1*dimen2))
        allocate(mpi_buf_output(mpi_total_nCells,dimen1,dimen2))
        allocate(mpi_indexes_all_ordered_counts_3d(mpi_procs))
        allocate(mpi_indexes_all_ordered_displs_3d(mpi_procs))
        do i = 1,mpi_procs
          mpi_indexes_all_ordered_counts_3d(i)= mpi_cell_indexes_all_ordered_counts(i)* dimen1*dimen2
        end do
        mpi_indexes_all_ordered_displs_3d(1)= 0
        do i = 2,mpi_procs
          mpi_indexes_all_ordered_displs_3d(i)= mpi_indexes_all_ordered_displs_3d(i - 1) + &
                                                     mpi_indexes_all_ordered_counts_3d(i - 1)
        end do
        call MPI_GATHERV(input,mpi_nCells*dimen1*dimen2,mpi_real_kind,mpi_buf_recv,mpi_indexes_all_ordered_counts_3d,&
                          mpi_indexes_all_ordered_displs_3d, mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)
        call mpi_oncell_output_real_adjust_3d(mpi_buf_recv,mpi_buf_output)
        call fiona_output(dataset_name, var_name, mpi_buf_output)
        deallocate(mpi_buf_recv)
        deallocate(mpi_buf_output)
        deallocate(mpi_indexes_all_ordered_counts_3d)
        deallocate(mpi_indexes_all_ordered_displs_3d)
      else
        call MPI_GATHERV(input,mpi_nCells*dimen1*dimen2,mpi_real_kind,mpi_buf_recv,mpi_indexes_all_ordered_counts_3d,&
                          mpi_indexes_all_ordered_displs_3d, mpi_real_kind,0,MPI_COMM_WORLD,mpi_err)

      end if
    end subroutine

    ! gather 1 dimension of data and output into netcdf files
    subroutine mpi_output_oncell_1d_integer(dataset_name,var_name,input)
      character(*), intent(IN) :: dataset_name
      character(*), intent(IN) :: var_name
      integer,intent(IN),dimension(:) :: input
      integer,allocatable,dimension(:) :: mpi_buf_recv
      integer,allocatable,dimension(:) :: mpi_buf_output

      if(mpi_rank==0) then
        allocate(mpi_buf_recv(mpi_total_nCells))
        allocate(mpi_buf_output(mpi_total_nCells))
        call MPI_GATHERV(input,mpi_nCells,MPI_INTEGER,mpi_buf_recv,mpi_cell_indexes_all_ordered_counts,&
                          mpi_cell_indexes_all_ordered_displs, MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
        call mpi_oncell_output_integer_adjust_1d(mpi_buf_recv,mpi_buf_output)
        call fiona_output(dataset_name, var_name, mpi_buf_output)
        deallocate(mpi_buf_recv)
        deallocate(mpi_buf_output)
      else
        call MPI_GATHERV(input,mpi_nCells,MPI_INTEGER,mpi_buf_recv,mpi_cell_indexes_all_ordered_counts,&
                          mpi_cell_indexes_all_ordered_displs, MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)

      end if

      
    end subroutine

    !adjust the buffer from received data order to output order
    subroutine mpi_oncell_output_real_adjust_1d(buf_recv,buf_output)
      real(kind=real_kind),dimension(:),intent(IN) :: buf_recv
      real(kind=real_kind),dimension(:),intent(INOUT) :: buf_output
      integer :: i

      do i=1,mpi_total_nCells
        buf_output(mpi_cell_indexes_all_ordered(i))=buf_recv(i)
      end do
        
    end subroutine

    !adjust the buffer from received data order to output order
    subroutine mpi_oncell_output_real_adjust_2d(buf_recv,buf_output)
      real(kind=real_kind),dimension(:),intent(IN) :: buf_recv
      real(kind=real_kind),dimension(:,:),intent(INOUT) :: buf_output
      integer :: i,j,k,dimen1,start_recv,start_output

      dimen1=SIZE(buf_output,2)
      start_recv=1
      do i=1,mpi_procs
        start_output=mpi_cell_indexes_all_ordered_displs(i)
        do j=1,dimen1
          do k=1,mpi_cell_indexes_all_ordered_counts(i)
            buf_output(mpi_cell_indexes_all_ordered(start_output+k),j)=buf_recv(start_recv)
            start_recv=start_recv+1
          end do
        end do
      end do
    end subroutine

    !adjust the buffer from received data order to output order
    subroutine mpi_oncell_output_real_adjust_3d(buf_recv,buf_output)
      real(kind=real_kind),dimension(:),intent(IN) :: buf_recv
      real(kind=real_kind),dimension(:,:,:),intent(INOUT) :: buf_output
      integer :: i,j,k,m,dimen1,dimen2,start_recv,start_output

      dimen1=SIZE(buf_output,2)
      dimen2=SIZE(buf_output,3)
      start_recv=1
      do i=1,mpi_procs
        start_output=mpi_cell_indexes_all_ordered_displs(i)
        do m=1,dimen2
          do j=1,dimen1
            do k=1,mpi_cell_indexes_all_ordered_counts(i)
              buf_output(mpi_cell_indexes_all_ordered(start_output+k),j,m)=buf_recv(start_recv)
              start_recv=start_recv+1
            end do
          end do
        end do
      end do
    end subroutine

    !adjust the buffer from received data order to output order
    subroutine mpi_oncell_input_real_adjust_2d(buf_input,buf_output)
      real(kind=real_kind),dimension(:,:),intent(IN) :: buf_input
      real(kind=real_kind),dimension(:),intent(INOUT) :: buf_output
      integer :: i,j,k,dimen1,start_input,start_output

      dimen1=SIZE(buf_input,2)
      start_output=1
      do i=1,mpi_procs
        start_input=mpi_cell_indexes_all_ordered_displs(i)
        do j=1,dimen1
          do k=1,mpi_cell_indexes_all_ordered_counts(i)
            buf_output(start_output)=buf_input(mpi_cell_indexes_all_ordered(start_input+k),j)
            start_output=start_output+1
          end do
        end do
      end do
    end subroutine

    !adjust the buffer from received data order to output order
    subroutine mpi_oncell_input_real_adjust_3d(buf_input,buf_output)
      real(kind=real_kind),dimension(:,:,:),intent(IN) :: buf_input
      real(kind=real_kind),dimension(:),intent(INOUT) :: buf_output
      integer :: i,j,k,m,dimen1,dimen2,start_input,start_output

      dimen1=SIZE(buf_input,2)
      dimen2=SIZE(buf_input,3)
      start_output=1
      do i=1,mpi_procs
        start_input=mpi_cell_indexes_all_ordered_displs(i)
        do m=1,dimen2
          do j=1,dimen1
            do k=1,mpi_cell_indexes_all_ordered_counts(i)
              buf_output(start_output)=buf_input(mpi_cell_indexes_all_ordered(start_input+k),j,m)
              start_output=start_output+1
            end do
          end do
        end do
      end do
    end subroutine

    !adjust the buffer from received data order to output order based on edge
    subroutine mpi_onedge_output_real_adjust_1d(buf_recv,buf_output)
      real(kind=real_kind),dimension(:),intent(IN) :: buf_recv
      real(kind=real_kind),dimension(:),intent(INOUT) :: buf_output
      integer :: i

      do i=1,mpi_edge_total_ordered_counts
        buf_output(mpi_edge_indexes_all_ordered(i))=buf_recv(i)
      end do
        
    end subroutine

    !adjust the buffer from received data order to output order based on input station 
    subroutine mpi_adjust_sta_1d_real(buf_recv,buf_output)
      real(kind=real_kind),dimension(:),intent(IN) :: buf_recv
      real(kind=real_kind),dimension(:),intent(INOUT) :: buf_output
      integer :: i

      do i=1,mpi_nSta_sum
        buf_output(mpi_nSta_index_all(i))=buf_recv(i)
      end do
        
    end subroutine

    !adjust the buffer from received data order to output order
    subroutine mpi_oncell_output_integer_adjust_1d(buf_recv,buf_output)
      integer,dimension(:),intent(IN) :: buf_recv
      integer,dimension(:),intent(INOUT) :: buf_output
      integer :: i

      do i=1,mpi_total_nCells
        buf_output(mpi_cell_indexes_all_ordered(i))=buf_recv(i)
      end do
        
    end subroutine

    ! sort edges or vertices on each process, delete duplicated edges or vertices
    subroutine mpi_get_index_edges_vertices(edges_vertices_oncell,dimen1_oncell,dimen2_oncell,size_edges_vertices,total_size_edges_vertices,mode)
      integer, intent(in),dimension(:,:) :: edges_vertices_oncell
      integer, intent(in) :: dimen1_oncell,dimen2_oncell
      integer, intent(inout) :: size_edges_vertices
      integer, intent(in) :: total_size_edges_vertices
      character(len = 6),intent(in) :: mode
      integer, allocatable,dimension(:) :: index_edges_vertices
      integer :: i,j,pos
      allocate(index_edges_vertices(total_size_edges_vertices))
      index_edges_vertices = 0
      size_edges_vertices = 0
      do i = 1,dimen2_oncell
        do j = 1,dimen1_oncell
          if (edges_vertices_oncell(j,i)> 0) then
            if(index_edges_vertices(edges_vertices_oncell(j,i)) == 0) then
              index_edges_vertices(edges_vertices_oncell(j,i))= edges_vertices_oncell(j,i)
              size_edges_vertices = size_edges_vertices + 1
            end if
          end if
        end do
      end do

      select case (mode)
        case ("onedge")
          mpi_nEdges = size_edges_vertices
          allocate(mpi_edge_indexes(size_edges_vertices))
          pos = 1
          do i = 1,total_size_edges_vertices
            if(index_edges_vertices(i)> 0) then
              mpi_edge_indexes(pos)= index_edges_vertices(i)
              pos = pos + 1
            end if
          end do

          ! collect index of edges to process 0
          if (mpi_rank == 0) then
            allocate(mpi_edge_indexes_all_ordered_counts(mpi_procs))
            allocate(mpi_edge_indexes_all_ordered_displs(mpi_procs))

            call MPI_GATHER(size_edges_vertices, 1, MPI_INTEGER, mpi_edge_indexes_all_ordered_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            mpi_edge_indexes_all_ordered_displs(1) = 0                                        

            do i = 2, mpi_procs                                             
              mpi_edge_indexes_all_ordered_displs(i) = mpi_edge_indexes_all_ordered_displs(i - 1) + &
                                                       mpi_edge_indexes_all_ordered_counts(i - 1) 
            end do
            mpi_edge_total_ordered_counts = mpi_edge_indexes_all_ordered_displs(mpi_procs)+ mpi_edge_indexes_all_ordered_counts(mpi_procs)

            allocate(mpi_edge_indexes_all_ordered(mpi_edge_total_ordered_counts))

            call MPI_GATHERV(mpi_edge_indexes, size_edges_vertices, MPI_INTEGER, mpi_edge_indexes_all_ordered, mpi_edge_indexes_all_ordered_counts, &
                             mpi_edge_indexes_all_ordered_displs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

          else
            call MPI_GATHER(size_edges_vertices, 1, MPI_INTEGER, mpi_edge_indexes_all_ordered_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            call MPI_GATHERV(mpi_edge_indexes, size_edges_vertices, MPI_INTEGER, mpi_edge_indexes_all_ordered, mpi_edge_indexes_all_ordered_counts, &
                             mpi_edge_indexes_all_ordered_displs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
          end if

        case ("onvert")
          allocate(mpi_vertex_indexes(size_edges_vertices))
          pos = 1
          do i = 1,total_size_edges_vertices
            if(index_edges_vertices(i)> 0) then
              mpi_vertex_indexes(pos)= index_edges_vertices(i)
              pos = pos + 1
            end if
          end do

          ! collect index of vertex to process 0
          if (mpi_rank == 0) then
            allocate(mpi_vertex_indexes_all_ordered_counts(mpi_procs))
            allocate(mpi_vertex_indexes_all_ordered_displs(mpi_procs))

            call MPI_GATHER(size_edges_vertices, 1, MPI_INTEGER, mpi_vertex_indexes_all_ordered_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            mpi_vertex_indexes_all_ordered_displs(1) = 0                                        

            do i = 2, mpi_procs                                             
              mpi_vertex_indexes_all_ordered_displs(i) = mpi_vertex_indexes_all_ordered_displs(i - 1) + &
                                                       mpi_vertex_indexes_all_ordered_counts(i - 1) 
            end do
            mpi_vertex_total_ordered_counts = mpi_vertex_indexes_all_ordered_displs(mpi_procs)+ mpi_vertex_indexes_all_ordered_counts(mpi_procs)

            allocate(mpi_vertex_indexes_all_ordered(mpi_vertex_total_ordered_counts))

            call MPI_GATHERV(mpi_vertex_indexes, size_edges_vertices, MPI_INTEGER, mpi_vertex_indexes_all_ordered, mpi_vertex_indexes_all_ordered_counts, &
                             mpi_vertex_indexes_all_ordered_displs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

          else
            call MPI_GATHER(size_edges_vertices, 1, MPI_INTEGER, mpi_vertex_indexes_all_ordered_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

            call MPI_GATHERV(mpi_vertex_indexes, size_edges_vertices, MPI_INTEGER, mpi_vertex_indexes_all_ordered, mpi_vertex_indexes_all_ordered_counts, &
                             mpi_vertex_indexes_all_ordered_displs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
          end if

        case DEFAULT
          write(*,*) "Error: wrong mode option for mpi_get_index_edges_vertices"
          stop
      end select

      deallocate(index_edges_vertices)
    end subroutine

    ! find halo cells on each process
    subroutine mpi_find_halo_cells(cellsOnCell)
      integer, dimension(:,:),intent(IN) ::  cellsOnCell
      integer :: maxEdges
      integer :: nCells
      integer :: i,j,k,l,m,cell_ID,proc_ID
      logical :: found
      integer,allocatable,dimension(:,:) :: cells_on_process
      integer,allocatable,dimension(:,:) :: halo_level
      integer,allocatable,dimension(:,:) :: halo_start
      integer,allocatable,dimension(:) :: cells_num_process
      maxEdges = SIZE(cellsOnCell,1)
      nCells = SIZE(cellsOnCell,2)
      allocate(cells_on_process(nCells,mpi_procs))
      allocate(halo_start(mpi_halo_level,mpi_procs))
      allocate(halo_level(nCells,mpi_procs))
      allocate(cells_num_process(mpi_procs))
      
      cells_num_process = 0
      halo_level = 0

!#ifdef MPI_DEBUG
!      ! only for debug
!      allocate(mpi_cellsOnCell(maxEdges,nCells))
!      mpi_cellsOnCell(:,:) = cellsOnCell(:,:)
!#endif

      ! add first layer of halo 
      do i = 1,nCells
        do j = 1,maxEdges
          cell_ID = cellsOnCell(j,i)
          if(cell_ID > 0) then
            proc_ID = mpi_cell_indexes_all(cell_ID)
            if(proc_ID .ne. mpi_rank) then
              found = .FALSE.
              do k = 1,cells_num_process(proc_ID + 1)
                ! find value, do nothing
                if(cell_ID == cells_on_process(k,proc_ID + 1)) then
                  found =.TRUE.
                  exit
                end if
              end do
              if(.not. found) then
                cells_num_process(proc_ID + 1)= cells_num_process(proc_ID + 1)+ 1
                cells_on_process(cells_num_process(proc_ID + 1),proc_ID + 1)= cell_ID
                halo_level(cells_num_process(proc_ID + 1),proc_ID + 1)= 1
              end if
            end if
          end if
        end do
      end do

      ! add other layeres of halo 
      halo_start = 1
      ! if mpi_halo_level >= 2,then keeping adding cell ID
      do i = 2,mpi_halo_level
        do j = 1,mpi_procs
        ! if mpi_halo_level = 2,adding cell ID from mpi_halo = 1 and if mpi_halo = 3,adding cell ID from mpi_halo = 2
          halo_start(i,j)= cells_num_process(j)+ 1
          do k = halo_start(i - 1,j),cells_num_process(j)
            do l = 1,maxEdges
              cell_ID = mpi_cellsOnCell_all(l,cells_on_process(k,j))
              if(cell_ID > 0) then
                proc_ID = mpi_cell_indexes_all(cell_ID)
                if( proc_ID .ne. mpi_rank) then
                  found =.FALSE.
                  do m = 1,cells_num_process(proc_ID + 1)
                    ! find value, do nothing
                    if(cell_ID == cells_on_process(m,proc_ID + 1)) then
                      found =.TRUE.
                      exit
                    end if
                  end do
                  if(.not. found) then
                    cells_num_process(proc_ID + 1)= cells_num_process(proc_ID + 1)+ 1
                    cells_on_process(cells_num_process(proc_ID + 1),proc_ID + 1)= cell_ID
                    halo_level(cells_num_process(proc_ID + 1),proc_ID + 1)= i
                  end if
                end if
              end if
            end do
          end do
        end do
      end do


      mpi_num_halocell = 0
      mpi_graph_indegree = 0
      do i = 1,mpi_procs
        if(cells_num_process(i) > 0) then
          mpi_graph_indegree = mpi_graph_indegree + 1
          mpi_num_halocell = mpi_num_halocell + cells_num_process(i)
        end if
      end do

      allocate(mpi_graph_sources(mpi_graph_indegree))
      allocate(mpi_cell_recv_indexes_1d_displs(mpi_graph_indegree))

      mpi_num_ncells_halo = mpi_num_halocell + nCells  
      call reallocate(mpi_cell_indexes,mpi_num_ncells_halo)
      allocate(mpi_halo_level_cell(mpi_num_halocell))
      allocate(mpi_cell_recv_indexes_1d_counts(mpi_graph_indegree))
      
      j = 1
      k = 1
      do i = 1,mpi_procs
        if(cells_num_process(i)> 0) then
          mpi_graph_sources(k)= i - 1
          mpi_cell_recv_indexes_1d_counts(k)= cells_num_process(i)
          k = k + 1
          call mpi_sort(cells_on_process(1:cells_num_process(i),i),halo_level(1:cells_num_process(i),i))
          mpi_cell_indexes(nCells + j:nCells + j + cells_num_process(i) - 1)= cells_on_process(1:cells_num_process(i),i)
          mpi_halo_level_cell(j:j + cells_num_process(i) - 1)= halo_level(1:cells_num_process(i),i)
          j = j + cells_num_process(i)
        end if
      end do

      mpi_cell_recv_indexes_1d_displs(1)=0
      do i=2,mpi_graph_indegree
        mpi_cell_recv_indexes_1d_displs(i)=mpi_cell_recv_indexes_1d_displs(i-1)+mpi_cell_recv_indexes_1d_counts(i-1)
      end do

    end subroutine

    ! find halo edges on each process
    subroutine mpi_find_halo_edges(edgesOnCell,nEdges)
      integer, dimension(:,:),intent(IN) ::  edgesOnCell
      integer, intent(IN) :: nEdges
      integer :: maxEdges,nCells
      integer :: i,j,k,edge_ID,proc_ID,cell_ID
      logical :: found
      integer,allocatable,dimension(:,:) :: edges_on_process
      integer,allocatable,dimension(:,:) :: halo_level
      integer,allocatable,dimension(:,:) :: halo_start
      integer,allocatable,dimension(:) :: edges_num_process
      maxEdges = SIZE(edgesOnCell,1)
      nCells = SIZE(edgesOnCell,2)
      allocate(edges_on_process(maxEdges * mpi_num_halocell,mpi_procs))
      allocate(halo_level(maxEdges * mpi_num_halocell,mpi_procs))
      allocate(edges_num_process(mpi_procs))
      edges_num_process = 0
      halo_level = 0

      ! add local indexes of edges

      ! add first layer of halo 
      do i = 1,mpi_num_halocell
        cell_ID = mpi_cell_indexes(nCells + i)
        proc_ID = mpi_cell_indexes_all(cell_ID)
        do j = 1,maxEdges
          edge_ID = mpi_edgesOnCell_all(j,cell_ID)
          if(edge_ID > 0) then
            if(mpi_findloc_ordered(mpi_edge_indexes,edge_ID) .eq. 0) then
              found =.FALSE.
              do k = 1,edges_num_process(proc_ID + 1)
                if(edge_ID == edges_on_process(k,proc_ID + 1)) then
                  found =.TRUE.
                  exit
                end if
              end do

              if(.not. found) then
                edges_num_process(proc_ID + 1)= edges_num_process(proc_ID + 1)+ 1
                edges_on_process(edges_num_process(proc_ID + 1),proc_ID + 1)= edge_ID
                halo_level(edges_num_process(proc_ID + 1),proc_ID + 1)= mpi_halo_level_cell(i)
              end if

            end if
          end if
        end do
      end do


      mpi_num_haloedge = 0
      mpi_edge_indegree=0
      do i = 1,mpi_procs
        if(edges_num_process(i) > 0) then
          mpi_edge_indegree=mpi_edge_indegree+1
          mpi_num_haloedge = mpi_num_haloedge + edges_num_process(i)
        end if
      end do

      mpi_num_nEdges_halo = nEdges + mpi_num_haloedge
      call reallocate(mpi_edge_indexes,mpi_num_nEdges_halo)
      allocate(mpi_halo_level_edge(mpi_num_haloedge))
      
      allocate(mpi_edge_recv_indexes_1d_counts(mpi_edge_indegree))
      allocate(mpi_edge_sources(mpi_edge_indegree))
      
      j = 1
      k = 1
      do i = 1,mpi_procs
        if(edges_num_process(i)> 0) then
          mpi_edge_recv_indexes_1d_counts(k)= edges_num_process(i)
          mpi_edge_sources(k)=i-1
          k = k + 1
          call mpi_sort(edges_on_process(1:edges_num_process(i),i),halo_level(1:edges_num_process(i),i))
          mpi_edge_indexes(nEdges + j:nEdges + j + edges_num_process(i) - 1)= edges_on_process(1:edges_num_process(i),i)
          mpi_halo_level_edge(j:j + edges_num_process(i) - 1)= halo_level(1:edges_num_process(i),i)
          j = j + edges_num_process(i)
        end if
      end do
      
      allocate(mpi_edge_recv_indexes_1d_displs(mpi_edge_indegree))
      mpi_edge_recv_indexes_1d_displs(1)= 0
      do i = 2,mpi_edge_indegree
        mpi_edge_recv_indexes_1d_displs(i)= mpi_edge_recv_indexes_1d_displs(i - 1)+ mpi_edge_recv_indexes_1d_counts(i - 1)
      end do

      deallocate(edges_on_process)
      deallocate(halo_level)
      deallocate(edges_num_process)

    end subroutine

    ! find index position of first value which is equal to the searching value    
    INTEGER FUNCTION mpi_findloc(a, val) RESULT(i)                                    
       INTEGER, DIMENSION(:) :: a                                                     
       INTEGER :: val                                                                 
       !$acc routine seq
       do i = 1, size(a)                                                              
          if (a(i) == val) exit                                                       
       enddo                                                                          
       if (i > size(a)) i = 0                                                         
       return                                                                         
    END FUNCTION mpi_findloc

    ! find index position of first value which is equal to the searching value    
    INTEGER FUNCTION mpi_findloc_ordered(a, val) RESULT(j)                                    
       INTEGER, DIMENSION(:) :: a                                                     
       INTEGER :: val,i                                                                 
       INTEGER :: test                                                                
       test = size(a)                                                                 
       do i = 1, size(a)                                                              
          if (a(i) == val) then
            j = i
            return
          end if
          if (a(i) > val) then
            j = 0
            return
          end if
       enddo                                                                          
       j = 0                                                         
       return                                                                         
    END FUNCTION 

    SUBROUTINE mpi_sort_one(a)                                                         
                                                                                      
      IMPLICIT NONE                                                                     
      INTEGER :: i, j, increment                                                        
      INTEGER :: temp                                                           
      INTEGER, INTENT(in out) :: a(:)                                             
                                                                                        
      increment = SIZE(a) / 2                                                           
      DO WHILE (increment > 0)                                                          
          DO i = increment + 1, SIZE(a)                                                   
             j = i                                                                      
             temp = a(i)                                                                
             DO WHILE (j >= increment + 1 .AND. a(j - increment) > temp)                    
                a(j) = a(j - increment)                                                   
                j = j - increment                                                       
             END DO                                                                     
             a(j) = temp                                                                
          END DO                                                                        
                                                                                        
          IF (increment == 2) THEN                                                      
             increment = 1                                                              
          ELSE                                                                          
             increment = increment * 5 / 11                                             
          END IF                                                                        
      END DO                                                                            
                                                                                        
    END SUBROUTINE 

    SUBROUTINE mpi_sort_two(a, b)                                                         
                                                                                        
      IMPLICIT NONE                                                                     
      INTEGER :: i, j, increment                                                        
      INTEGER :: temp, temp_b                                                           
      INTEGER, INTENT(inout) :: a(:), b(:)                                             
                                                                                        
      increment = SIZE(a) / 2                                                           
      DO WHILE (increment > 0)                                                          
          DO i = increment + 1, SIZE(a)                                                   
             j = i                                                                      
             temp = a(i)                                                                
             temp_b = b(i)                                                              
!             DO WHILE (j >= increment + 1 .AND. a(j - increment) > temp)                    
             DO WHILE (j >= increment + 1)                    
               IF(a(j - increment) > temp) THEN
                 a(j) = a(j - increment)                                                   
                 b(j) = b(j-increment)                                                   
                 j = j - increment             
               ELSE
                 EXIT
               END IF
             END DO                                                                     
             a(j) = temp                                                                
             b(j) = temp_b                                                              
          END DO                                                                        
                                                                                        
          IF (increment == 2) THEN                                                      
             increment = 1                                                              
          ELSE                                                                          
             increment = increment * 5 / 11                                             
          END IF                                                                        
      END DO                                                                            
                                                                                        
    END SUBROUTINE 

    subroutine mpi_bcast_2d(array)
      integer,intent(INOUT),dimension(:,:) :: array
      integer :: dimen1,dimen2
      dimen1 = SIZE(array,1)
      dimen2 = SIZE(array,2)
      call MPI_BCAST(array, dimen1 * dimen2, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    end subroutine

    subroutine mpi_get_edge_destination(nEdges)
      integer, intent(IN) :: nEdges
!      integer :: i
!      allocate(mpi_edge_send_indexes_1d_counts(mpi_graph_outdegree))
!      CALL MPI_NEIGHBOR_ALLTOALL(mpi_edge_recv_indexes_1d_counts, 1, MPI_INTEGER, mpi_edge_send_indexes_1d_counts, &
!                                  1, MPI_INTEGER, mpi_graph_comm, mpi_err) 
!      
!      allocate(mpi_edge_send_indexes_1d_displs(mpi_graph_outdegree))
!        mpi_edge_send_indexes_1d_displs(1) = 0
!        mpi_edge_send_indexes_1d_counts_sum = mpi_edge_send_indexes_1d_counts(1)
!      do i = 2,mpi_graph_outdegree
!        mpi_edge_send_indexes_1d_displs(i)= mpi_edge_send_indexes_1d_displs(i - 1)+ mpi_edge_send_indexes_1d_counts(i - 1)
!        mpi_edge_send_indexes_1d_counts_sum = mpi_edge_send_indexes_1d_counts_sum + mpi_edge_send_indexes_1d_counts(i)
!      end do
!      
!      allocate(mpi_edge_send_indexes_1d(mpi_edge_send_indexes_1d_counts_sum))
!      
!      CALL MPI_NEIGHBOR_ALLTOALLV(mpi_edge_indexes(nEdges + 1), mpi_edge_recv_indexes_1d_counts, mpi_edge_recv_indexes_1d_displs, &
!                                  MPI_INTEGER, mpi_edge_send_indexes_1d, mpi_edge_send_indexes_1d_counts, mpi_edge_send_indexes_1d_displs, &
!                                  MPI_INTEGER, mpi_graph_comm, mpi_err)
!      
      INTEGER :: i, k, dest_process,dest_position,start,dest_postion_end,source_start,source_end
      INTEGER :: mpi_edge_indegree_counts_sum
      INTEGER :: mpi_edge_outdegree_counts_sum
      INTEGER :: mpi_edge_total_recv_indexes_all_current_p
      INTEGER :: mpi_edge_recv_indexes_1d_counts_sum_all_sum
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_indegree_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_recv_indexes_1d_counts_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_outdegree_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_outdegree_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_send_indexes_1d_counts_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_dests_all_current_p
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_send_indexes_1d_all_current_p
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_send_indexes_1d_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_recv_indexes_1d_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_recv_indexes_1d_counts_sum_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_recv_indexes_1d_counts_sum_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_send_indexes_1d_counts_sum_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_edge_send_indexes_1d_counts_sum_displs_all

      IF (mpi_rank == 0) THEN

         ! collecting edge information in process 0, then scatter them to other processes
         ALLOCATE (mpi_edge_indegree_all(mpi_procs))
         ALLOCATE (mpi_edge_indegree_displs_all(mpi_procs))
         CALL MPI_GATHER(mpi_edge_indegree, 1, MPI_INTEGER, mpi_edge_indegree_all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         ! collecting neighbouring process number for each process
         mpi_edge_indegree_counts_sum = mpi_edge_indegree_all(1)
         mpi_edge_indegree_displs_all(1) = 0
         DO i = 2, mpi_procs
            mpi_edge_indegree_counts_sum = mpi_edge_indegree_counts_sum + mpi_edge_indegree_all(i)
            mpi_edge_indegree_displs_all(i) = mpi_edge_indegree_displs_all(i - 1) + mpi_edge_indegree_all(i - 1)
         END DO
         ALLOCATE (mpi_edge_sources_all(mpi_edge_indegree_counts_sum))
         CALL MPI_GATHERV(mpi_edge_sources, mpi_edge_indegree, MPI_INTEGER, mpi_edge_sources_all, mpi_edge_indegree_all, &
                          mpi_edge_indegree_displs_all, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         ALLOCATE (mpi_edge_recv_indexes_1d_counts_all(mpi_edge_indegree_counts_sum))
         CALL MPI_GATHERV(mpi_edge_recv_indexes_1d_counts, mpi_edge_indegree, MPI_INTEGER, mpi_edge_recv_indexes_1d_counts_all, &
                          mpi_edge_indegree_all, mpi_edge_indegree_displs_all, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         ! get number of destitaton processes
         ALLOCATE (mpi_edge_outdegree_all(mpi_procs))
         mpi_edge_outdegree_all = 0
         ! total indegrees are equal to total outdegrees
         mpi_edge_outdegree_counts_sum = mpi_edge_indegree_counts_sum
         DO i = 1, mpi_edge_outdegree_counts_sum
           mpi_edge_outdegree_all(mpi_edge_sources_all(i) + 1) = mpi_edge_outdegree_all(mpi_edge_sources_all(i)+ 1)+ 1
         END DO
         CALL MPI_SCATTER(mpi_edge_outdegree_all, 1, MPI_INTEGER, mpi_edge_outdegree, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,mpi_err)

         ! according to process order, put all destitaton process index and sending index counts together in the following codes
         ALLOCATE (mpi_edge_outdegree_displs_all(mpi_procs))
         mpi_edge_outdegree_displs_all(1) = 0
         DO i = 2, mpi_procs
            mpi_edge_outdegree_displs_all(i) = mpi_edge_outdegree_displs_all(i - 1) + mpi_edge_outdegree_all(i - 1)
         END DO
         ALLOCATE (mpi_edge_dests_all(mpi_edge_outdegree_counts_sum))
         ALLOCATE (mpi_edge_dests_all_current_p(mpi_procs))

         ! according to received process ID and received number of grids, get sending process ID and sending number of grids
         ALLOCATE (mpi_edge_send_indexes_1d_counts_all(mpi_edge_outdegree_counts_sum))
         ALLOCATE(mpi_edge_recv_indexes_1d_counts_sum_all(mpi_procs))
         mpi_edge_recv_indexes_1d_counts_sum_all = 0
         mpi_edge_dests_all_current_p = 1
         start = 1
         DO i = 1,mpi_procs
           DO k = 1,mpi_edge_indegree_all(i)
             dest_process = mpi_edge_sources_all(start) + 1
             dest_position = mpi_edge_outdegree_displs_all(dest_process) + mpi_edge_dests_all_current_p(dest_process)
             mpi_edge_dests_all(dest_position) = i - 1
             mpi_edge_send_indexes_1d_counts_all(dest_position) = mpi_edge_recv_indexes_1d_counts_all(start)
             mpi_edge_recv_indexes_1d_counts_sum_all(i) = mpi_edge_recv_indexes_1d_counts_sum_all(i) + mpi_edge_recv_indexes_1d_counts_all(start)
             start = start + 1
             mpi_edge_dests_all_current_p(dest_process) = mpi_edge_dests_all_current_p(dest_process) + 1
           END DO
         END DO

         ALLOCATE (mpi_edge_dests(mpi_edge_outdegree))
         CALL MPI_SCATTERV(mpi_edge_dests_all, mpi_edge_outdegree_all, mpi_edge_outdegree_displs_all, MPI_INTEGER, &
                            mpi_edge_dests, mpi_edge_outdegree, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         ALLOCATE (mpi_edge_send_indexes_1d_counts(mpi_edge_outdegree))
         CALL MPI_SCATTERV(mpi_edge_send_indexes_1d_counts_all, mpi_edge_outdegree_all, mpi_edge_outdegree_displs_all, &
                            MPI_INTEGER, mpi_edge_send_indexes_1d_counts, mpi_edge_outdegree, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

        ! gather receiving index ID
        ALLOCATE(mpi_edge_recv_indexes_1d_counts_sum_displs_all(mpi_procs))
        mpi_edge_recv_indexes_1d_counts_sum_displs_all(1) = 0
        mpi_edge_recv_indexes_1d_counts_sum_all_sum = mpi_edge_recv_indexes_1d_counts_sum_all(1)
        DO i = 2,mpi_procs
          mpi_edge_recv_indexes_1d_counts_sum_displs_all(i)= mpi_edge_recv_indexes_1d_counts_sum_displs_all(i - 1)+ &
            mpi_edge_recv_indexes_1d_counts_sum_all(i - 1)
          mpi_edge_recv_indexes_1d_counts_sum_all_sum = mpi_edge_recv_indexes_1d_counts_sum_all_sum + mpi_edge_recv_indexes_1d_counts_sum_all(i)
        END DO
        ALLOCATE(mpi_edge_recv_indexes_1d_all(mpi_edge_recv_indexes_1d_counts_sum_all_sum))
        CALL MPI_GATHERV(mpi_edge_indexes(nedges + 1), mpi_num_haloedge, MPI_INTEGER, mpi_edge_recv_indexes_1d_all, &
                          mpi_edge_recv_indexes_1d_counts_sum_all, mpi_edge_recv_indexes_1d_counts_sum_displs_all, MPI_INTEGER, 0, &
                          MPI_COMM_WORLD, mpi_err)


        ! according to received process ID and index ID to get sending index ID
        ALLOCATE(mpi_edge_send_indexes_1d_counts_sum_all(mpi_procs))
        ALLOCATE(mpi_edge_send_indexes_1d_counts_sum_displs_all(mpi_procs))
        start = 1
        mpi_edge_send_indexes_1d_counts_sum_all = 0
        DO k = 1,mpi_edge_outdegree_all(1)
          mpi_edge_send_indexes_1d_counts_sum_all(1) = mpi_edge_send_indexes_1d_counts_sum_all(1)+ mpi_edge_send_indexes_1d_counts_all(start)
          start = start + 1
        END DO
        mpi_edge_send_indexes_1d_counts_sum_displs_all(1) = 0
        DO i = 2,mpi_procs
          DO k = 1,mpi_edge_outdegree_all(i)
            mpi_edge_send_indexes_1d_counts_sum_all(i) = mpi_edge_send_indexes_1d_counts_sum_all(i)+ mpi_edge_send_indexes_1d_counts_all(start)
            start = start + 1
          END DO
          mpi_edge_send_indexes_1d_counts_sum_displs_all(i)= mpi_edge_send_indexes_1d_counts_sum_displs_all(i - 1) + &
            mpi_edge_send_indexes_1d_counts_sum_all(i - 1)
        END DO
        ALLOCATE(mpi_edge_send_indexes_1d_all(mpi_edge_recv_indexes_1d_counts_sum_all_sum))
        ALLOCATE(mpi_edge_send_indexes_1d_all_current_p(mpi_procs))
        mpi_edge_send_indexes_1d_all_current_p = 1
        start = 1
        source_start = 1
        DO i = 1,mpi_procs
          DO k = 1,mpi_edge_indegree_all(i)
            dest_process = mpi_edge_sources_all(start) + 1
            dest_position = mpi_edge_send_indexes_1d_counts_sum_displs_all(dest_process) + &
              mpi_edge_send_indexes_1d_all_current_p(dest_process)
            dest_postion_end = dest_position + mpi_edge_recv_indexes_1d_counts_all(start) - 1
            source_end = source_start + mpi_edge_recv_indexes_1d_counts_all(start) - 1
            mpi_edge_send_indexes_1d_all(dest_position:dest_postion_end)= mpi_edge_recv_indexes_1d_all(source_start:source_end)
            mpi_edge_send_indexes_1d_all_current_p(dest_process) = mpi_edge_send_indexes_1d_all_current_p(dest_process) + &
              mpi_edge_recv_indexes_1d_counts_all(start)
            start = start + 1
            source_start = source_end + 1
          END DO
        END DO
        ! initialize number and displacement of send indexes
        ALLOCATE(mpi_edge_send_indexes_1d_displs(mpi_edge_outdegree))
        mpi_edge_send_indexes_1d_displs(1) = 0
        mpi_edge_send_indexes_1d_counts_sum = mpi_edge_send_indexes_1d_counts(1)
        DO i = 2, mpi_edge_outdegree
           mpi_edge_send_indexes_1d_displs(i) = mpi_edge_send_indexes_1d_displs(i - 1) + mpi_edge_send_indexes_1d_counts(i - 1)
           mpi_edge_send_indexes_1d_counts_sum = mpi_edge_send_indexes_1d_counts_sum + mpi_edge_send_indexes_1d_counts(i)
        END DO
        ! scatter send index ID
        ALLOCATE(mpi_edge_send_indexes_1d(mpi_edge_send_indexes_1d_counts_sum))
        CALL MPI_SCATTERV(mpi_edge_send_indexes_1d_all, mpi_edge_send_indexes_1d_counts_sum_all, mpi_edge_send_indexes_1d_counts_sum_displs_all, &
                           MPI_INTEGER, mpi_edge_send_indexes_1d, mpi_edge_send_indexes_1d_counts_sum, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         DEALLOCATE(mpi_edge_outdegree_all)
         DEALLOCATE(mpi_edge_outdegree_displs_all)
         DEALLOCATE(mpi_edge_send_indexes_1d_counts_all)
         DEALLOCATE(mpi_edge_dests_all_current_p)
         DEALLOCATE(mpi_edge_indegree_displs_all)
         DEALLOCATE(mpi_edge_recv_indexes_1d_counts_all)
         DEALLOCATE(mpi_edge_recv_indexes_1d_counts_sum_all)
         DEALLOCATE(mpi_edge_recv_indexes_1d_counts_sum_displs_all)
         DEALLOCATE(mpi_edge_recv_indexes_1d_all)
         DEALLOCATE(mpi_edge_send_indexes_1d_all_current_p)
         DEALLOCATE(mpi_edge_send_indexes_1d_all)
         DEALLOCATE(mpi_edge_send_indexes_1d_counts_sum_all)

      ELSE

         CALL MPI_GATHER(mpi_edge_indegree, 1, MPI_INTEGER, mpi_edge_indegree_all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         CALL MPI_GATHERV(mpi_edge_sources, mpi_edge_indegree, MPI_INTEGER, mpi_edge_sources_all, mpi_edge_indegree_all, &
                          mpi_edge_indegree_displs_all, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         CALL MPI_GATHERV(mpi_edge_recv_indexes_1d_counts, mpi_edge_indegree, MPI_INTEGER, mpi_edge_recv_indexes_1d_counts_all, &
                          mpi_edge_indegree_all, mpi_edge_indegree_displs_all, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         CALL MPI_SCATTER(mpi_edge_outdegree_all, mpi_procs, MPI_INTEGER, mpi_edge_outdegree, 1, &
                          MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         ALLOCATE (mpi_edge_dests(mpi_edge_outdegree))
         CALL MPI_SCATTERV(mpi_edge_dests_all, mpi_edge_outdegree_all, mpi_edge_outdegree_displs_all, MPI_INTEGER, &
                            mpi_edge_dests, mpi_edge_outdegree, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         ALLOCATE (mpi_edge_send_indexes_1d_counts(mpi_edge_outdegree))
         CALL MPI_SCATTERV(mpi_edge_send_indexes_1d_counts_all, mpi_edge_outdegree_all, &
                           mpi_edge_outdegree_displs_all, MPI_INTEGER, mpi_edge_send_indexes_1d_counts, &
                           mpi_edge_outdegree, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         CALL MPI_GATHERV(mpi_edge_indexes(nedges + 1), mpi_num_haloedge, MPI_INTEGER, mpi_edge_recv_indexes_1d_all, &
                          mpi_edge_recv_indexes_1d_counts_sum_all, mpi_edge_recv_indexes_1d_counts_sum_displs_all, MPI_INTEGER, 0, &
                          MPI_COMM_WORLD, mpi_err)

         ! initialize number and displacement of send indexes
         ALLOCATE(mpi_edge_send_indexes_1d_displs(mpi_edge_outdegree))
         mpi_edge_send_indexes_1d_displs(1) = 0
         mpi_edge_send_indexes_1d_counts_sum = mpi_edge_send_indexes_1d_counts(1)
         DO i = 2, mpi_edge_outdegree
            mpi_edge_send_indexes_1d_displs(i) = mpi_edge_send_indexes_1d_displs(i - 1) + mpi_edge_send_indexes_1d_counts(i - 1)
            mpi_edge_send_indexes_1d_counts_sum = mpi_edge_send_indexes_1d_counts_sum + mpi_edge_send_indexes_1d_counts(i)
         END DO
        ! scatter send index ID
        ALLOCATE(mpi_edge_send_indexes_1d(mpi_edge_send_indexes_1d_counts_sum))
        CALL MPI_SCATTERV(mpi_edge_send_indexes_1d_all, mpi_edge_send_indexes_1d_counts_sum_all, mpi_edge_send_indexes_1d_counts_sum_displs_all, &
                           MPI_INTEGER, mpi_edge_send_indexes_1d, mpi_edge_send_indexes_1d_counts_sum, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      END IF
    end subroutine
    
    ! adjust global cell id to local cell id
    subroutine mpi_adjust_cell_indexes(array)
      integer,intent(INOUT),dimension(:,:) :: array
      integer :: dimen1
      integer :: dimen2
      integer :: i,j,new_value
      dimen1 = SIZE(array,1)
      dimen2 = SIZE(array,2)
      !$acc kernels loop collapse(2) present(array) copyin(mpi_cell_indexes)
      do i = 1,dimen2
        do j = 1,dimen1
          if(array(j,i)> 0) then
            new_value = mpi_findloc(mpi_cell_indexes,array(j,i))
            if(new_value > 0) then
              array(j,i)= new_value
            else
!              write(*,*)"Error: cannot find cell id =",array(j,i)," in local and halo buffer in process =",mpi_rank
              stop
            end if
          end if
        end do
      end do
      !$acc end kernels
    end subroutine

    ! adjust global edge id to edge cell id
    subroutine mpi_adjust_edge_indexes(array)
      integer,intent(INOUT),dimension(:,:) :: array
      integer :: dimen1
      integer :: dimen2
      integer :: i,j,new_value
      dimen1 = SIZE(array,1)
      dimen2 = SIZE(array,2)
      !$acc kernels loop collapse(2) present(array) copyin(mpi_edge_indexes)
      do i = 1,dimen2
        do j = 1,dimen1
          if(array(j,i)> 0) then
            new_value = mpi_findloc(mpi_edge_indexes,array(j,i))
            if(new_value > 0) then
              array(j,i)= new_value
            else
!              write(*,*)"Error: cannot find edge id =",array(j,i)," in local and halo buffer in process =",mpi_rank
              stop
            end if
          end if
        end do
      end do
      !$acc end kernels
    end subroutine

    ! searching destination process ID for all computing processes
    SUBROUTINE mpi_get_cell_destination(nCells)
      integer,intent(in) :: nCells
      INTEGER :: i, k, dest_process,dest_position,start,dest_postion_end,source_start,source_end
      INTEGER :: mpi_graph_indegree_counts_sum
      INTEGER :: mpi_graph_outdegree_counts_sum
      INTEGER :: mpi_graph_total_recv_indexes_all_current_p
      INTEGER :: mpi_cell_recv_indexes_1d_counts_sum_all_sum
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_graph_indegree_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_recv_indexes_1d_counts_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_graph_outdegree_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_graph_outdegree_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d_counts_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_graph_dests_all_current_p
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d_all_current_p
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_recv_indexes_1d_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_recv_indexes_1d_counts_sum_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_recv_indexes_1d_counts_sum_displs_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d_counts_sum_all
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d_counts_sum_displs_all

      IF (mpi_rank == 0) THEN

         ! collecting graph information in process 0, then scatter them to other processes
         ALLOCATE (mpi_graph_indegree_all(mpi_procs))
         ALLOCATE (mpi_graph_indegree_displs_all(mpi_procs))
         CALL MPI_GATHER(mpi_graph_indegree, 1, MPI_INTEGER, mpi_graph_indegree_all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         ! collecting neighbouring process number for each process
         mpi_graph_indegree_counts_sum = mpi_graph_indegree_all(1)
         mpi_graph_indegree_displs_all(1) = 0
         DO i = 2, mpi_procs
            mpi_graph_indegree_counts_sum = mpi_graph_indegree_counts_sum + mpi_graph_indegree_all(i)
            mpi_graph_indegree_displs_all(i) = mpi_graph_indegree_displs_all(i - 1) + mpi_graph_indegree_all(i - 1)
         END DO
         ALLOCATE (mpi_graph_sources_all(mpi_graph_indegree_counts_sum))
         CALL MPI_GATHERV(mpi_graph_sources, mpi_graph_indegree, MPI_INTEGER, mpi_graph_sources_all, mpi_graph_indegree_all, &
                          mpi_graph_indegree_displs_all, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         ALLOCATE (mpi_cell_recv_indexes_1d_counts_all(mpi_graph_indegree_counts_sum))
         CALL MPI_GATHERV(mpi_cell_recv_indexes_1d_counts, mpi_graph_indegree, MPI_INTEGER, mpi_cell_recv_indexes_1d_counts_all, &
                          mpi_graph_indegree_all, mpi_graph_indegree_displs_all, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         ! get number of destitaton processes
         ALLOCATE (mpi_graph_outdegree_all(mpi_procs))
         mpi_graph_outdegree_all = 0
         ! total indegrees are equal to total outdegrees
         mpi_graph_outdegree_counts_sum = mpi_graph_indegree_counts_sum
         DO i = 1, mpi_graph_outdegree_counts_sum
           mpi_graph_outdegree_all(mpi_graph_sources_all(i) + 1) = mpi_graph_outdegree_all(mpi_graph_sources_all(i)+ 1)+ 1
         END DO
         CALL MPI_SCATTER(mpi_graph_outdegree_all, 1, MPI_INTEGER, mpi_graph_outdegree, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,mpi_err)

         ! according to process order, put all destitaton process index and sending index counts together in the following codes
         ALLOCATE (mpi_graph_outdegree_displs_all(mpi_procs))
         mpi_graph_outdegree_displs_all(1) = 0
         DO i = 2, mpi_procs
            mpi_graph_outdegree_displs_all(i) = mpi_graph_outdegree_displs_all(i - 1) + mpi_graph_outdegree_all(i - 1)
         END DO
         ALLOCATE (mpi_graph_dests_all(mpi_graph_outdegree_counts_sum))
         ALLOCATE (mpi_graph_dests_all_current_p(mpi_procs))

         ! according to received process ID and received number of grids, get sending process ID and sending number of grids
         ALLOCATE (mpi_cell_send_indexes_1d_counts_all(mpi_graph_outdegree_counts_sum))
         ALLOCATE(mpi_cell_recv_indexes_1d_counts_sum_all(mpi_procs))
         mpi_cell_recv_indexes_1d_counts_sum_all = 0
         mpi_graph_dests_all_current_p = 1
         start = 1
         DO i = 1,mpi_procs
           DO k = 1,mpi_graph_indegree_all(i)
             dest_process = mpi_graph_sources_all(start) + 1
             dest_position = mpi_graph_outdegree_displs_all(dest_process) + mpi_graph_dests_all_current_p(dest_process)
             mpi_graph_dests_all(dest_position) = i - 1
             mpi_cell_send_indexes_1d_counts_all(dest_position) = mpi_cell_recv_indexes_1d_counts_all(start)
             mpi_cell_recv_indexes_1d_counts_sum_all(i) = mpi_cell_recv_indexes_1d_counts_sum_all(i) + mpi_cell_recv_indexes_1d_counts_all(start)
             start = start + 1
             mpi_graph_dests_all_current_p(dest_process) = mpi_graph_dests_all_current_p(dest_process) + 1
           END DO
         END DO

         ALLOCATE (mpi_graph_dests(mpi_graph_outdegree))
         CALL MPI_SCATTERV(mpi_graph_dests_all, mpi_graph_outdegree_all, mpi_graph_outdegree_displs_all, MPI_INTEGER, &
                            mpi_graph_dests, mpi_graph_outdegree, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         ALLOCATE (mpi_cell_send_indexes_1d_counts(mpi_graph_outdegree))
         CALL MPI_SCATTERV(mpi_cell_send_indexes_1d_counts_all, mpi_graph_outdegree_all, mpi_graph_outdegree_displs_all, &
                            MPI_INTEGER, mpi_cell_send_indexes_1d_counts, mpi_graph_outdegree, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

        ! gather receiving index ID
        ALLOCATE(mpi_cell_recv_indexes_1d_counts_sum_displs_all(mpi_procs))
        mpi_cell_recv_indexes_1d_counts_sum_displs_all(1) = 0
        mpi_cell_recv_indexes_1d_counts_sum_all_sum = mpi_cell_recv_indexes_1d_counts_sum_all(1)
        DO i = 2,mpi_procs
          mpi_cell_recv_indexes_1d_counts_sum_displs_all(i)= mpi_cell_recv_indexes_1d_counts_sum_displs_all(i - 1)+ &
            mpi_cell_recv_indexes_1d_counts_sum_all(i - 1)
          mpi_cell_recv_indexes_1d_counts_sum_all_sum = mpi_cell_recv_indexes_1d_counts_sum_all_sum + mpi_cell_recv_indexes_1d_counts_sum_all(i)
        END DO
        ALLOCATE(mpi_cell_recv_indexes_1d_all(mpi_cell_recv_indexes_1d_counts_sum_all_sum))
        CALL MPI_GATHERV(mpi_cell_indexes(nCells + 1), mpi_num_halocell, MPI_INTEGER, mpi_cell_recv_indexes_1d_all, &
                          mpi_cell_recv_indexes_1d_counts_sum_all, mpi_cell_recv_indexes_1d_counts_sum_displs_all, MPI_INTEGER, 0, &
                          MPI_COMM_WORLD, mpi_err)


        ! according to received process ID and index ID to get sending index ID
        ALLOCATE(mpi_cell_send_indexes_1d_counts_sum_all(mpi_procs))
        ALLOCATE(mpi_cell_send_indexes_1d_counts_sum_displs_all(mpi_procs))
        start = 1
        mpi_cell_send_indexes_1d_counts_sum_all = 0
        DO k = 1,mpi_graph_outdegree_all(1)
          mpi_cell_send_indexes_1d_counts_sum_all(1) = mpi_cell_send_indexes_1d_counts_sum_all(1)+ mpi_cell_send_indexes_1d_counts_all(start)
          start = start + 1
        END DO
        mpi_cell_send_indexes_1d_counts_sum_displs_all(1) = 0
        DO i = 2,mpi_procs
          DO k = 1,mpi_graph_outdegree_all(i)
            mpi_cell_send_indexes_1d_counts_sum_all(i) = mpi_cell_send_indexes_1d_counts_sum_all(i)+ mpi_cell_send_indexes_1d_counts_all(start)
            start = start + 1
          END DO
          mpi_cell_send_indexes_1d_counts_sum_displs_all(i)= mpi_cell_send_indexes_1d_counts_sum_displs_all(i - 1) + &
            mpi_cell_send_indexes_1d_counts_sum_all(i - 1)
        END DO
        ALLOCATE(mpi_cell_send_indexes_1d_all(mpi_cell_recv_indexes_1d_counts_sum_all_sum))
        ALLOCATE(mpi_cell_send_indexes_1d_all_current_p(mpi_procs))
        mpi_cell_send_indexes_1d_all_current_p = 1
        start = 1
        source_start = 1
        DO i = 1,mpi_procs
          DO k = 1,mpi_graph_indegree_all(i)
            dest_process = mpi_graph_sources_all(start) + 1
            dest_position = mpi_cell_send_indexes_1d_counts_sum_displs_all(dest_process) + &
              mpi_cell_send_indexes_1d_all_current_p(dest_process)
            dest_postion_end = dest_position + mpi_cell_recv_indexes_1d_counts_all(start) - 1
            source_end = source_start + mpi_cell_recv_indexes_1d_counts_all(start) - 1
            mpi_cell_send_indexes_1d_all(dest_position:dest_postion_end)= mpi_cell_recv_indexes_1d_all(source_start:source_end)
            mpi_cell_send_indexes_1d_all_current_p(dest_process) = mpi_cell_send_indexes_1d_all_current_p(dest_process) + &
              mpi_cell_recv_indexes_1d_counts_all(start)
            start = start + 1
            source_start = source_end + 1
          END DO
        END DO
        ! initialize number and displacement of send indexes
        ALLOCATE(mpi_cell_send_indexes_1d_displs(mpi_graph_outdegree))
        mpi_cell_send_indexes_1d_displs(1) = 0
        mpi_cell_send_indexes_1d_counts_sum = mpi_cell_send_indexes_1d_counts(1)
        DO i = 2, mpi_graph_outdegree
           mpi_cell_send_indexes_1d_displs(i) = mpi_cell_send_indexes_1d_displs(i - 1) + mpi_cell_send_indexes_1d_counts(i - 1)
           mpi_cell_send_indexes_1d_counts_sum = mpi_cell_send_indexes_1d_counts_sum + mpi_cell_send_indexes_1d_counts(i)
        END DO
        ! scatter send index ID
        ALLOCATE(mpi_cell_send_indexes_1d(mpi_cell_send_indexes_1d_counts_sum))
        CALL MPI_SCATTERV(mpi_cell_send_indexes_1d_all, mpi_cell_send_indexes_1d_counts_sum_all, mpi_cell_send_indexes_1d_counts_sum_displs_all, &
                           MPI_INTEGER, mpi_cell_send_indexes_1d, mpi_cell_send_indexes_1d_counts_sum, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         DEALLOCATE(mpi_graph_outdegree_all)
         DEALLOCATE(mpi_graph_outdegree_displs_all)
         DEALLOCATE(mpi_cell_send_indexes_1d_counts_all)
         DEALLOCATE(mpi_graph_dests_all_current_p)
         DEALLOCATE(mpi_graph_indegree_displs_all)
         DEALLOCATE(mpi_cell_recv_indexes_1d_counts_all)
         DEALLOCATE(mpi_cell_recv_indexes_1d_counts_sum_all)
         DEALLOCATE(mpi_cell_recv_indexes_1d_counts_sum_displs_all)
         DEALLOCATE(mpi_cell_recv_indexes_1d_all)
         DEALLOCATE(mpi_cell_send_indexes_1d_all_current_p)
         DEALLOCATE(mpi_cell_send_indexes_1d_all)
         DEALLOCATE(mpi_cell_send_indexes_1d_counts_sum_all)

      ELSE

         CALL MPI_GATHER(mpi_graph_indegree, 1, MPI_INTEGER, mpi_graph_indegree_all, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         CALL MPI_GATHERV(mpi_graph_sources, mpi_graph_indegree, MPI_INTEGER, mpi_graph_sources_all, mpi_graph_indegree_all, &
                          mpi_graph_indegree_displs_all, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         CALL MPI_GATHERV(mpi_cell_recv_indexes_1d_counts, mpi_graph_indegree, MPI_INTEGER, mpi_cell_recv_indexes_1d_counts_all, &
                          mpi_graph_indegree_all, mpi_graph_indegree_displs_all, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         CALL MPI_SCATTER(mpi_graph_outdegree_all, mpi_procs, MPI_INTEGER, mpi_graph_outdegree, 1, &
                          MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         ALLOCATE (mpi_graph_dests(mpi_graph_outdegree))
         CALL MPI_SCATTERV(mpi_graph_dests_all, mpi_graph_outdegree_all, mpi_graph_outdegree_displs_all, MPI_INTEGER, &
                            mpi_graph_dests, mpi_graph_outdegree, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
         ALLOCATE (mpi_cell_send_indexes_1d_counts(mpi_graph_outdegree))
         CALL MPI_SCATTERV(mpi_cell_send_indexes_1d_counts_all, mpi_graph_outdegree_all, &
                           mpi_graph_outdegree_displs_all, MPI_INTEGER, mpi_cell_send_indexes_1d_counts, &
                           mpi_graph_outdegree, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)

         CALL MPI_GATHERV(mpi_cell_indexes(nCells + 1), mpi_num_halocell, MPI_INTEGER, mpi_cell_recv_indexes_1d_all, &
                          mpi_cell_recv_indexes_1d_counts_sum_all, mpi_cell_recv_indexes_1d_counts_sum_displs_all, MPI_INTEGER, 0, &
                          MPI_COMM_WORLD, mpi_err)

         ! initialize number and displacement of send indexes
         ALLOCATE(mpi_cell_send_indexes_1d_displs(mpi_graph_outdegree))
         mpi_cell_send_indexes_1d_displs(1) = 0
         mpi_cell_send_indexes_1d_counts_sum = mpi_cell_send_indexes_1d_counts(1)
         DO i = 2, mpi_graph_outdegree
            mpi_cell_send_indexes_1d_displs(i) = mpi_cell_send_indexes_1d_displs(i - 1) + mpi_cell_send_indexes_1d_counts(i - 1)
            mpi_cell_send_indexes_1d_counts_sum = mpi_cell_send_indexes_1d_counts_sum + mpi_cell_send_indexes_1d_counts(i)
         END DO
        ! scatter send index ID
        ALLOCATE(mpi_cell_send_indexes_1d(mpi_cell_send_indexes_1d_counts_sum))
        CALL MPI_SCATTERV(mpi_cell_send_indexes_1d_all, mpi_cell_send_indexes_1d_counts_sum_all, mpi_cell_send_indexes_1d_counts_sum_displs_all, &
                           MPI_INTEGER, mpi_cell_send_indexes_1d, mpi_cell_send_indexes_1d_counts_sum, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
      END IF


    END SUBROUTINE 

    subroutine mpi_init_graph_index

    ! create graph topology                                                    
    CALL MPI_DIST_GRAPH_CREATE_ADJACENT(MPI_COMM_WORLD, mpi_graph_indegree, mpi_graph_sources, MPI_UNWEIGHTED, &
                                        mpi_graph_outdegree, mpi_graph_dests, MPI_UNWEIGHTED, &
                                        MPI_INFO_NULL, .true., mpi_graph_comm, mpi_err)
      
    end subroutine

   ! blocking data exchange for wp of 1 dimension
   SUBROUTINE mpi_real_1d_block_exchange_oncell(mpi_recv_data)
      real(kind=real_kind), DIMENSION(:), INTENT(INOUT) :: mpi_recv_data
      CALL mpi_real_1d_prepare_sendbuf(mpi_recv_data, mpi_real_cell_sendbuf_1d)
      CALL MPI_NEIGHBOR_ALLTOALLV(mpi_real_cell_sendbuf_1d, mpi_cell_send_indexes_1d_counts, mpi_cell_send_indexes_1d_displs, mpi_real_kind, &
                                  mpi_recv_data(mpi_nCells + 1), mpi_cell_recv_indexes_1d_counts, mpi_cell_recv_indexes_1d_displs, &
                                  mpi_real_kind, mpi_graph_comm, mpi_err)

   END SUBROUTINE
   
   ! blocking data exchange for integer of 1 dimension
   SUBROUTINE mpi_int_1d_block_exchange_oncell(mpi_recv_data)
      integer, DIMENSION(:), INTENT(INOUT) :: mpi_recv_data
      CALL mpi_int_1d_prepare_sendbuf(mpi_recv_data, mpi_int_cell_sendbuf_1d)
      CALL MPI_NEIGHBOR_ALLTOALLV(mpi_int_cell_sendbuf_1d, mpi_cell_send_indexes_1d_counts, mpi_cell_send_indexes_1d_displs, MPI_INTEGER, &
                                  mpi_recv_data(mpi_nCells + 1), mpi_cell_recv_indexes_1d_counts, mpi_cell_recv_indexes_1d_displs, &
                                  MPI_INTEGER, mpi_graph_comm, mpi_err)

   END SUBROUTINE

   ! blocking data exchange for integer of 1 dimension
   SUBROUTINE mpi_int_1d_block_exchange_onedge(mpi_recv_data)
      integer, DIMENSION(:), INTENT(INOUT) :: mpi_recv_data
      CALL mpi_int_1d_prepare_sendbuf(mpi_recv_data, mpi_int_cell_sendbuf_1d)
      CALL MPI_NEIGHBOR_ALLTOALLV(mpi_int_cell_sendbuf_1d, mpi_cell_send_indexes_1d_counts, mpi_cell_send_indexes_1d_displs, MPI_INTEGER, &
                                  mpi_recv_data(mpi_nCells + 1), mpi_cell_recv_indexes_1d_counts, mpi_cell_recv_indexes_1d_displs, &
                                  MPI_INTEGER, mpi_graph_comm, mpi_err)

   END SUBROUTINE

   ! blocking data exchange for real type of 2 dimension
   SUBROUTINE mpi_real_2d_block_exchange_oncell(mpi_recv_data)
      real(kind=real_kind), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_data
      INTEGER :: dimen1
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs_tmp
      real(kind=real_kind), ALLOCATABLE, DIMENSION(:) :: mpi_send_buf_tmp
      real(kind=real_kind), ALLOCATABLE, DIMENSION(:) :: mpi_recv_buf_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_displs_tmp
      dimen1=SIZE(mpi_recv_data,2)
      ALLOCATE (mpi_send_indexes_1d_counts_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_indexes_1d_displs_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_buf_tmp(mpi_cell_send_indexes_1d_counts_sum * dimen1))
      ALLOCATE (mpi_recv_buf_tmp(mpi_num_halocell * dimen1))
      ALLOCATE (mpi_recv_indexes_counts_tmp(mpi_graph_indegree))
      ALLOCATE (mpi_recv_indexes_displs_tmp(mpi_graph_indegree))

      CALL mpi_real_2d_prepare_sendbuf(mpi_recv_data, mpi_send_buf_tmp, dimen1, mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp)
      CALL mpi_2d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1)
      CALL MPI_NEIGHBOR_ALLTOALLV(mpi_send_buf_tmp, mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp, mpi_real_kind, &
                                  mpi_recv_buf_tmp, mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, mpi_real_kind, &
                                  mpi_graph_comm, mpi_err)
      CALL mpi_real_2d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1)

      DEALLOCATE (mpi_send_indexes_1d_counts_tmp)
      DEALLOCATE (mpi_send_indexes_1d_displs_tmp)
      DEALLOCATE (mpi_recv_buf_tmp)
      DEALLOCATE (mpi_send_buf_tmp)
      DEALLOCATE (mpi_recv_indexes_counts_tmp)
      DEALLOCATE (mpi_recv_indexes_displs_tmp)

   END SUBROUTINE 

   ! blocking data exchange for real type of 3 dimension
   SUBROUTINE mpi_real_3d_block_exchange_oncell(mpi_recv_data)
      real(kind=real_kind), DIMENSION(:, :,:), INTENT(INOUT) :: mpi_recv_data
      INTEGER :: dimen1,dimen2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs_tmp
      real(kind=real_kind), ALLOCATABLE, DIMENSION(:) :: mpi_send_buf_tmp
      real(kind=real_kind), ALLOCATABLE, DIMENSION(:) :: mpi_recv_buf_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_displs_tmp
      dimen1=SIZE(mpi_recv_data,2)
      dimen2=SIZE(mpi_recv_data,3)
      ALLOCATE (mpi_send_indexes_1d_counts_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_indexes_1d_displs_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_buf_tmp(mpi_cell_send_indexes_1d_counts_sum * dimen1* dimen2))
      ALLOCATE (mpi_recv_buf_tmp(mpi_num_halocell * dimen1*dimen2))
      ALLOCATE (mpi_recv_indexes_counts_tmp(mpi_graph_indegree))
      ALLOCATE (mpi_recv_indexes_displs_tmp(mpi_graph_indegree))

      CALL mpi_real_3d_prepare_sendbuf(mpi_recv_data, mpi_send_buf_tmp, dimen1, dimen2,mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp)
      CALL mpi_3d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1,dimen2)
      CALL SYSTEM_CLOCK(mpi_clock_neigh_commu_start)                                    
      CALL MPI_NEIGHBOR_ALLTOALLV(mpi_send_buf_tmp, mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp, mpi_real_kind, &
                                  mpi_recv_buf_tmp, mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, mpi_real_kind, &
                                  mpi_graph_comm, mpi_err)
      CALL SYSTEM_CLOCK(mpi_clock_neigh_commu_end) 
      mpi_clock_neigh_commu=mpi_clock_neigh_commu+mpi_clock_neigh_commu_end-mpi_clock_neigh_commu_start
      CALL mpi_real_3d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1,dimen2)
      DEALLOCATE (mpi_send_indexes_1d_counts_tmp)
      DEALLOCATE (mpi_send_indexes_1d_displs_tmp)
      DEALLOCATE (mpi_recv_buf_tmp)
      DEALLOCATE (mpi_send_buf_tmp)
      DEALLOCATE (mpi_recv_indexes_counts_tmp)
      DEALLOCATE (mpi_recv_indexes_displs_tmp)

   END SUBROUTINE 

   ! blocking data exchange for real type of 3 dimension
   SUBROUTINE mpi_real_3d_block_exchange_oncell_onebyone(mpi_recv_data)
      real(kind=real_kind), DIMENSION(:, :,:), INTENT(INOUT) :: mpi_recv_data
      INTEGER :: dimen1,dimen2,i
      INTEGER(int_double) :: start
      INTEGER :: mpi_status(MPI_STATUS_SIZE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_send_indexes_1d_displs_tmp
      real(kind=real_kind), ALLOCATABLE, DIMENSION(:) :: mpi_send_buf_tmp
      real(kind=real_kind), ALLOCATABLE, DIMENSION(:) :: mpi_recv_buf_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_counts_tmp
      INTEGER, ALLOCATABLE, DIMENSION(:) :: mpi_recv_indexes_displs_tmp
      INTEGER, ALLOCATABLE,DIMENSION(:) ::  mpi_req
      dimen1=SIZE(mpi_recv_data,2)
      dimen2=SIZE(mpi_recv_data,3)
      ALLOCATE (mpi_send_indexes_1d_counts_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_send_indexes_1d_displs_tmp(mpi_graph_outdegree))
      ALLOCATE (mpi_req(mpi_graph_outdegree))
      ALLOCATE (mpi_send_buf_tmp(INT(mpi_cell_send_indexes_1d_counts_sum * dimen1* dimen2,int_double)))
      ALLOCATE (mpi_recv_indexes_counts_tmp(mpi_graph_indegree))
      ALLOCATE (mpi_recv_indexes_displs_tmp(mpi_graph_indegree))

!      write(*,*)"end mpi_real_3d_prepare_sendbuf"
      CALL mpi_real_3d_prepare_sendbuf(mpi_recv_data, mpi_send_buf_tmp, dimen1, dimen2,mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp)
      CALL mpi_3d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1,dimen2)
!      write(*,*)"end mpi_3d_prepare_recv_counts_displs"
      ALLOCATE (mpi_recv_buf_tmp(maxval(mpi_cell_recv_indexes_1d_counts) * dimen1*dimen2))
      mpi_clock_prep_recv=0
      CALL SYSTEM_CLOCK(mpi_clock_neigh_commu_start)                                    
      DO i=1,mpi_graph_outdegree
        start=mpi_send_indexes_1d_displs_tmp(i)+1
        CALL MPI_ISEND(mpi_send_buf_tmp(start),mpi_send_indexes_1d_counts_tmp(i),mpi_real_kind,mpi_graph_dests(i),&
                        100,MPI_COMM_WORLD,mpi_req(i),mpi_err)
      END DO
!      write(*,*)"end MPI_ISEND"
      DO i=1,mpi_graph_indegree
        CALL MPI_RECV(mpi_recv_buf_tmp,mpi_cell_recv_indexes_1d_counts(i)*dimen1*dimen2,mpi_real_kind,&
                        mpi_graph_sources(i),100,MPI_COMM_WORLD,mpi_status,mpi_err)
        CALL SYSTEM_CLOCK(mpi_clock_prep_recv_start)                                    
        CALL mpi_real_3d_prepare_recvbuf_onebyone(mpi_recv_buf_tmp, mpi_recv_data, dimen1,dimen2,i)
        CALL SYSTEM_CLOCK(mpi_clock_prep_recv_end)                                    
        mpi_clock_prep_recv=mpi_clock_prep_recv+mpi_clock_prep_recv_end-mpi_clock_prep_recv_start
      END DO
!      write(*,*)"end MPI_RECV and mpi_real_3d_prepare_recvbuf_onebyone"
       CALL SYSTEM_CLOCK(mpi_clock_neigh_commu_end) 
       mpi_clock_neigh_commu=mpi_clock_neigh_commu+mpi_clock_neigh_commu_end-mpi_clock_neigh_commu_start-mpi_clock_prep_recv
!      CALL MPI_NEIGHBOR_ALLTOALLV(mpi_send_buf_tmp, mpi_send_indexes_1d_counts_tmp, mpi_send_indexes_1d_displs_tmp, mpi_real_kind, &
!                                  mpi_recv_buf_tmp, mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, mpi_real_kind, &
!                                  mpi_graph_comm, mpi_err)
!
      DEALLOCATE (mpi_send_indexes_1d_counts_tmp)
      DEALLOCATE (mpi_send_indexes_1d_displs_tmp)
      DEALLOCATE (mpi_recv_buf_tmp)
      DEALLOCATE (mpi_send_buf_tmp)
      DEALLOCATE (mpi_recv_indexes_counts_tmp)
      DEALLOCATE (mpi_recv_indexes_displs_tmp)

   END SUBROUTINE 

   ! prepare sending buffer of dp for each process,for all arrays which have same position for collecting sending buffer
   SUBROUTINE mpi_real_1d_prepare_sendbuf(mpi_source_buf, mpi_send_buf)
      real(kind=real_kind), DIMENSION(:), INTENT(IN) :: mpi_source_buf
      real(kind=real_kind), DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER :: i
      DO i = 1, mpi_cell_send_indexes_1d_counts_sum
         mpi_send_buf(i) = mpi_source_buf(mpi_cell_send_indexes_1d(i))
      END DO
   END SUBROUTINE mpi_real_1d_prepare_sendbuf

   ! prepare sending buffer of integer for each process,for all arrays which have same position for collecting sending buffer
   SUBROUTINE mpi_int_1d_prepare_sendbuf(mpi_source_buf, mpi_send_buf)
      integer, DIMENSION(:), INTENT(IN) :: mpi_source_buf
      integer, DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER :: i
      DO i = 1, mpi_cell_send_indexes_1d_counts_sum
         mpi_send_buf(i) = mpi_source_buf(mpi_cell_send_indexes_1d(i))
      END DO
   END SUBROUTINE 

   ! prepare sending buffer of dp with 2 dimensions for each process
   SUBROUTINE mpi_real_2d_prepare_sendbuf(mpi_source_buf, mpi_send_buf, dimen1, mpi_send_indexes_1d_counts_tmp, &
                                        mpi_send_indexes_1d_displs_tmp)
      real(kind=real_kind), DIMENSION(:, :), INTENT(IN) :: mpi_source_buf
      real(kind=real_kind), DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_displs_tmp
      INTEGER :: i, j, k
      INTEGER :: current_p

      mpi_send_indexes_1d_counts_tmp(1) = mpi_cell_send_indexes_1d_counts(1)* dimen1
      mpi_send_indexes_1d_displs_tmp(1) = 0
      DO i = 2, mpi_graph_outdegree
         mpi_send_indexes_1d_counts_tmp(i) = mpi_cell_send_indexes_1d_counts(i)* dimen1
         mpi_send_indexes_1d_displs_tmp(i) = mpi_send_indexes_1d_displs_tmp(i - 1) + mpi_send_indexes_1d_counts_tmp(i - 1)
      END DO

      current_p=0
      DO i = 1, mpi_graph_outdegree
         DO j = 1, dimen1
           DO k = 1, mpi_cell_send_indexes_1d_counts(i)
              current_p=current_p+1
              mpi_send_buf(current_p) = mpi_source_buf(mpi_cell_send_indexes_1d(mpi_cell_send_indexes_1d_displs(i)+k),j)
           END DO
         END DO
      END DO

   END SUBROUTINE 

   ! prepare sending buffer of real type with 3 dimensions for each process
   SUBROUTINE mpi_real_3d_prepare_sendbuf(mpi_source_buf, mpi_send_buf, dimen1, dimen2,mpi_send_indexes_1d_counts_tmp, &
                                        mpi_send_indexes_1d_displs_tmp)
      real(kind=real_kind), DIMENSION(:, :,:), INTENT(IN) :: mpi_source_buf
      real(kind=real_kind), DIMENSION(:), INTENT(INOUT) :: mpi_send_buf
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) :: mpi_send_indexes_1d_displs_tmp
      INTEGER :: i, j, k,m
      INTEGER :: current_p

      mpi_send_indexes_1d_counts_tmp(1) = mpi_cell_send_indexes_1d_counts(1)* dimen1*dimen2
      mpi_send_indexes_1d_displs_tmp(1) = 0
      DO i = 2, mpi_graph_outdegree
         mpi_send_indexes_1d_counts_tmp(i) = mpi_cell_send_indexes_1d_counts(i)* dimen1*dimen2
         mpi_send_indexes_1d_displs_tmp(i) = mpi_send_indexes_1d_displs_tmp(i - 1) + mpi_send_indexes_1d_counts_tmp(i - 1)
      END DO

      current_p=0
      DO i = 1, mpi_graph_outdegree
        DO m=1,dimen2
          DO j = 1, dimen1
            DO k = 1, mpi_cell_send_indexes_1d_counts(i)
               current_p=current_p+1
               mpi_send_buf(current_p) = mpi_source_buf(mpi_cell_send_indexes_1d(mpi_cell_send_indexes_1d_displs(i)+k),j,m)
            END DO
          END DO
        END DO
      END DO

   END SUBROUTINE 

   ! prepare index counts and displacement for data exchange in 2 dimension
   SUBROUTINE mpi_2d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1)
      INTEGER, DIMENSION(:), INTENT(INOUT) ::mpi_recv_indexes_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) ::mpi_recv_indexes_displs_tmp
      INTEGER, INTENT(IN) :: dimen1
      INTEGER :: i
      mpi_recv_indexes_counts_tmp(1) = mpi_cell_recv_indexes_1d_counts(1)* dimen1
      mpi_recv_indexes_displs_tmp(1) = 0
      DO i = 2, mpi_graph_indegree
         mpi_recv_indexes_counts_tmp(i) = mpi_cell_recv_indexes_1d_counts(i)* dimen1
         mpi_recv_indexes_displs_tmp(i) = mpi_recv_indexes_displs_tmp(i - 1) + mpi_recv_indexes_counts_tmp(i - 1)
      END DO

   END SUBROUTINE
   
   ! prepare index counts and displacement for data exchange in 3 dimension
   SUBROUTINE mpi_3d_prepare_recv_counts_displs(mpi_recv_indexes_counts_tmp, mpi_recv_indexes_displs_tmp, dimen1,dimen2)
      INTEGER, DIMENSION(:), INTENT(INOUT) ::mpi_recv_indexes_counts_tmp
      INTEGER, DIMENSION(:), INTENT(INOUT) ::mpi_recv_indexes_displs_tmp
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER :: i
      mpi_recv_indexes_counts_tmp(1) = mpi_cell_recv_indexes_1d_counts(1)* dimen1*dimen2
      mpi_recv_indexes_displs_tmp(1) = 0
      DO i = 2, mpi_graph_indegree
         mpi_recv_indexes_counts_tmp(i) = mpi_cell_recv_indexes_1d_counts(i)* dimen1*dimen2
         mpi_recv_indexes_displs_tmp(i) = mpi_recv_indexes_displs_tmp(i - 1) + mpi_recv_indexes_counts_tmp(i - 1)
      END DO

   END SUBROUTINE

   ! adjust receiving data of dp in the right order in 2 dimension
   SUBROUTINE mpi_real_2d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1)
      real(kind=real_kind), DIMENSION(:), INTENT(IN) :: mpi_recv_buf_tmp
      real(kind=real_kind), DIMENSION(:, :), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER :: i, k, current_p
      current_p = 1
      DO i = 1, mpi_graph_indegree
         DO k = 1, dimen1
            mpi_recv_data(mpi_nCells + mpi_cell_recv_indexes_1d_displs(i) + 1 &
                          :mpi_nCells + mpi_cell_recv_indexes_1d_displs(i) + mpi_cell_recv_indexes_1d_counts(i), k) = &
               mpi_recv_buf_tmp(current_p:current_p + mpi_cell_recv_indexes_1d_counts(i) - 1)
            current_p = current_p + mpi_cell_recv_indexes_1d_counts(i)
         END DO
      END DO
   END SUBROUTINE

   ! adjust receiving data of dp in the right order in 2 dimension
   SUBROUTINE mpi_real_3d_prepare_recvbuf(mpi_recv_buf_tmp, mpi_recv_data, dimen1,dimen2)
      real(kind=real_kind), DIMENSION(:), INTENT(IN) :: mpi_recv_buf_tmp
      real(kind=real_kind), DIMENSION(:, :,:), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER :: i, j,k, current_p
      current_p = 1
      DO i = 1, mpi_graph_indegree
        DO j=1,dimen2
          DO k = 1, dimen1
             mpi_recv_data(mpi_nCells + mpi_cell_recv_indexes_1d_displs(i) + 1 &
                           :mpi_nCells + mpi_cell_recv_indexes_1d_displs(i) + mpi_cell_recv_indexes_1d_counts(i), k,j) = &
                mpi_recv_buf_tmp(current_p:current_p + mpi_cell_recv_indexes_1d_counts(i) - 1)
             current_p = current_p + mpi_cell_recv_indexes_1d_counts(i)
          END DO
        END DO
      END DO
   END SUBROUTINE

   ! adjust receiving data of dp in the right order in 2 dimension
   SUBROUTINE mpi_real_3d_prepare_recvbuf_onebyone(mpi_recv_buf_tmp, mpi_recv_data, dimen1,dimen2,index_degree)
      real(kind=real_kind), DIMENSION(:), INTENT(IN) :: mpi_recv_buf_tmp
      real(kind=real_kind), DIMENSION(:, :,:), INTENT(INOUT) :: mpi_recv_data
      INTEGER, INTENT(IN) :: dimen1
      INTEGER, INTENT(IN) :: dimen2
      INTEGER, INTENT(IN) :: index_degree
      INTEGER :: j,k, current_p
      current_p = 1
      DO j=1,dimen2
        DO k = 1, dimen1
           mpi_recv_data(mpi_nCells + mpi_cell_recv_indexes_1d_displs(index_degree) + 1 &
                         :mpi_nCells + mpi_cell_recv_indexes_1d_displs(index_degree) + mpi_cell_recv_indexes_1d_counts(index_degree), k,j) = &
              mpi_recv_buf_tmp(current_p:current_p + mpi_cell_recv_indexes_1d_counts(index_degree) - 1)
           current_p = current_p + mpi_cell_recv_indexes_1d_counts(index_degree)
        END DO
      END DO
   END SUBROUTINE

   ! initialize arrays used in mpi
   subroutine mpi_init_arrays(nCells,nEdges)
     integer, intent(IN) :: nCells
     integer, intent(IN) :: nEdges
     integer :: i,new_value
     ! initialize buffer for data exchange
     allocate(mpi_real_cell_sendbuf_1d(mpi_cell_send_indexes_1d_counts_sum))
     allocate(mpi_int_cell_sendbuf_1d(mpi_cell_send_indexes_1d_counts_sum))

     do i = 1, mpi_cell_send_indexes_1d_counts_sum
       new_value = mpi_findloc_ordered(mpi_cell_indexes(1:nCells),mpi_cell_send_indexes_1d(i))
       if(new_value > 0) then
         mpi_cell_send_indexes_1d(i)= new_value
       else
         write(*,*) "Error: can not find cell_ID =",mpi_cell_send_indexes_1d(i),"for sending data at processID =",mpi_rank
         stop
       end if
     end do
  
     do i = 1, mpi_edge_send_indexes_1d_counts_sum
       new_value = mpi_findloc_ordered(mpi_edge_indexes(1:nEdges),mpi_edge_send_indexes_1d(i))
       if(new_value > 0) then
         mpi_edge_send_indexes_1d(i)= new_value
       else
         write(*,*)"Error: can not find edge_ID =",mpi_edge_send_indexes_1d(i),"for sending data at processID =",mpi_rank
         stop
       end if
     end do
   end subroutine
   
   ! broadcast key parameters from process 0
   subroutine mpi_bcast_key_parameters(lonTC,latTC,rmwTC,presTC,utTC,vtTC)
     real(kind=real_kind), intent(INOUT), dimension(:) :: lonTC
     real(kind=real_kind), intent(INOUT), dimension(:) :: latTC
     real(kind=real_kind), intent(INOUT), dimension(:) :: rmwTC
     real(kind=real_kind), intent(INOUT), dimension(:) :: presTC
     real(kind=real_kind), intent(INOUT), dimension(:) :: utTC
     real(kind=real_kind), intent(INOUT), dimension(:) :: vtTC
     real(kind=real_kind), allocatable, dimension(:) :: mpi_buf
     integer :: dimen,mpi_buf_size,mpi_pos
     
     dimen=SIZE(lonTC,1)
     
     mpi_buf_size=6*SIZE(lonTC,1)*real_kind

     allocate(mpi_buf(mpi_buf_size))
     mpi_pos=0

     if(mpi_rank==0) then
      call MPI_PACK(lonTC, dimen, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_pos, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(latTC, dimen, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_pos, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(rmwTC, dimen, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_pos, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(presTC, dimen, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_pos, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(utTC, dimen, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_pos, MPI_COMM_WORLD, mpi_err)
      call MPI_PACK(vtTC, dimen, mpi_real_kind, mpi_buf, mpi_buf_size, mpi_pos, MPI_COMM_WORLD, mpi_err)
      
      call MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, 0, MPI_COMM_WORLD, mpi_err)
      deallocate(mpi_buf) 
     else
      
      call MPI_BCAST(mpi_buf, mpi_buf_size, MPI_PACKED, 0, MPI_COMM_WORLD, mpi_err)

      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_pos, lonTC, dimen, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_pos, latTC, dimen, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_pos, rmwTC, dimen, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_pos, presTC, dimen, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_pos, utTC, dimen, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      call MPI_UNPACK(mpi_buf, mpi_buf_size, mpi_pos, vtTC, dimen, mpi_real_kind, MPI_COMM_WORLD, mpi_err)
      
      deallocate(mpi_buf) 
     end if
     
   end subroutine


end module
