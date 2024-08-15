module mod_mpi_variables
  use const_mod
  use mpi
  implicit none

  integer,parameter :: mpi_halo_level = 2 ! how many levels in halo
  INTEGER,PARAMETER  :: DI = MPI_OFFSET_KIND       ! double integer

  integer :: mpi_err ! mpi error parameter for all processes
  integer :: mpi_rank ! process id in all processes
  integer :: mpi_procs ! total number of all processes
  integer :: mpi_real_kind ! total number of all processes
  integer :: mpi_total_nCells ! total number of all cells
  integer :: mpi_total_nEdges ! total number of all edges
  integer :: mpi_total_nVertices ! total number of all vertices
  INTEGER :: mpi_graph_comm ! mpi communicator for all computation in graph for neighbouring communication
  
  integer :: mpi_nCells ! equal to nCells
  integer :: mpi_nEdges ! equal to nCells

  ! source process ID of graph neighbour communication                             
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_sources                  
  ! total number of source process ID of graph neighbour communication,only for process 0
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_total_sources            
  ! total destination process ID of graph neighbour communication,only for process 0
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_total_dests              
  ! destination process ID of graph neighbour communication                        
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_dests                    
  ! number of mpi_graph_sources                                                    
  INTEGER, PUBLIC :: mpi_graph_indegree                                            
  ! number of mpi_graph_total_dests                                                
  INTEGER, PUBLIC :: mpi_graph_outdegree                                           
  ! all indegree in mpi_rank order,only for process 0                              
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_total_indegree           
  ! all indegree in mpi_rank order,only for process 0                              
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_total_outdegree
  ! all indegree in mpi_rank order,only for process 0                           
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_indegree_all          
  ! all indegree in mpi_rank order,only for process 0                           
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_outdegree_all
  ! total number of source process ID of graph neighbour communication,only for process 0
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_sources_all           
  ! total destination process ID of graph neighbour communication,only for process 0
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_graph_dests_all
      
  ! source process ID of edge neighbour communication                             
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_sources                  
  ! total number of source process ID of edge neighbour communication,only for process 0
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_total_sources            
  ! total destination process ID of edge neighbour communication,only for process 0
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_total_dests              
  ! destination process ID of edge neighbour communication                        
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_dests                    
  ! number of mpi_edge_sources                                                    
  INTEGER, PUBLIC :: mpi_edge_indegree                                            
  ! number of mpi_edge_total_dests                                                
  INTEGER, PUBLIC :: mpi_edge_outdegree                                           
  ! all indegree in mpi_rank order,only for process 0                              
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_total_indegree           
  ! all indegree in mpi_rank order,only for process 0                              
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_total_outdegree
  ! all indegree in mpi_rank order,only for process 0                           
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_indegree_all          
  ! all indegree in mpi_rank order,only for process 0                           
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_outdegree_all
  ! total number of source process ID of edge neighbour communication,only for process 0
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_sources_all           
  ! total destination process ID of edge neighbour communication,only for process 0
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_dests_all

  ! global cell index number for cells on each process
  ! will be reallocated from nCells to mpi_num_ncells_halo
  integer, allocatable, dimension(:) :: mpi_cell_indexes 
  integer, allocatable, dimension(:) :: mpi_cell_indexes_all ! global cell index number for cells on each process
  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered ! global cell index number for cells ordered by processes,only for process 0
  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_counts ! number of global cell index for cells on each process,only for process 0
  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_displs ! displacement of global cell index for cells on each process,only for process 0
  integer :: mpi_nSta_sum !sum of nSta, only for process 0
  integer, allocatable, dimension(:) :: mpi_nSta_all !all nSta, only for process 0
  integer, allocatable, dimension(:) :: mpi_nSta_all_displs ! displacement of nSta each process,only for process 0
  integer, allocatable, dimension(:) :: mpi_nSta_index_all ! displacement of nSta each process,only for process 0
 

  ! used for indexes for halo of cells and verticies
  integer :: mpi_num_halocell
  integer :: mpi_num_ncells_halo
  integer, allocatable, dimension(:,:) :: mpi_cellsOnCell_all 
  integer, allocatable, dimension(:,:) :: mpi_edgesOnCell_all 
  ! index starts from nCells + 1
  integer,allocatable,dimension(:) :: mpi_halo_level_cell
  ! pay attention, it starts from nEdges + 1
  integer,allocatable,dimension(:) :: mpi_halo_level_edge
  integer :: mpi_array_size(2)

  ! variables for index of data exchange
  INTEGER, PUBLIC :: mpi_cell_send_indexes_1d_counts_sum ! all number of mpi_send_indexes_1d_counts
  ! INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d ! local 1 dimension indexes of sending
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_recv_total_indexes          
  ! all global 1 dimension indexes of receiving,only for process 0              
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d ! sending counts for each process
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d_counts ! sending counts for each process
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_recv_indexes_1d_counts ! receiving counts for each process
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d_displs ! sending displs for each process
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_recv_indexes_1d_displs ! receiving displs for each process
  ! receiving counts for each process,only for process 0                        
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_recv_total_indexes_counts

  integer :: mpi_edge_total_ordered_counts ! sum of mpi_edge_indexes_all_ordered_counts
  integer, allocatable, dimension(:) :: mpi_edge_indexes ! global edge index number for cells on each process
  integer, allocatable, dimension(:) :: mpi_edge_indexes_all_ordered ! global edge index number for cells ordered by processes ,only for process 0
  integer, allocatable, dimension(:) :: mpi_edge_indexes_all_ordered_counts ! number of global edge index for cells on each process,only for process 0
  integer, allocatable, dimension(:) :: mpi_edge_indexes_all_ordered_displs ! displacement of global edge index for cells on each process,only for process 0

  ! variables for index of data exchange,sum of mpi_edge_send_indexes_1d_counts
  integer :: mpi_num_haloedge
  integer :: mpi_num_nEdges_halo
  ! INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_cell_send_indexes_1d ! local 1 dimension indexes of sending
!  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_recv_total_indexes          
  ! all global 1 dimension indexes of receiving,only for process 0              
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_send_indexes_1d ! sending counts for each process
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_send_indexes_1d_counts ! sending counts for each process
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_recv_indexes_1d_counts ! receiving counts for each process
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_send_indexes_1d_displs ! sending displs for each process
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: mpi_edge_recv_indexes_1d_displs ! receiving displs for each process
  INTEGER, PUBLIC :: mpi_edge_send_indexes_1d_counts_sum ! all number of mpi_send_indexes_1d_counts

  integer :: mpi_vertex_total_ordered_counts ! sum of mpi_vertex_indexes_all_ordered_counts
  integer, allocatable, dimension(:) :: mpi_vertex_indexes ! global edge index number for cells on each process
  integer, allocatable, dimension(:) :: mpi_vertex_indexes_all_ordered ! global edge index number for cells ordered by processes,only for process 0
  integer, allocatable, dimension(:) :: mpi_vertex_indexes_all_ordered_counts ! number of global edge index for cells on each process,only for process 0
  integer, allocatable, dimension(:) :: mpi_vertex_indexes_all_ordered_displs ! displacement of global edge index for cells on each process,only for process 0
!  ! number of global cell index for cells on each process, the last dimension is maxedges
!  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_counts_2d_maxedges 
!  ! displacement of global cell index for cells on each process. the last dimension is maxedges
!  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_displs_2d_maxedges 
!  ! number of global cell index for cells on each process,the last dimension is 3
!  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_counts_2d_3size
!  ! displacement of global cell index for cells on each process,the last dimension is 3
!  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_displs_2d_3size
!  ! number of global cell index for cells on each process,the last two dimension is 3 and maxedges
!  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_counts_3d_3size_maxedges
!  ! displacement of global cell index for cells on each process,the last two dimension is 3 and maxedges
!  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_displs_3d_3size_maxedges
!  ! number of global cell index for cells on each process,the last two dimension is 3 and 2
!  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_counts_3d_3size_2size
!  ! displacement of global cell index for cells on each process,the last two dimension is 3 and 2
!  integer, allocatable, dimension(:) :: mpi_cell_indexes_all_ordered_displs_3d_3size_2size


  integer, parameter :: MPI_MIDIUM_LEN = 256 ! medium length for message
  integer, parameter :: MPI_INPUT_UNIT = 88  ! input unit ID for reading input file
               

  ! variables used in reading netcdf input
  integer,allocatable,dimension(:) :: mpi_integer_buf_nCells
  integer,allocatable,dimension(:,:) :: mpi_integer_buf_nCells_2d
  real(kind=real_kind),allocatable,dimension(:) :: mpi_real_buf_nCells
  real(kind=real_kind),allocatable,dimension(:) :: mpi_real_buf_nCells_2
  real(kind=real_kind),allocatable,dimension(:) :: mpi_real_buf_nEdges
  integer,allocatable,dimension(:) :: mpi_integer_buf_nEdges
  integer,allocatable,dimension(:,:) :: mpi_integer_buf_nEdges_2d_2
  real(kind=real_kind),allocatable,dimension(:,:) :: mpi_real_buf_nEdges_2d_2
  integer,allocatable,dimension(:,:) :: mpi_integer_buf_nEdges_2d_max
  real(kind = real_kind),allocatable,dimension(:) :: mpi_real_buf_vertices
  real(kind = real_kind),allocatable,dimension(:,:) :: mpi_real_buf_vertices_2d
  integer,allocatable,dimension(:) :: mpi_integer_buf_vertices
  integer,allocatable,dimension(:,:) :: mpi_integer_buf_vertices_2d


  ! for mpi communication
  integer,allocatable,dimension(:) :: mpi_int_cell_sendbuf_1d
  integer,allocatable,dimension(:,:) :: mpi_int_cell_sendbuf_2d
  integer,allocatable,dimension(:) :: mpi_int_edge_sendbuf_1d
  integer,allocatable,dimension(:,:) :: mpi_int_edge_sendbuf_2d
  real(kind = real_kind),allocatable,dimension(:) :: mpi_real_cell_sendbuf_1d
  real(kind = real_kind),allocatable,dimension(:,:) :: mpi_real_cell_sendbuf_2d
  real(kind = real_kind),allocatable,dimension(:) :: mpi_real_edge_sendbuf_1d
  real(kind = real_kind),allocatable,dimension(:,:) :: mpi_real_edge_sendbuf_2d


end module mod_mpi_variables
