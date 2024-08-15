module history_mod

  use params_mod
  use const_mod
  use mesh_mod
  use time_mod
  use fiona_mod
  use log_mod
  use string
  use state_mod
  use forcing_mod
  use mask_mod
  use fre_dir_mod
  
  use wave_propag_mod
  use initial_mod
  use mod_mpi_interfaces
  use mod_mpi_variables
  use wam_source_module

  implicit none
  
  ! output
  integer :: nSta
  real(real_kind), allocatable :: SWH(:)         ! Significant Wave Height, meter
  real(real_kind), allocatable :: MWP(:)         ! Mean Wave Period, second
  real(real_kind), allocatable :: MWD(:)         ! Mean Wave Direction, degree
      
  ! station output
  integer, allocatable       :: staCellGlobal(:)

  integer, allocatable :: staCellLocal(:)
  real(real_kind), allocatable :: staFL(:,:,:)      ! Frequency Spectrum, (nSta,nDir,nFre)

  real(real_kind), allocatable :: staSWH(:)         ! Station's Significant Wave Height, meter (nSta)
  real(real_kind), allocatable :: staMWP(:)         ! Station's Mean Wave Period, second       (nSta)
  real(real_kind), allocatable :: staMWD(:)         ! Station's Mean Wave Direction, degree    (nSta)

  

  public history_init
  public history_write
  public history_create_new_file
  public calc_output
  public station_output_init
  public station_output_write

  interface history_write
    module procedure history_write_state
  end interface history_write

contains

  subroutine history_init()

    character(10) time_value, time_units
    real(real_kind) seconds, seconds2
    character(19) :: curr_date_format = 'YYYY-MM-DD_HH_MM_SS'
    
    allocate(SWH(nCells))
    allocate(MWP(nCells))
    allocate(MWD(nCells))
    SWH(:) = 0.0
    MWP(:) = 0.0
    MWD(:) = 0.0


    curr_date_format(1:4)   = curr_time_format(1:4)
    curr_date_format(6:7)   = curr_time_format(6:7)
    curr_date_format(9:10)  = curr_time_format(9:10)
    curr_date_format(12:13) = curr_time_format(12:13)
    curr_date_format(15:16) = curr_time_format(15:16)
    curr_date_format(18:19) = curr_time_format(18:19)
    
    if (station_output_option)  call station_output_init(curr_date_format)

    ! history_interval
    time_value = split_string(history_interval(1), ' ',  1)
    time_units = split_string(history_interval(1), ' ',  2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      if(mpi_rank==0) then
      call log_error('Invalid history interval ' // trim(history_interval(1)) // '!')
      end if
    end select

    ! frames_per_file
    time_value = split_string(frames_per_file, ' ',  1)
    time_units = split_string(frames_per_file, ' ',  2)
    read(time_value, *) seconds2
    select case (time_units)
    case ('days')
      seconds2 = seconds2 * 86400
    case ('hours')
      seconds2 = seconds2 * 3600
    case ('minutes')
      seconds2 = seconds2 * 60
    case ('seconds')
      seconds2 = seconds2
    case default
      if(mpi_rank==0) then
      call log_error('Invalid frames_per_file unit ' // trim(frames_per_file) // '!')
      end if
    end select

    call time_add_alert('history_write',    seconds=seconds )
    call time_add_alert('history_new_file', seconds=seconds2)

    if(mpi_rank==0) then

        if (output_file_prefix /= 'N/A') then
          call fiona_create_dataset('h0', desc=case_name, file_prefix=trim(adjustL(output_file_prefix))//'_'//curr_date_format)
        else
          call fiona_create_dataset('h0', desc=case_name, file_prefix=trim(adjustL(case_name))//'_'//curr_date_format)
        end if

        call fiona_add_att('h0', 'source',         'wave model')
        call fiona_add_att('h0', 'dt',             dt)
        call fiona_add_att('h0', 'author',         'NMEFC')
        call fiona_add_att('h0', 'on_a_sphere',    'YES')
        call fiona_add_att('h0', 'sphere_radius',  radius)
        ! Dimensions
        call fiona_add_dim('h0', 'Time',           add_var=.true.)
        !update to the global size, not the local size
        call fiona_add_dim('h0', 'nCells',         size=mpi_total_nCells)
        call fiona_add_dim('h0', 'nEdges',         size=mpi_total_nEdges)
        call fiona_add_dim('h0', 'nVertices',      size=mpi_total_nVertices)
        call fiona_add_dim('h0', 'TWO',            size=2)
        call fiona_add_dim('h0', 'vertexDegree',   size=vertexDegree)
        call fiona_add_dim('h0', 'maxEdges',       size=maxEdges)
        call fiona_add_dim('h0', 'maxEdges2',      size=maxEdges2)
        ! Mesh parameters
        call fiona_add_var('h0', 'lonCell',        long_name='Longitude on the cell',                       units='radian', dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'latCell',        long_name='Latitude on the cell',                        units='radian', dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'xCell',          long_name='Cartesian X on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'yCell',          long_name='Cartesian Y on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'zCell',          long_name='Cartesian Z on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'indexToCellID',  long_name='Global cell ID',                              units='1',      dim_names=['nCells      '],                 data_type='integer')
        call fiona_add_var('h0', 'lonEdge',        long_name='Longitude on the edge',                       units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'latEdge',        long_name='Latitude on the edge',                        units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'xEdge',          long_name='Cartesian X on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'yEdge',          long_name='Cartesian Y on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'zEdge',          long_name='Cartesian Z on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'indexToEdgeID',  long_name='Global edge ID',                              units='1',      dim_names=['nEdges      '],                 data_type='integer')
        call fiona_add_var('h0', 'lonVertex',      long_name='Longitude on the vertex',                     units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'latVertex',      long_name='Latitude on the vertex',                      units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'xVertex',        long_name='Cartesian X on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'yVertex',        long_name='Cartesian Y on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'zVertex',        long_name='Cartesian Z on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'indexToVertexID',long_name='Global vertex ID',                            units='1',      dim_names=['nVertices   '],                 data_type='integer')
        call fiona_add_var('h0', 'areaCell',       long_name='Primary cell area',                           units='m2',     dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'areaTriangle',   long_name='Dual cell area',                              units='m2',     dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'areaEdge',       long_name='Defined edge area',                           units='m2',     dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'nEdgesOnCell',   long_name='Edge number on the cell',                     units='1',      dim_names=['nCells      '],                 data_type='integer')
        call fiona_add_var('h0', 'nEdgesOnEdge',   long_name='Edge number to reconstruct tangent velocity', units='1',      dim_names=['nEdges      '],                 data_type='integer')
        call fiona_add_var('h0', 'cellsOnCell',    long_name='Cell indices that surround cell',             units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
        call fiona_add_var('h0', 'cellsOnEdge',    long_name='Cell indices that saddle cell',               units='1',      dim_names=['TWO         ', 'nEdges      '], data_type='integer')
!        call fiona_add_var('h0', 'cellsOnVertex',  long_name='Cell indices that surround vertex',           units='1',      dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
        call fiona_add_var('h0', 'edgesOnCell',    long_name='Edge indices on the cell',                    units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
!        call fiona_add_var('h0', 'edgesOnEdge',    long_name='Edge indices to reconstruct tangent velocity',units='1',      dim_names=['maxEdges2   ', 'nEdges      '], data_type='integer')
!        call fiona_add_var('h0', 'edgesOnVertex',  long_name='Edge indices on the vertex',                  units='1',      dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
        call fiona_add_var('h0', 'verticesOnCell', long_name='Vertex indices on the cell',                  units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
        call fiona_add_var('h0', 'verticesOnEdge', long_name='Vertex indices on the edge',                  units='1',      dim_names=['TWO         ', 'nEdges      '], data_type='integer')

        ! Dynamical variables
        call fiona_add_var('h0', 'bottomDepth',    long_name='Bottom depth on the cell',                    units='m',      dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'maskCell',       long_name='Mask, land=0, water=1',                       units='none',   dim_names=['nCells      '],                 data_type='integer')
        call fiona_add_var('h0', 'maskBdy',        long_name='Mask, land=0, water=1',                       units='none',   dim_names=['nCells      '],                 data_type='integer')
        call fiona_add_var('h0', 'u10Spd',         long_name='Wind speed at 10m on the cell',               units='m',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
        call fiona_add_var('h0', 'SWH',            long_name='Significant wave height',                     units='m',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
        call fiona_add_var('h0', 'MWP',            long_name='Mean Wave Period',                            units='s',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
        call fiona_add_var('h0', 'MWD',            long_name='Mean Wave Direction',                         units='deg',    dim_names=['nCells   ', 'Time     '],       data_type='real(8)')

        call log_notice('History module is initialized.')

    end if
    
    if(mpi_procs>1) then
      if (mpi_rank==0) then

        call fiona_open_dataset('mesh', file_path=mesh_file_path)

        call fiona_start_output('h0', time_elapsed_seconds(), new_file=.True.)
        
        allocate(mpi_real_buf_nCells(mpi_total_nCells))
        call fiona_input('mesh', 'lonCell',           mpi_real_buf_nCells)
        call fiona_output('h0', 'lonCell',    mpi_real_buf_nCells)

        call fiona_input('mesh', 'latCell',           mpi_real_buf_nCells)
        call fiona_output('h0', 'latCell',         mpi_real_buf_nCells)


        allocate(mpi_integer_buf_nCells(mpi_total_nCells))
        call fiona_input('mesh', 'indexToCellID',     mpi_integer_buf_nCells)
        call fiona_output('h0', 'indexToCellID',      mpi_integer_buf_nCells)

        call fiona_input('mesh', 'nEdgesOnCell',      mpi_integer_buf_nCells)
        call fiona_output('h0', 'nEdgesOnCell',    mpi_integer_buf_nCells)


        allocate(mpi_real_buf_nEdges(mpi_total_nEdges))
        call fiona_input('mesh', 'lonEdge',           mpi_real_buf_nEdges)
        call fiona_output('h0', 'lonEdge',         mpi_real_buf_nEdges)

        call fiona_input('mesh', 'latEdge',           mpi_real_buf_nEdges)
        call fiona_output('h0', 'latEdge',         mpi_real_buf_nEdges)

        allocate(mpi_integer_buf_nEdges(mpi_total_nEdges))
        call fiona_input('mesh', 'indexToEdgeID',     mpi_integer_buf_nEdges)
        call fiona_output('h0', 'indexToEdgeID',   mpi_integer_buf_nEdges)

        call fiona_input('mesh', 'nEdgesOnEdge',      mpi_integer_buf_nEdges)
        call fiona_output('h0', 'nEdgesOnEdge',    mpi_integer_buf_nEdges)

        allocate(mpi_real_buf_vertices(mpi_total_nVertices))
        call fiona_input('mesh', 'lonVertex',         mpi_real_buf_vertices)
        call fiona_output('h0', 'lonVertex',       mpi_real_buf_vertices)

        call fiona_input('mesh', 'latVertex',         mpi_real_buf_vertices)
        call fiona_output('h0', 'latVertex',       mpi_real_buf_vertices)

        !update value on vertex
        call fiona_input('mesh', 'xVertex',         mpi_real_buf_vertices)
        mpi_real_buf_vertices(1:mpi_total_nVertices)=mpi_real_buf_vertices(1:mpi_total_nVertices) * radius
        call fiona_output('h0', 'xVertex',       mpi_real_buf_vertices)

        call fiona_input('mesh', 'yVertex',         mpi_real_buf_vertices)
        mpi_real_buf_vertices(1:mpi_total_nVertices)=mpi_real_buf_vertices(1:mpi_total_nVertices) * radius
        call fiona_output('h0', 'yVertex',       mpi_real_buf_vertices)

        call fiona_input('mesh', 'zVertex',         mpi_real_buf_vertices)
        mpi_real_buf_vertices(1:mpi_total_nVertices)=mpi_real_buf_vertices(1:mpi_total_nVertices) * radius
        call fiona_output('h0', 'zVertex',       mpi_real_buf_vertices)

        call fiona_input('mesh', 'areaTriangle',         mpi_real_buf_vertices)
        mpi_real_buf_vertices(1:mpi_total_nVertices)=mpi_real_buf_vertices(1:mpi_total_nVertices) * radius**2
        call fiona_output('h0', 'areaTriangle',       mpi_real_buf_vertices)

        allocate(mpi_integer_buf_vertices(mpi_total_nVertices))
        call fiona_input('mesh', 'indexToVertexID',   mpi_integer_buf_vertices)
        call fiona_output('h0', 'indexToVertexID', mpi_integer_buf_vertices)

        allocate(mpi_integer_buf_nCells_2d(maxEdges,mpi_total_nCells))
        call fiona_input('mesh', 'cellsOnCell',       mpi_integer_buf_nCells_2d)
        call fiona_output('h0', 'cellsOnCell',     mpi_integer_buf_nCells_2d)  

        call fiona_input('mesh', 'edgesOnCell',       mpi_integer_buf_nCells_2d)
        call fiona_output('h0', 'edgesOnCell',     mpi_integer_buf_nCells_2d)

        call fiona_input('mesh', 'verticesOnCell',    mpi_integer_buf_nCells_2d)
        call fiona_output('h0', 'verticesOnCell',  mpi_integer_buf_nCells_2d)

        allocate(mpi_integer_buf_nEdges_2d_2(2,mpi_total_nEdges))
        call fiona_input('mesh', 'cellsOnEdge',       mpi_integer_buf_nEdges_2d_2)
        call fiona_output('h0', 'cellsOnEdge',     mpi_integer_buf_nEdges_2d_2)  

        call fiona_input('mesh', 'verticesOnEdge',    mpi_integer_buf_nEdges_2d_2)
        call fiona_output('h0', 'verticesOnEdge',  mpi_integer_buf_nEdges_2d_2)

!        allocate(mpi_integer_buf_vertices_2d(vertexDegree,mpi_total_nVertices))
!        call fiona_input('mesh', 'cellsOnVertex',     mpi_integer_buf_vertices_2d)
!        call fiona_output('h0', 'cellsOnVertex',   mpi_integer_buf_vertices_2d) 

!        call fiona_input('mesh', 'edgesOnVertex',     mpi_integer_buf_vertices_2d)
!        call fiona_output('h0', 'edgesOnVertex',   mpi_integer_buf_vertices_2d)

!        allocate(mpi_integer_buf_nEdges_2d_max(maxEdges2,mpi_total_nEdges))
!        call fiona_input('mesh', 'edgesOnEdge',       mpi_integer_buf_nEdges_2d_max)
!        call fiona_output('h0', 'edgesOnEdge',     mpi_integer_buf_nEdges_2d_max)

        deallocate(mpi_real_buf_nCells)
        deallocate(mpi_integer_buf_nCells)
        deallocate(mpi_real_buf_nEdges)
        deallocate(mpi_integer_buf_nEdges)
        deallocate(mpi_real_buf_vertices)
        deallocate(mpi_integer_buf_vertices)
        deallocate(mpi_integer_buf_nCells_2d)
        deallocate(mpi_integer_buf_nEdges_2d_2)
!        deallocate(mpi_integer_buf_vertices_2d)
!        deallocate(mpi_integer_buf_nEdges_2d_max)

      end if
      ! following variables are not changed after initialization
     ! call mpi_output_onedge('h0', 'areaEdge',areaEdge(1:nEdges))
     ! call mpi_output_onedge('h0', 'xEdge',xEdge(1:nEdges))
     ! call mpi_output_onedge('h0', 'yEdge',yEdge(1:nEdges))
     ! call mpi_output_onedge('h0', 'zEdge',zEdge(1:nEdges))

     ! call mpi_output_oncell('h0', 'xCell',xCell(1:nCells))
     ! call mpi_output_oncell('h0', 'yCell',yCell(1:nCells))
     ! call mpi_output_oncell('h0', 'zCell',zCell(1:nCells))
     ! call mpi_output_oncell('h0', 'areaCell',areaCell(1:nCells))
     ! call mpi_output_oncell('h0', 'maskCell',maskCell(1:nCells))
     ! call mpi_output_oncell('h0', 'maskBdy',maskBdy(1:nCells))

     ! call mpi_output_oncell('h0', 'bottomDepth',bottomDepth(1:nCells))
    else ! mpirun -np 1

      call fiona_start_output('h0', time_elapsed_seconds(), new_file=.True.)
      call fiona_output('h0', 'lonCell',    lonCell)
      call fiona_output('h0', 'latCell',    latCell)
      call fiona_output('h0', 'indexToCellID',      indexToCellID)
      call fiona_output('h0', 'nEdgesOnCell',    nEdgesOnCell)
      call fiona_output('h0', 'lonEdge',         lonEdge)
      call fiona_output('h0', 'latEdge',         latEdge)
      call fiona_output('h0', 'indexToEdgeID',   indexToEdgeID)
      call fiona_output('h0', 'nEdgesOnEdge',    nEdgesOnEdge)
      call fiona_output('h0', 'lonVertex',       lonVertex)
      call fiona_output('h0', 'latVertex',       latVertex)
      call fiona_output('h0', 'xVertex',       xVertex)
      call fiona_output('h0', 'yVertex',      yVertex)
      call fiona_output('h0', 'zVertex',       zVertex)
      call fiona_output('h0', 'areaTriangle',       areaTriangle)
      call fiona_output('h0', 'indexToVertexID', indexToVertexID)
      call fiona_output('h0', 'cellsOnCell',     cellsOnCell)   
      call fiona_output('h0', 'edgesOnCell',     edgesOncell)
      call fiona_output('h0', 'verticesOnCell',  verticesOnCell)
      call fiona_output('h0', 'cellsOnEdge',     cellsOnEdge)   
      call fiona_output('h0', 'verticesOnEdge',  verticesOnEdge)
!      call fiona_output('h0', 'cellsOnVertex',   cellsOnVertex) 
!      call fiona_output('h0', 'edgesOnVertex',   edgesOnVertex) 
!      call fiona_output('h0', 'edgesOnEdge',     edgesOnEdge)

      call fiona_output('h0', 'areaEdge',areaEdge(1:nEdges))
      call fiona_output('h0', 'xEdge',xEdge(1:nEdges))
      call fiona_output('h0', 'yEdge',yEdge(1:nEdges))
      call fiona_output('h0', 'zEdge',zEdge(1:nEdges))
      call fiona_output('h0', 'xCell',xCell(1:nCells))
      call fiona_output('h0', 'yCell',yCell(1:nCells))
      call fiona_output('h0', 'zCell',zCell(1:nCells))
      call fiona_output('h0', 'areaCell',areaCell(1:nCells))
      call fiona_output('h0', 'maskCell',maskCell(1:nCells))
      call fiona_output('h0', 'maskBdy',maskBdy(1:nCells))
    
      call fiona_output('h0', 'bottomDepth',bottomDepth(1:nCells))
    end if

    if (mpi_rank==0) then
      call fiona_end_output('h0')
    end if

  end subroutine history_init
  
  subroutine history_create_new_file()
    implicit none
    character(19) :: curr_date_format = 'YYYY-MM-DD_HH_MM_SS'

    curr_date_format(1:4)   = curr_time_format(1:4)
    curr_date_format(6:7)   = curr_time_format(6:7)
    curr_date_format(9:10)  = curr_time_format(9:10)
    curr_date_format(12:13) = curr_time_format(12:13)
    curr_date_format(15:16) = curr_time_format(15:16)
    curr_date_format(18:19) = curr_time_format(18:19)
    if(mpi_rank==0) then

        if (output_file_prefix /= 'N/A') then
          call fiona_create_dataset('h0', desc=case_name, file_prefix=trim(adjustL(output_file_prefix))//'_'//curr_date_format)
        else
          call fiona_create_dataset('h0', desc=case_name,file_prefix=trim(adjustL(case_name))//'_'//curr_date_format)
        end if

        call fiona_add_att('h0', 'source',         'wave model')
        call fiona_add_att('h0', 'dt',             dt)
        call fiona_add_att('h0', 'author',         'NMEFC')
        call fiona_add_att('h0', 'on_a_sphere',    'YES')
        call fiona_add_att('h0', 'sphere_radius',  radius)
        ! Dimensions
        call fiona_add_dim('h0', 'Time',           add_var=.true.)
        !update to the global size, not the local size
        call fiona_add_dim('h0', 'nCells',         size=mpi_total_nCells)
        call fiona_add_dim('h0', 'nEdges',         size=mpi_total_nEdges)
        call fiona_add_dim('h0', 'nVertices',      size=mpi_total_nVertices)
        call fiona_add_dim('h0', 'TWO',            size=2)
        call fiona_add_dim('h0', 'vertexDegree',   size=vertexDegree)
        call fiona_add_dim('h0', 'maxEdges',       size=maxEdges)
        call fiona_add_dim('h0', 'maxEdges2',      size=maxEdges2)
        ! Mesh parameters
        call fiona_add_var('h0', 'lonCell',        long_name='Longitude on the cell',                       units='radian', dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'latCell',        long_name='Latitude on the cell',                        units='radian', dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'xCell',          long_name='Cartesian X on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'yCell',          long_name='Cartesian Y on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'zCell',          long_name='Cartesian Z on the cell',                     units='m',      dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'indexToCellID',  long_name='Global cell ID',                              units='1',      dim_names=['nCells      '],                 data_type='integer')
        call fiona_add_var('h0', 'lonEdge',        long_name='Longitude on the edge',                       units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'latEdge',        long_name='Latitude on the edge',                        units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'xEdge',          long_name='Cartesian X on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'yEdge',          long_name='Cartesian Y on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'zEdge',          long_name='Cartesian Z on the edge',                     units='m',      dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'indexToEdgeID',  long_name='Global edge ID',                              units='1',      dim_names=['nEdges      '],                 data_type='integer')
        call fiona_add_var('h0', 'lonVertex',      long_name='Longitude on the vertex',                     units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'latVertex',      long_name='Latitude on the vertex',                      units='radian', dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'xVertex',        long_name='Cartesian X on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'yVertex',        long_name='Cartesian Y on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'zVertex',        long_name='Cartesian Z on the vertex',                   units='m',      dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'indexToVertexID',long_name='Global vertex ID',                            units='1',      dim_names=['nVertices   '],                 data_type='integer')
        call fiona_add_var('h0', 'areaCell',       long_name='Primary cell area',                           units='m2',     dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'areaTriangle',   long_name='Dual cell area',                              units='m2',     dim_names=['nVertices   '],                 data_type='real(8)')
        call fiona_add_var('h0', 'areaEdge',       long_name='Defined edge area',                           units='m2',     dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'nEdgesOnCell',   long_name='Edge number on the cell',                     units='1',      dim_names=['nCells      '],                 data_type='integer')
        call fiona_add_var('h0', 'nEdgesOnEdge',   long_name='Edge number to reconstruct tangent velocity', units='1',      dim_names=['nEdges      '],                 data_type='integer')
        call fiona_add_var('h0', 'cellsOnCell',    long_name='Cell indices that surround cell',             units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
        call fiona_add_var('h0', 'cellsOnEdge',    long_name='Cell indices that saddle cell',               units='1',      dim_names=['TWO         ', 'nEdges      '], data_type='integer')
!        call fiona_add_var('h0', 'cellsOnVertex',  long_name='Cell indices that surround vertex',           units='1',      dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
        call fiona_add_var('h0', 'edgesOnCell',    long_name='Edge indices on the cell',                    units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
!        call fiona_add_var('h0', 'edgesOnEdge',    long_name='Edge indices to reconstruct tangent velocity',units='1',      dim_names=['maxEdges2   ', 'nEdges      '], data_type='integer')
!        call fiona_add_var('h0', 'edgesOnVertex',  long_name='Edge indices on the vertex',                  units='1',      dim_names=['vertexDegree', 'nVertices   '], data_type='integer')
        call fiona_add_var('h0', 'verticesOnCell', long_name='Vertex indices on the cell',                  units='1',      dim_names=['maxEdges    ', 'nCells      '], data_type='integer')
        call fiona_add_var('h0', 'verticesOnEdge', long_name='Vertex indices on the edge',                  units='1',      dim_names=['TWO         ', 'nEdges      '], data_type='integer')

        ! Dynamical variables
        call fiona_add_var('h0', 'bottomDepth',    long_name='Bottom depth on the cell',                    units='m',      dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('h0', 'maskCell',       long_name='Mask, land=0, water=1',                       units='none',   dim_names=['nCells      '],                 data_type='integer')
        call fiona_add_var('h0', 'maskBdy',        long_name='Mask, land=0, water=1',                       units='none',   dim_names=['nCells      '],                 data_type='integer')
        call fiona_add_var('h0', 'N',              long_name='Wave action on the cell',                     units='m',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
        call fiona_add_var('h0', 'u10Spd',         long_name='Wind speed at 10m on the cell',               units='m',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
        call fiona_add_var('h0', 'SWH',            long_name='Significant wave height',                     units='m',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
        call fiona_add_var('h0', 'MWP',            long_name='Mean Wave Period',                            units='s',      dim_names=['nCells   ', 'Time     '],       data_type='real(8)')
        call fiona_add_var('h0', 'MWD',            long_name='Mean Wave Direction',                         units='deg',    dim_names=['nCells   ', 'Time     '],       data_type='real(8)')

        call log_notice('History new file!')

    end if
    
    if(mpi_procs>1) then
      if (mpi_rank==0) then

        call fiona_start_output('h0', time_elapsed_seconds(), new_file=.True.)
        call fiona_open_dataset('mesh', file_path=mesh_file_path)
        
        allocate(mpi_real_buf_nCells(mpi_total_nCells))
        call fiona_input('mesh', 'lonCell',           mpi_real_buf_nCells)
        call fiona_output('h0', 'lonCell',    mpi_real_buf_nCells)

        call fiona_input('mesh', 'latCell',           mpi_real_buf_nCells)
        call fiona_output('h0', 'latCell',         mpi_real_buf_nCells)


        allocate(mpi_integer_buf_nCells(mpi_total_nCells))
        call fiona_input('mesh', 'indexToCellID',     mpi_integer_buf_nCells)
        call fiona_output('h0', 'indexToCellID',      mpi_integer_buf_nCells)

        call fiona_input('mesh', 'nEdgesOnCell',      mpi_integer_buf_nCells)
        call fiona_output('h0', 'nEdgesOnCell',    mpi_integer_buf_nCells)


        allocate(mpi_real_buf_nEdges(mpi_total_nEdges))
        call fiona_input('mesh', 'lonEdge',           mpi_real_buf_nEdges)
        call fiona_output('h0', 'lonEdge',         mpi_real_buf_nEdges)

        call fiona_input('mesh', 'latEdge',           mpi_real_buf_nEdges)
        call fiona_output('h0', 'latEdge',         mpi_real_buf_nEdges)

        allocate(mpi_integer_buf_nEdges(mpi_total_nEdges))
        call fiona_input('mesh', 'indexToEdgeID',     mpi_integer_buf_nEdges)
        call fiona_output('h0', 'indexToEdgeID',   mpi_integer_buf_nEdges)

        call fiona_input('mesh', 'nEdgesOnEdge',      mpi_integer_buf_nEdges)
        call fiona_output('h0', 'nEdgesOnEdge',    mpi_integer_buf_nEdges)

        allocate(mpi_real_buf_vertices(mpi_total_nVertices))
        call fiona_input('mesh', 'lonVertex',         mpi_real_buf_vertices)
        call fiona_output('h0', 'lonVertex',       mpi_real_buf_vertices)

        call fiona_input('mesh', 'latVertex',         mpi_real_buf_vertices)
        call fiona_output('h0', 'latVertex',       mpi_real_buf_vertices)

        !update value on vertex
        call fiona_input('mesh', 'xVertex',         mpi_real_buf_vertices)
        mpi_real_buf_vertices(1:mpi_total_nVertices)=mpi_real_buf_vertices(1:mpi_total_nVertices) * radius
        call fiona_output('h0', 'xVertex',       mpi_real_buf_vertices)

        call fiona_input('mesh', 'yVertex',         mpi_real_buf_vertices)
        mpi_real_buf_vertices(1:mpi_total_nVertices)=mpi_real_buf_vertices(1:mpi_total_nVertices) * radius
        call fiona_output('h0', 'yVertex',       mpi_real_buf_vertices)

        call fiona_input('mesh', 'zVertex',         mpi_real_buf_vertices)
        mpi_real_buf_vertices(1:mpi_total_nVertices)=mpi_real_buf_vertices(1:mpi_total_nVertices) * radius
        call fiona_output('h0', 'zVertex',       mpi_real_buf_vertices)

        call fiona_input('mesh', 'areaTriangle',         mpi_real_buf_vertices)
        mpi_real_buf_vertices(1:mpi_total_nVertices)=mpi_real_buf_vertices(1:mpi_total_nVertices) * radius**2
        call fiona_output('h0', 'areaTriangle',       mpi_real_buf_vertices)

        allocate(mpi_integer_buf_vertices(mpi_total_nVertices))
        call fiona_input('mesh', 'indexToVertexID',   mpi_integer_buf_vertices)
        call fiona_output('h0', 'indexToVertexID', mpi_integer_buf_vertices)

        allocate(mpi_integer_buf_nCells_2d(maxEdges,mpi_total_nCells))
        call fiona_input('mesh', 'cellsOnCell',       mpi_integer_buf_nCells_2d)
        call fiona_output('h0', 'cellsOnCell',     mpi_integer_buf_nCells_2d)   

        call fiona_input('mesh', 'edgesOnCell',       mpi_integer_buf_nCells_2d)
        call fiona_output('h0', 'edgesOnCell',     mpi_integer_buf_nCells_2d)

        call fiona_input('mesh', 'verticesOnCell',    mpi_integer_buf_nCells_2d)
        call fiona_output('h0', 'verticesOnCell',  mpi_integer_buf_nCells_2d)

        allocate(mpi_integer_buf_nEdges_2d_2(2,mpi_total_nEdges))
        call fiona_input('mesh', 'cellsOnEdge',       mpi_integer_buf_nEdges_2d_2)
        call fiona_output('h0', 'cellsOnEdge',     mpi_integer_buf_nEdges_2d_2)   

        call fiona_input('mesh', 'verticesOnEdge',    mpi_integer_buf_nEdges_2d_2)
        call fiona_output('h0', 'verticesOnEdge',  mpi_integer_buf_nEdges_2d_2)

!        allocate(mpi_integer_buf_vertices_2d(vertexDegree,mpi_total_nVertices))
!        call fiona_input('mesh', 'cellsOnVertex',     mpi_integer_buf_vertices_2d)
!        call fiona_output('h0', 'cellsOnVertex',   mpi_integer_buf_vertices_2d)

!        call fiona_input('mesh', 'edgesOnVertex',     mpi_integer_buf_vertices_2d)
!        call fiona_output('h0', 'edgesOnVertex',   mpi_integer_buf_vertices_2d)  

!        allocate(mpi_integer_buf_nEdges_2d_max(maxEdges2,mpi_total_nEdges))
!        call fiona_input('mesh', 'edgesOnEdge',       mpi_integer_buf_nEdges_2d_max)
!        call fiona_output('h0', 'edgesOnEdge',     mpi_integer_buf_nEdges_2d_max)

        deallocate(mpi_real_buf_nCells)
        deallocate(mpi_integer_buf_nCells)
        deallocate(mpi_real_buf_nEdges)
        deallocate(mpi_integer_buf_nEdges)
        deallocate(mpi_real_buf_vertices)
        deallocate(mpi_integer_buf_vertices)
        deallocate(mpi_integer_buf_nCells_2d)
        deallocate(mpi_integer_buf_nEdges_2d_2)
!        deallocate(mpi_integer_buf_vertices_2d)
!        deallocate(mpi_integer_buf_nEdges_2d_max)

      end if
      ! following variables are not changed after initialization
      call mpi_output_onedge('h0', 'areaEdge',areaEdge(1:nEdges))
      call mpi_output_onedge('h0', 'xEdge',xEdge(1:nEdges))
      call mpi_output_onedge('h0', 'yEdge',yEdge(1:nEdges))
      call mpi_output_onedge('h0', 'zEdge',zEdge(1:nEdges))

      call mpi_output_oncell('h0', 'xCell',xCell(1:nCells))
      call mpi_output_oncell('h0', 'yCell',yCell(1:nCells))
      call mpi_output_oncell('h0', 'zCell',zCell(1:nCells))
      call mpi_output_oncell('h0', 'areaCell',areaCell(1:nCells))
      call mpi_output_oncell('h0', 'maskCell',maskCell(1:nCells))
      call mpi_output_oncell('h0', 'maskBdy',maskBdy(1:nCells))

      call mpi_output_oncell('h0', 'bottomDepth',bottomDepth(1:nCells))
    else

      call fiona_start_output('h0', time_elapsed_seconds(), new_file=.True.)
      call fiona_output('h0', 'lonCell',    lonCell)
      call fiona_output('h0', 'latCell',    latCell)
      call fiona_output('h0', 'indexToCellID',      indexToCellID)
      call fiona_output('h0', 'nEdgesOnCell',    nEdgesOnCell)
      call fiona_output('h0', 'lonEdge',         lonEdge)
      call fiona_output('h0', 'latEdge',         latEdge)
      call fiona_output('h0', 'indexToEdgeID',   indexToEdgeID)
      call fiona_output('h0', 'nEdgesOnEdge',    nEdgesOnEdge)
      call fiona_output('h0', 'lonVertex',       lonVertex)
      call fiona_output('h0', 'latVertex',       latVertex)
      call fiona_output('h0', 'xVertex',       xVertex)
      call fiona_output('h0', 'yVertex',      yVertex)
      call fiona_output('h0', 'zVertex',       zVertex)
      call fiona_output('h0', 'areaTriangle',       areaTriangle)
      call fiona_output('h0', 'indexToVertexID', indexToVertexID)
      call fiona_output('h0', 'cellsOnCell',     cellsOnCell)   
      call fiona_output('h0', 'edgesOnCell',     edgesOncell)
      call fiona_output('h0', 'verticesOnCell',  verticesOnCell)
      call fiona_output('h0', 'cellsOnEdge',     cellsOnEdge)   
      call fiona_output('h0', 'verticesOnEdge',  verticesOnEdge)
!      call fiona_output('h0', 'cellsOnVertex',   cellsOnVertex) 
!      call fiona_output('h0', 'edgesOnVertex',   edgesOnVertex) 
!      call fiona_output('h0', 'edgesOnEdge',     edgesOnEdge)

      call fiona_output('h0', 'areaEdge',areaEdge(1:nEdges))
      call fiona_output('h0', 'xEdge',xEdge(1:nEdges))
      call fiona_output('h0', 'yEdge',yEdge(1:nEdges))
      call fiona_output('h0', 'zEdge',zEdge(1:nEdges))
      call fiona_output('h0', 'xCell',xCell(1:nCells))
      call fiona_output('h0', 'yCell',yCell(1:nCells))
      call fiona_output('h0', 'zCell',zCell(1:nCells))
      call fiona_output('h0', 'areaCell',areaCell(1:nCells))
      call fiona_output('h0', 'maskCell',maskCell(1:nCells))
      call fiona_output('h0', 'maskBdy',maskBdy(1:nCells))
    
      call fiona_output('h0', 'bottomDepth',bottomDepth(1:nCells))
    end if

    if (mpi_rank==0) then
      call fiona_end_output('h0')
    end if

  end subroutine history_create_new_file

  subroutine history_final()

    call log_notice('History module is finalized.')

  end subroutine history_final

  subroutine history_write_state()

!    if(mpi_procs>1) then
!      if (mpi_rank==0) then
!         call fiona_start_output('h0', time_elapsed_seconds(), new_file=.False.)
!      end if
!      call mpi_output_oncell('h0', 'u10Spd',          forcing%u10Spd(1:nCells))
!      call mpi_output_oncell('h0', 'SWH',             SWH(1:nCells))
!      call mpi_output_oncell('h0', 'MWP',             MWP(1:nCells))
!      call mpi_output_oncell('h0', 'MWD',             MWD(1:nCells))
!    else
!      call fiona_start_output('h0', time_elapsed_seconds(), new_file=.False.)
!      call fiona_output('h0', 'u10Spd',          forcing%u10Spd(1:nCells))
!      call fiona_output('h0', 'SWH',             SWH(1:nCells))
!      call fiona_output('h0', 'MWP',             MWP(1:nCells))
!      call fiona_output('h0', 'MWD',             MWD(1:nCells))
!    end if
!
!    if (mpi_rank==0) then
!      call fiona_end_output('h0')
!    end if

  end subroutine history_write_state
  
  subroutine calc_output(N, SWH, MWP, MWD)
     real(real_kind), intent(in) :: N(0:,:,:)
     real(real_kind), intent(out):: SWH(:) ! Significant Wave Height, meter 
     real(real_kind), intent(out):: MWP(:) ! Mean Wave Period, second 
     real(real_kind), intent(out):: MWD(:) ! Mean Wave Direction, degree

     integer:: iCell, k, m
     
     real(real_kind) :: TEMP(nCells, nFre) !
     real(real_kind) :: EMEAN(nCells)      ! Total energy
     real(real_kind) :: Etail(nCells)      ! Tail energy
     
     real(real_kind) :: temps1, temps2, SI, CI ! for calculate wave direction
      
     !----------- Total energe and SWH 
     do m = 1, nFre
        do iCell = 1, nCells
           TEMP(iCell, m) = 0.0
           do k = 1, nDir
              TEMP(iCell, m) = TEMP(iCell, m) + N(iCell,k,m)
           end do
        end do 
     end do 

     do iCell = 1, nCells
        EMEAN(iCell) = MO_TAIL*TEMP(iCell,nFre)
        do m = 1, nFre
           EMEAN(iCell) = EMEAN(iCell) + TEMP(iCell,m)*delFre(m)*delTheta
        end do
        EMEAN(iCell) = max(EMEAN(iCell), EMIN) ! Total energy
 
        SWH(iCell) = 4.0*sqrt(EMEAN(iCell))  ! Significant Wave Height

     end do

     !----------- MWP, Mean Wave Period
     ! Tail Energy
     do iCell = 1, nCells
         Etail(iCell) = MM1_TAIL*TEMP(iCell, nFre)
     end do
     
     ! INTEGRATE OVER FREQUENCIES
     do iCell = 1, nCells
        do m = 1, nFre
           Etail(iCell) = Etail(iCell) + TEMP(iCell,m)*DFIMOFR(m)
        end do
     end do

     do iCell = 1, nCells
        MWP(iCell) = EMEAN(iCell)/MAX(Etail(iCell), EMIN) !MEAN FREQUENCY 
        MWP(iCell) = 1.0_rk/MWP(iCell) 
     end do
     

     !----------- MWD, Mean Wave Direction
     do iCell = 1, nCells
        SI = 0.0_rk
        CI = 0.0_rk
        temps2 = 0.0_rk
        do k = 1, nDir
           temps1 = 0.0_rk
           do m = 1, nFre
              temps1 = temps1 + N(iCell,k,m)*DFIM(m)
           end do
           SI = SI + sinTheta(k)*temps1
           CI = CI + cosTheta(k)*temps1
        end do
       
        if (CI .eq. 0.0_rk) CI = 0.1E-30
        MWD(iCell) = ATAN2(SI,CI)
        if (MWD(iCell) .lt. 0.0) MWD(iCell) = MWD(iCell) + pi2
        if (MWD(iCell) .gt. pi2-0.001) MWD(iCell) = 0.0
        MWD(iCell) = MWD(iCell)*deg
        
     end do
  end subroutine calc_output

  subroutine station_output_init(curr_date_format)
      character(19), intent(in) :: curr_date_format 

      integer :: i, iCellGlobal
      integer :: funit, stat
      
      integer, allocatable, dimension(:) :: staCellLocal_tmp
      integer, allocatable, dimension(:) :: mpi_nSta_index_tmp
      real(real_kind), allocatable :: latSta(:)
      real(real_kind), allocatable :: lonSta(:)
      real(real_kind), allocatable :: mpi_buf_nSta_sum(:)
      real(real_kind), allocatable :: mpi_adjust_buf_nSta_sum(:)
      !staName save all name of stations
      character(max_name_len), allocatable :: staName(:)

      ! open the file
      if(mpi_rank==0) then
        funit = 199
        open(funit, file=trim(adjustL(station_info_path)),status='old',iostat=stat)
        if (stat .ne. 0) then
          print*, "Error: "//station_info_path//' cannot be opened!'
          stop
        end if
        
        read(funit, *, iostat=stat) mpi_nSta_sum  
        if (stat .ne. 0) then
          print*, "Error: "//station_info_path//' reading error!'
          stop
        end if
        
        allocate(staCellGlobal(mpi_nSta_sum))
        allocate(staName(mpi_nSta_sum))

        do i =1 , mpi_nSta_sum
           read(funit, *, iostat=stat) staCellGlobal(i),  staName(i) 
           if (stat .ne. 0) then
             print*, "Error: "//station_info_path//' reading error!'
             stop
           end if
        end do
        close(funit) 

        if(mpi_procs>1) then
          call MPI_BCAST(mpi_nSta_sum, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, mpi_err)
          call MPI_BCAST(staCellGlobal, mpi_nSta_sum, MPI_INTEGER, 0,MPI_COMM_WORLD, mpi_err)
          call MPI_BCAST(staName, max_name_len*mpi_nSta_sum, MPI_CHARACTER, 0,MPI_COMM_WORLD, mpi_err)
        end if
      else

        if(mpi_procs>1) then
          call MPI_BCAST(mpi_nSta_sum, 1, MPI_INTEGER, 0,MPI_COMM_WORLD, mpi_err)
          allocate(staCellGlobal(mpi_nSta_sum))
          call MPI_BCAST(staCellGlobal, mpi_nSta_sum, MPI_INTEGER, 0,MPI_COMM_WORLD, mpi_err)
          allocate(staName(mpi_nSta_sum))
          call MPI_BCAST(staName, max_name_len*mpi_nSta_sum, MPI_CHARACTER, 0,MPI_COMM_WORLD, mpi_err)
        end if
      end if

      nSta=0
      allocate(staCellLocal_tmp(mpi_nSta_sum))
      allocate(mpi_nSta_index_tmp(mpi_nSta_sum))
      do i =  1, mpi_nSta_sum
         iCellGlobal = staCellGlobal(i)
         if (mpi_procs>1) then
            if ( mpi_rank == mpi_cell_indexes_all(iCellGlobal)) then
               nSta=nSta+1
               staCellLocal_tmp(nSta) = mpi_findloc(mpi_cell_indexes, iCellGlobal)
               mpi_nSta_index_tmp(nSta) = i
            end if
         else
            nSta=nSta+1
            staCellLocal_tmp(nSta) = iCellGlobal
         end if
      end do

      !there is possibility that nSta in some processes is 0
      if(nSta>0) then
        allocate(latSta(nSta))
        allocate(lonSta(nSta))
        allocate(staCellLocal(nSta))
        allocate(staSWH(nSta))
        allocate(staMWP(nSta))
        allocate(staMWD(nSta))
        allocate(staFL(nSta,nDir,nFre))

        do i=1,nSta
          staCellLocal(i)=staCellLocal_tmp(i)
          latSta(i) = latCell(staCellLocal(i))
          lonSta(i) = lonCell(staCellLocal(i))
        end do
        staSWH = 0.0
        staMWP = 0.0
        staMWD = 0.0
        staFL  = 0.0
      end if

      deallocate(staCellLocal_tmp)
      
      if(mpi_procs>1) then
        if(mpi_rank==0) then
          allocate(mpi_nSta_all(mpi_procs))
          call MPI_GATHER(nSta, 1, MPI_INTEGER, mpi_nSta_all, 1,MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
           
          allocate(mpi_nSta_all_displs(mpi_procs))
          mpi_nSta_all_displs(1)=0                                     
          do i = 2, mpi_procs                                             
            mpi_nSta_all_displs(i) = mpi_nSta_all_displs(i - 1) + mpi_nSta_all(i - 1) 
          end do

          allocate(mpi_nSta_index_all(mpi_nSta_sum))
          call MPI_GATHERV(mpi_nSta_index_tmp, nSta, MPI_INTEGER, mpi_nSta_index_all,mpi_nSta_all, &
                            mpi_nSta_all_displs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
        else
          call MPI_GATHER(nSta, 1, MPI_INTEGER, mpi_nSta_all, 1,MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
          call MPI_GATHERV(mpi_nSta_index_tmp, nSta, MPI_INTEGER, mpi_nSta_index_all,mpi_nSta_all, &
                            mpi_nSta_all_displs, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
        end if
      end if

      deallocate(mpi_nSta_index_tmp)

      
      if (mpi_rank==0) then
          call fiona_create_dataset('s0', desc=case_name, file_prefix=trim(adjustL(output_file_prefix))//'_station_'//curr_date_format)

          call fiona_add_att('s0', 'source',    'wave model')
          call fiona_add_att('s0', 'author',    'NMEFC')
          ! Dimensions
          call fiona_add_dim('s0', 'Time',      add_var=.true.)
          !nSta is local station size,mpi_nSta_sum is all station size
          call fiona_add_dim('s0', 'nSta',      size=mpi_nSta_sum)
          call fiona_add_dim('s0', 'nFre',      size=nFre)
          call fiona_add_dim('s0', 'nDir',      size=nDir)
          call fiona_add_dim('s0', 'nStr',      size=max_name_len)
          ! Mesh parameters
          call fiona_add_var('s0', 'staCell',   long_name='Station cell index',                      units='degree', dim_names=['nSta'],          data_type='integer')
          call fiona_add_var('s0', 'staName',   long_name='Station name',                            units='None',   dim_names=['nStr', 'nSta'],  data_type='character')
          call fiona_add_var('s0', 'fre',       long_name='Frequency',                               units='Hz',     dim_names=['nFre'],          data_type='real(8)')
          call fiona_add_var('s0', 'sigma',     long_name='Frequency',                               units='radian', dim_names=['nFre'],          data_type='real(8)')
          call fiona_add_var('s0', 'theta',     long_name='Direction, True north is 0, clockwise',   units='degree', dim_names=['nDir'],          data_type='real(8)')

          call fiona_add_var('s0', 'lonSta',    long_name='Longitude of the stations',               units='degree', dim_names=['nSta'],          data_type='real(8)')
          call fiona_add_var('s0', 'latSta',    long_name='Latitude of the stations',                units='degree', dim_names=['nSta'],          data_type='real(8)')
          call fiona_add_var('s0', 'staSWH',    long_name='Significant wave height',                 units='m',      dim_names=['nSta', 'Time'],  data_type='real(8)')
          call fiona_add_var('s0', 'staMWP',    long_name='Mean Wave Period',                        units='s',      dim_names=['nSta', 'Time'],  data_type='real(8)')
          call fiona_add_var('s0', 'staMWD',    long_name='Mean Wave Direction',                     units='deg',    dim_names=['nSta', 'Time'],  data_type='real(8)')
          if (spectrum_output_option) then
             call fiona_add_var('s0', 'staFL',     long_name='Frequency Spectrum',                      units='m^2 s',  dim_names=['nSta', 'nDir', 'nFre', 'Time'],  data_type='real(8)')
          end if
      end if

      if (mpi_procs>1) then
          ! mutiple process output
          if (mpi_rank==0) then
              ! output
              call fiona_start_output('s0', time_elapsed_seconds(), new_file=.True.)
              ! no need communication
              call fiona_output('s0', 'fre',       fre)
              call fiona_output('s0', 'sigma',     sigma)
              call fiona_output('s0', 'theta',     theta*deg)

              call fiona_output('s0', 'staCell',   staCellGlobal)
              call fiona_output('s0', 'staName',   staName)
              
              ! need communication
              allocate(mpi_buf_nSta_sum(mpi_nSta_sum))
              allocate(mpi_adjust_buf_nSta_sum(mpi_nSta_sum))
              call MPI_GATHERV(lonSta, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                               mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
              !adjust array to the order same with input station file
              call mpi_adjust_sta(mpi_buf_nSta_sum,mpi_adjust_buf_nSta_sum)
              call fiona_output('s0', 'lonSta',    mpi_adjust_buf_nSta_sum*deg)
              call MPI_GATHERV(latSta, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                               mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
              call mpi_adjust_sta(mpi_buf_nSta_sum,mpi_adjust_buf_nSta_sum)
              call fiona_output('s0', 'latSta',    mpi_adjust_buf_nSta_sum*deg)
              
              call fiona_end_output('s0')
              deallocate(mpi_buf_nSta_sum)
              deallocate(mpi_adjust_buf_nSta_sum)
          else
              call MPI_GATHERV(lonSta, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                               mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
              call MPI_GATHERV(latSta, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                               mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
          end if

      else
          ! single process output
          call fiona_start_output('s0', time_elapsed_seconds(), new_file=.True.)
          call fiona_output('s0', 'fre',       fre)
          call fiona_output('s0', 'sigma',     sigma)
          call fiona_output('s0', 'theta',     theta*deg)

          call fiona_output('s0', 'staCell',   staCellGlobal)
          call fiona_output('s0', 'staName',   staName)
          call fiona_output('s0', 'lonSta',    lonSta*deg)
          call fiona_output('s0', 'latSta',    latSta*deg)
          if(spectrum_output_option) call fiona_output('s0', 'staFL',     staFL)
          
          call fiona_end_output('s0')
      end if

  end subroutine station_output_init

  subroutine station_output_write(N)
     real(real_kind), intent(in) :: N(0:,:,:)
     real(real_kind), allocatable :: mpi_buf_nSta_sum(:)
     real(real_kind), allocatable :: mpi_adjust_buf_nSta_sum(:)

     integer :: i, iCell, k, m

    if (spectrum_output_option) then
        do i = 1, size(staCellLocal)
            iCell = staCellLocal(i)
             do k = 1, nDir
                do m = 1, nFre
                   staFL(i,k,m) = N(iCell,k,m)
                end do
             end do
        end do
     end if
 
     do i =  1, size(staCellLocal)
         iCell = staCellLocal(i)
           staSWH(i) = SWH(iCell)
           staMWP(i) = MWP(iCell)
           staMWD(i) = MWD(iCell)
     end do
     
     if (mpi_procs>1) then
         ! mutiple process output
         if(mpi_rank==0) then
             call fiona_start_output('s0', time_elapsed_seconds(), new_file=.False.)
             allocate(mpi_buf_nSta_sum(mpi_nSta_sum))
             allocate(mpi_adjust_buf_nSta_sum(mpi_nSta_sum))
             call MPI_GATHERV(staSWH, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                              mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
             !adjust array to the order same with input station file
             call mpi_adjust_sta(mpi_buf_nSta_sum,mpi_adjust_buf_nSta_sum)
             call fiona_output('s0', 'staSWH', mpi_adjust_buf_nSta_sum)

             call MPI_GATHERV(staMWP, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                              mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
             call mpi_adjust_sta(mpi_buf_nSta_sum,mpi_adjust_buf_nSta_sum)
             call fiona_output('s0', 'staMWP', mpi_adjust_buf_nSta_sum)

             call MPI_GATHERV(staMWD, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                              mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
             call mpi_adjust_sta(mpi_buf_nSta_sum,mpi_adjust_buf_nSta_sum)
             call fiona_output('s0', 'staMWD', mpi_adjust_buf_nSta_sum)
             
             call fiona_end_output('s0')
             deallocate(mpi_buf_nSta_sum)
             deallocate(mpi_adjust_buf_nSta_sum)
         else
             call MPI_GATHERV(staSWH, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                              mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
             call MPI_GATHERV(staMWP, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                              mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)
             call MPI_GATHERV(staMWD, nSta, mpi_real_kind, mpi_buf_nSta_sum,mpi_nSta_all, &
                              mpi_nSta_all_displs, mpi_real_kind, 0, MPI_COMM_WORLD, mpi_err)

         end if
     else
         ! single process output
         call fiona_start_output('s0', time_elapsed_seconds(), new_file=.False.)
         call fiona_output('s0', 'staSWH', staSWH)
         call fiona_output('s0', 'staMWP', staMWP)
         call fiona_output('s0', 'staMWD', staMWD)
         if (spectrum_output_option) call fiona_output('s0', 'staFL' , staFL )
     end if

    if (mpi_rank==0) then
      call fiona_end_output('s0')
    end if

  end subroutine station_output_write

end module history_mod
