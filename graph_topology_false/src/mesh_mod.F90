module mesh_mod

  use const_mod
  use params_mod
  use log_mod
  use fiona_mod
  use mod_mpi_interfaces

  implicit none

  ! Dimension sizes
  integer nCells
  integer nEdges
  integer nVertices
  integer vertexDegree                                    ! The maximum number of cells connected to a dual cell (or the number of corners in a dual cell)
  integer maxEdges
  integer maxEdges2
  integer, allocatable :: nEdgesOnCell(:)                 ! Number of edges on a given cell
  integer, allocatable :: nEdgesOnEdge(:)                 ! Number of edges on a given edge to reconstruct tangential velocities
  integer, allocatable :: nCellsOnVertex(:)               ! Number of cells connected with a given vertex
  ! Coordinates

  real(real_kind), allocatable :: latCell(:)
  real(real_kind), allocatable :: lonCell(:)
  real(real_kind), allocatable :: xCell(:)
  real(real_kind), allocatable :: yCell(:)
  real(real_kind), allocatable :: zCell(:)
  real(real_kind), allocatable :: latEdge(:)
  real(real_kind), allocatable :: lonEdge(:)
  real(real_kind), allocatable :: xEdge(:)
  real(real_kind), allocatable :: yEdge(:)
  real(real_kind), allocatable :: zEdge(:)
  real(real_kind), allocatable :: latVertex(:)
  real(real_kind), allocatable :: lonVertex(:)
  real(real_kind), allocatable :: xVertex(:)
  real(real_kind), allocatable :: yVertex(:)
  real(real_kind), allocatable :: zVertex(:)
  ! Geometric measures
  real(real_kind), allocatable :: dvEdge(:)                     ! Distance in meters between the vertices that saddle a given edge
  real(real_kind), allocatable :: dcEdge(:)                     ! Distance in meters between the cells that saddle a given edge
  real(real_kind), allocatable :: areaCell(:)                   ! Area in square meters for a given cell of the primary mesh
  real(real_kind), allocatable :: areaEdge(:)                   ! Area in square meters for a given edge point
  real(real_kind), allocatable :: areaTriangle(:)               ! Area in square meters for a given triangle of the dual mesh
!  real(real_kind), allocatable :: kiteAreasOnVertex(:,:)        ! The intersection area of areaTriangle with each cell that radiates from a given vertex
  real(real_kind), allocatable :: angleEdge(:)                  ! Angle in radians an edge’s normal vector makes with the local eastward direction
  real(real_kind), allocatable :: edgeNormalVectors(:,:  )      ! The unit vector normal to the edge and tangent to the sphere
  real(real_kind), allocatable :: cellTangentPlane(:,:,:)       ! 2 orthogonal unit vectors in the tangent plane of each cell. The first unit vector is chosen to point toward the center of the first edge on the cell.
  real(real_kind), allocatable :: localVerticalUnitVectors(:,:) ! the unit normal vector of the tangent plane at the center of each cell
  
  integer, allocatable :: nSignEdge(:,:)
  real(real_kind) totalArea
  ! Indices
  integer, allocatable :: indexToCellID(:)                ! Global cell ID for all cell centers
  integer, allocatable :: indexToEdgeID(:)                ! Global edge ID for all edge locations
  integer, allocatable :: indexToVertexID(:)              ! Global vertex ID for all cell vertices
  integer, allocatable :: cellsOnCell(:,:)                ! Cell indices that surround a given cell
  integer, allocatable :: cellsOnEdge(:,:)                ! Cell indices that saddle a given edge
!  integer, allocatable :: cellsOnVertex(:,:)              ! Cell indices that radiate from a given vertex
  integer, allocatable :: edgesOnCell(:,:)                ! Edge indices that surround a given cell
!  integer, allocatable :: edgesOnEdge(:,:)                ! Edge indices that are used to reconstruct tangential velocities
!  integer, allocatable :: edgesOnVertex(:,:)              ! Edge indices that radiate from a given vertex
  integer, allocatable :: verticesOnCell(:,:)             ! Vertex indices that surround a given cell
  integer, allocatable :: verticesOnEdge(:,:)             ! Vertex indices that saddle a given edge
  ! Weights
!  real(real_kind), allocatable :: weightsOnEdge(:,:)        ! Weights to reconstruct tangential velocities
  real(real_kind), allocatable :: coeffs_reconstruct(:,:,:) ! Weights to reconstruct vector on cells
  
  real(real_kind), allocatable :: fCell(:)                ! Coriolis coefficients on a given cell
  real(real_kind), allocatable :: fEdge(:)                ! Coriolis coefficients on a given Edge
  real(real_kind), allocatable :: fVertex(:)              ! Coriolis coefficients on a given vertex
  real(real_kind), allocatable :: tanLat(:)               ! tan(latitude)
  
  ! size is expanded to mpi_num_ncells_halo
  real(real_kind), allocatable :: bottomDepth(:)          ! bottom depth, down is positive
  !bottomDepth is reallocated, can not use acc declare create
  
  real(real_kind), allocatable :: angleEdgeCos(:)                  
  real(real_kind), allocatable :: angleEdgeSin(:)                  

  real(real_kind), allocatable :: dHdx_read(:)
  real(real_kind), allocatable :: dHdy_read(:)

contains

  subroutine mesh_init()

    integer iCell, iVertex, i

    if(mpi_rank==0) then
      call fiona_open_dataset('mesh', file_path=mesh_file_path)
      call fiona_get_dim('mesh', 'nCells',       size=mpi_total_nCells)
      call fiona_get_dim('mesh', 'nEdges',       size=mpi_total_nEdges)
      call fiona_get_dim('mesh', 'nVertices',    size=mpi_total_nVertices)
      call fiona_get_dim('mesh', 'vertexDegree', size=vertexDegree)
      call fiona_get_dim('mesh', 'maxEdges',     size=maxEdges)
      call fiona_get_dim('mesh', 'maxEdges2',    size=maxEdges2)
    end if

    if(mpi_procs>1) then
      call mpi_pack_bcast(vertexDegree,maxEdges,maxEdges2)
      call mpi_bcast_process_partition(nCells,maxEdges)
    else
      nCells=mpi_total_nCells
      nEdges=mpi_total_nEdges
      nVertices=mpi_total_nVertices
    end if
    
    allocate(nEdgesOnCell(nCells))
    allocate(latCell(nCells))
    allocate(lonCell(nCells))
    allocate(tanLat(nCells))
    allocate(xCell(nCells))
    allocate(yCell(nCells))
    allocate(zCell(nCells))
    allocate(areaCell(nCells))
    allocate(indexToCellID(nCells))
    allocate(cellsOnCell(maxEdges,nCells))
    allocate(edgesOnCell(maxEdges,nCells))
    allocate(verticesOnCell(maxEdges,nCells))
    allocate(coeffs_reconstruct(3,maxEdges,nCells))
    allocate(cellTangentPlane(3,2,nCells))
    allocate(localVerticalUnitVectors(3,nCells))
    allocate(bottomDepth(nCells))
    allocate(dHdx_read(nCells))
    allocate(dHdy_read(nCells))
    
    if(mpi_procs>1) then
      if(mpi_rank==0) then
        write(*,*) "mpi_total_nCells=",mpi_total_nCells
  
        allocate(mpi_integer_buf_nCells(mpi_total_nCells))
  
        call fiona_input('mesh', 'nEdgesOnCell',      mpi_integer_buf_nCells)
        call mpi_input_exchange_root(mpi_integer_buf_nCells,nEdgesOnCell,nCells,'oncell')
  
        call fiona_input('mesh', 'indexToCellID',     mpi_integer_buf_nCells)
        call mpi_input_exchange_root(mpi_integer_buf_nCells,indexToCellID,nCells,'oncell')

        deallocate(mpi_integer_buf_nCells)

        allocate(mpi_real_buf_nCells(mpi_total_nCells))
  
        call fiona_input('mesh', 'latCell',           mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,latCell,nCells,'oncell')
  
        call fiona_input('mesh', 'lonCell',           mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,lonCell,nCells,'oncell')
  
        call fiona_input('mesh', 'xCell',             mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,xCell,nCells,'oncell')
  
        call fiona_input('mesh', 'yCell',             mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,yCell,nCells,'oncell')
  
        call fiona_input('mesh', 'zCell',             mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,zCell,nCells,'oncell')
  
        call fiona_input('mesh', 'areaCell',          mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,areaCell,nCells,'oncell')
  
        
        call fiona_input('mesh', 'bottomDepth',       mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,bottomDepth,nCells,'oncell')
  
        call fiona_input('mesh', 'dHdx',              mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,dHdx_read,nCells,'oncell')
  
        call fiona_input('mesh', 'dHdy',              mpi_real_buf_nCells)
        call mpi_input_exchange_root(mpi_real_buf_nCells,dHdy_read,nCells,'oncell')

        deallocate(mpi_real_buf_nCells)
  
        allocate(mpi_integer_buf_nCells_2d(maxEdges,mpi_total_nCells))
  
        call fiona_input('mesh', 'cellsOnCell',       mpi_integer_buf_nCells_2d)
        call mpi_input_exchange_root(mpi_integer_buf_nCells_2d,cellsOnCell,maxEdges,nCells,'oncell')
        
        ! for mesh partition, broadcast all cellsOncell
        allocate(mpi_cellsOnCell_all(maxEdges,mpi_total_nCells))
        mpi_cellsOnCell_all(:,:)=reshape(mpi_integer_buf_nCells_2d,(/maxEdges,mpi_total_nCells/))
        call mpi_bcast_2d(mpi_cellsOnCell_all)
  
        call fiona_input('mesh', 'edgesOnCell',       mpi_integer_buf_nCells_2d)
        call mpi_input_exchange_root(mpi_integer_buf_nCells_2d,edgesOnCell,maxEdges,nCells,'oncell')
  
        ! for edge partition, broadcast all edgesOncell
        allocate(mpi_edgesOnCell_all(maxEdges,mpi_total_nCells))
        mpi_edgesOnCell_all(:,:)=reshape(mpi_integer_buf_nCells_2d,(/maxEdges,mpi_total_nCells/))
        call mpi_bcast_2d(mpi_edgesOnCell_all)
  
        call fiona_input('mesh', 'verticesOnCell',    mpi_integer_buf_nCells_2d)
        call mpi_input_exchange_root(mpi_integer_buf_nCells_2d,verticesOnCell,maxEdges,nCells,'oncell')
      
        deallocate(mpi_integer_buf_nCells_2d)
  
      else
        call mpi_input_exchange_nonroot(mpi_integer_buf_nCells,nEdgesOnCell,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_integer_buf_nCells,indexToCellID,nCells,'oncell')

        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,latCell,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,lonCell,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,xCell,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,yCell,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,zCell,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,areaCell,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,bottomDepth,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,dHdx_read,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_real_buf_nCells,dHdy_read,nCells,'oncell')
        call mpi_input_exchange_nonroot(mpi_integer_buf_nCells_2d,cellsOnCell,maxEdges,nCells,'oncell')
  
        allocate(mpi_cellsOnCell_all(maxEdges,mpi_total_nCells))
        call mpi_bcast_2d(mpi_cellsOnCell_all)
  
        call mpi_input_exchange_nonroot(mpi_integer_buf_nCells_2d,edgesOnCell,maxEdges,nCells,'oncell')
        
        ! for edge partition, broadcast all edgesOncell
        allocate(mpi_edgesOnCell_all(maxEdges,mpi_total_nCells))
        call mpi_bcast_2d(mpi_edgesOnCell_all)
  
        call mpi_input_exchange_nonroot(mpi_integer_buf_nCells_2d,verticesOnCell,maxEdges,nCells,'oncell')
      end if

    else
      !execute serial code
  
        call fiona_input('mesh', 'nEdgesOnCell',      nEdgesOnCell)
  
        call fiona_input('mesh', 'latCell',           latCell)
  
        call fiona_input('mesh', 'lonCell',           lonCell)
  
        call fiona_input('mesh', 'xCell',             xCell)
  
        call fiona_input('mesh', 'yCell',             yCell)
  
        call fiona_input('mesh', 'zCell',             zCell)
  
        call fiona_input('mesh', 'areaCell',          areaCell)
  
        call fiona_input('mesh', 'indexToCellID',     indexToCellID)
        
        call fiona_input('mesh', 'bottomDepth',       bottomDepth)
  
        call fiona_input('mesh', 'dHdx',              dHdx_read)
  
        call fiona_input('mesh', 'dHdy',              dHdy_read)
  
        call fiona_input('mesh', 'cellsOnCell',       cellsOnCell)
        
        call fiona_input('mesh', 'edgesOnCell',       edgesOnCell)
  
        call fiona_input('mesh', 'verticesOnCell',    verticesOnCell)
    end if



    if(mpi_procs>1) then
      call mpi_get_index_edges_vertices(edgesOnCell,maxEdges,nCells,nEdges,mpi_total_nEdges,'onedge')
    end if


    allocate(nEdgesOnEdge(nEdges))
    allocate(latEdge(nEdges))
    allocate(lonEdge(nEdges))
    allocate(xEdge(nEdges))
    allocate(yEdge(nEdges))
    allocate(zEdge(nEdges))
    allocate(dvEdge(nEdges))
    allocate(dcEdge(nEdges))
    allocate(angleEdge(nEdges))
    allocate(indexToEdgeID(nEdges))
    allocate(cellsOnEdge(2,nEdges))
!    allocate(edgesOnEdge(maxEdges2,nEdges))
    allocate(verticesOnEdge(2,nEdges))
!    allocate(weightsOnEdge(maxEdges2,nEdges))
    allocate(edgeNormalVectors(3,nEdges))

    if(mpi_procs>1) then
      if(mpi_rank==0) then

        write(*,*) "mpi_total_nEdges=",mpi_total_nEdges
        allocate(mpi_real_buf_nEdges(mpi_total_nEdges))
        call fiona_input('mesh', 'lonEdge',           mpi_real_buf_nEdges)
        call mpi_input_exchange_root(mpi_real_buf_nEdges,lonEdge,nEdges,'onedge')

        call fiona_input('mesh', 'latEdge',           mpi_real_buf_nEdges)
        call mpi_input_exchange_root(mpi_real_buf_nEdges,latEdge,nEdges,'onedge')

        call fiona_input('mesh', 'xEdge',             mpi_real_buf_nEdges)
        call mpi_input_exchange_root(mpi_real_buf_nEdges,xEdge,nEdges,'onedge')

        call fiona_input('mesh', 'yEdge',             mpi_real_buf_nEdges)
        call mpi_input_exchange_root(mpi_real_buf_nEdges,yEdge,nEdges,'onedge')

        call fiona_input('mesh', 'zEdge',             mpi_real_buf_nEdges)
        call mpi_input_exchange_root(mpi_real_buf_nEdges,zEdge,nEdges,'onedge')

        call fiona_input('mesh', 'dvEdge',            mpi_real_buf_nEdges)
        call mpi_input_exchange_root(mpi_real_buf_nEdges,dvEdge,nEdges,'onedge')

        call fiona_input('mesh', 'dcEdge',            mpi_real_buf_nEdges)
        call mpi_input_exchange_root(mpi_real_buf_nEdges,dcEdge,nEdges,'onedge')

        call fiona_input('mesh', 'angleEdge',         mpi_real_buf_nEdges)
        call mpi_input_exchange_root(mpi_real_buf_nEdges,angleEdge,nEdges,'onedge')
        deallocate(mpi_real_buf_nEdges)

        allocate(mpi_integer_buf_nEdges(mpi_total_nEdges))
        call fiona_input('mesh', 'indexToEdgeID',     mpi_integer_buf_nEdges)
        call mpi_input_exchange_root(mpi_integer_buf_nEdges,indexToEdgeID,nEdges,'onedge')
      
        call fiona_input('mesh', 'nEdgesOnEdge',      mpi_integer_buf_nEdges)
        call mpi_input_exchange_root(mpi_integer_buf_nEdges,nEdgesOnEdge,nEdges,'onedge')
        deallocate(mpi_integer_buf_nEdges)

        allocate(mpi_integer_buf_nEdges_2d_2(2,mpi_total_nEdges))
        call fiona_input('mesh', 'cellsOnEdge',       mpi_integer_buf_nEdges_2d_2)
        call mpi_input_exchange_root(mpi_integer_buf_nEdges_2d_2,cellsOnEdge,2,nEdges,'onedge')

        call fiona_input('mesh', 'verticesOnEdge',    mpi_integer_buf_nEdges_2d_2)
        call mpi_input_exchange_root(mpi_integer_buf_nEdges_2d_2,verticesOnEdge,2,nEdges,'onedge')

!        write(*,*) "maxEdges2=",maxEdges2,",nEdges=",nEdges,",mpi_rank=",mpi_rank
!        allocate(mpi_real_buf_nEdges_2d_2(maxEdges2,mpi_total_nEdges))
!        call fiona_input('mesh', 'weightsOnEdge',     mpi_real_buf_nEdges_2d_2)
!        call mpi_input_exchange_root(mpi_real_buf_nEdges_2d_2,weightsOnEdge,maxEdges2,nEdges,'onedge')
       
!        allocate(mpi_integer_buf_nEdges_2d_max(maxEdges2,mpi_total_nEdges))
!        call fiona_input('mesh', 'edgesOnEdge',       mpi_integer_buf_nEdges_2d_max)
!        call mpi_input_exchange_root(mpi_integer_buf_nEdges_2d_max,edgesOnEdge,maxEdges2,nEdges,'onedge')

        deallocate(mpi_integer_buf_nEdges_2d_2)
!        deallocate(mpi_real_buf_nEdges_2d_2)
!        deallocate(mpi_integer_buf_nEdges_2d_max)

      else
        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges,lonEdge,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges,latEdge,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges,xEdge,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges,yEdge,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges,zEdge,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges,dvEdge,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges,dcEdge,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges,angleEdge,nEdges,'onedge')

        call mpi_input_exchange_nonroot(mpi_integer_buf_nEdges,indexToEdgeID,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_integer_buf_nEdges,nEdgesOnEdge,nEdges,'onedge')

        call mpi_input_exchange_nonroot(mpi_integer_buf_nEdges_2d_2,cellsOnEdge,2,nEdges,'onedge')
        call mpi_input_exchange_nonroot(mpi_integer_buf_nEdges_2d_2,verticesOnEdge,2,nEdges,'onedge')

!        write(*,*) "maxEdges2=",maxEdges2,",nEdges=",nEdges,",mpi_rank=",mpi_rank
!        call mpi_input_exchange_nonroot(mpi_real_buf_nEdges_2d_2,weightsOnEdge,maxEdges2,nEdges,'onedge')

!        call mpi_input_exchange_nonroot(mpi_integer_buf_nEdges_2d_max,edgesOnEdge,maxEdges2,nEdges,'onedge')
  
      end if
    else
      !excute serial code
        call fiona_input('mesh', 'lonEdge',           lonEdge)

        call fiona_input('mesh', 'latEdge',           latEdge)

        call fiona_input('mesh', 'xEdge',             xEdge)

        call fiona_input('mesh', 'yEdge',             yEdge)

        call fiona_input('mesh', 'zEdge',             zEdge)

        call fiona_input('mesh', 'dvEdge',            dvEdge)

        call fiona_input('mesh', 'dcEdge',            dcEdge)

        call fiona_input('mesh', 'angleEdge',         angleEdge)

        call fiona_input('mesh', 'indexToEdgeID',     indexToEdgeID)
      
        call fiona_input('mesh', 'nEdgesOnEdge',      nEdgesOnEdge)

        call fiona_input('mesh', 'cellsOnEdge',       cellsOnEdge)

        call fiona_input('mesh', 'verticesOnEdge',    verticesOnEdge)

!        call fiona_input('mesh', 'weightsOnEdge',     weightsOnEdge)
       
!        call fiona_input('mesh', 'edgesOnEdge',       edgesOnEdge)
    end if

    if(mpi_procs>1) then
      call mpi_get_index_edges_vertices(verticesOnCell,maxEdges,nCells,nVertices,mpi_total_nVertices,'onvert')
    end if

    allocate(latVertex(nVertices))
    allocate(lonVertex(nVertices))
    allocate(xVertex(nVertices))
    allocate(yVertex(nVertices))
    allocate(zVertex(nVertices))
    allocate(areaTriangle(nVertices))
!    allocate(kiteAreasOnVertex(vertexDegree,nVertices))
    allocate(indexToVertexID(nVertices))
!    allocate(cellsOnVertex(vertexDegree,nVertices))
!    allocate(edgesOnVertex(vertexDegree,nVertices))
  
    if(mpi_procs>1) then
      if(mpi_rank==0) then

        allocate(mpi_real_buf_vertices(mpi_total_nVertices))
        write(*,*) "mpi_total_nVertices=",mpi_total_nVertices
        call fiona_input('mesh', 'lonVertex',         mpi_real_buf_vertices)
        call mpi_input_exchange_root(mpi_real_buf_vertices,lonVertex,nVertices,'onvert')

        call fiona_input('mesh', 'latVertex',         mpi_real_buf_vertices)
        call mpi_input_exchange_root(mpi_real_buf_vertices,latVertex,nVertices,'onvert')

        call fiona_input('mesh', 'xVertex',           mpi_real_buf_vertices)
        call mpi_input_exchange_root(mpi_real_buf_vertices,xVertex,nVertices,'onvert')

        call fiona_input('mesh', 'yVertex',           mpi_real_buf_vertices)
        call mpi_input_exchange_root(mpi_real_buf_vertices,yVertex,nVertices,'onvert')

        call fiona_input('mesh', 'zVertex',           mpi_real_buf_vertices)
        call mpi_input_exchange_root(mpi_real_buf_vertices,zVertex,nVertices,'onvert')
        
        call fiona_input('mesh', 'areaTriangle',      mpi_real_buf_vertices)
        call mpi_input_exchange_root(mpi_real_buf_vertices,areaTriangle,nVertices,'onvert')
        deallocate(mpi_real_buf_vertices)

        allocate(mpi_integer_buf_vertices(mpi_total_nVertices))
        call fiona_input('mesh', 'indexToVertexID',   mpi_integer_buf_vertices)
        call mpi_input_exchange_root(mpi_integer_buf_vertices,indexToVertexID,nVertices,'onvert')

!        allocate(mpi_integer_buf_vertices_2d(vertexDegree,mpi_total_nVertices))
!        call fiona_input('mesh', 'edgesOnVertex',     mpi_integer_buf_vertices_2d)
!        call mpi_input_exchange_root(mpi_integer_buf_vertices_2d,edgesOnVertex,vertexDegree,nVertices,'onvert')

!        call fiona_input('mesh', 'cellsOnVertex',     mpi_integer_buf_vertices_2d)
!        call mpi_input_exchange_root(mpi_integer_buf_vertices_2d,cellsOnVertex,vertexDegree,nVertices,'onvert')
!
!        allocate(mpi_real_buf_vertices_2d(vertexDegree,mpi_total_nVertices))
!        call fiona_input('mesh', 'kiteAreasOnVertex', mpi_real_buf_vertices_2d)
!        call mpi_input_exchange_root(mpi_real_buf_vertices_2d,kiteAreasOnVertex,vertexDegree,nVertices,'onvert')
        
        deallocate(mpi_integer_buf_vertices)
!        deallocate(mpi_integer_buf_vertices_2d)
!        deallocate(mpi_real_buf_vertices_2d)
      else
        call mpi_input_exchange_nonroot(mpi_real_buf_vertices,lonVertex,nVertices,'onvert')
        call mpi_input_exchange_nonroot(mpi_real_buf_vertices,latVertex,nVertices,'onvert')
        call mpi_input_exchange_nonroot(mpi_real_buf_vertices,xVertex,nVertices,'onvert')
        call mpi_input_exchange_nonroot(mpi_real_buf_vertices,yVertex,nVertices,'onvert')
        call mpi_input_exchange_nonroot(mpi_real_buf_vertices,zVertex,nVertices,'onvert')
        call mpi_input_exchange_nonroot(mpi_real_buf_vertices,areaTriangle,nVertices,'onvert')
        call mpi_input_exchange_nonroot(mpi_integer_buf_vertices,indexToVertexID,nVertices,'onvert')
!        call mpi_input_exchange_nonroot(mpi_integer_buf_vertices_2d,edgesOnVertex,vertexDegree,nVertices,'onvert')
!        call mpi_input_exchange_nonroot(mpi_integer_buf_vertices_2d,cellsOnVertex,vertexDegree,nVertices,'onvert')
!        call mpi_input_exchange_nonroot(mpi_real_buf_vertices_2d,kiteAreasOnVertex,vertexDegree,nVertices,'onvert')
      end if
    else

        call fiona_input('mesh', 'lonVertex',         lonVertex)

        call fiona_input('mesh', 'latVertex',         latVertex)

        call fiona_input('mesh', 'xVertex',           xVertex)

        call fiona_input('mesh', 'yVertex',           yVertex)

        call fiona_input('mesh', 'zVertex',           zVertex)
        
        call fiona_input('mesh', 'areaTriangle',      areaTriangle)

        call fiona_input('mesh', 'indexToVertexID',   indexToVertexID)

!        call fiona_input('mesh', 'edgesOnVertex',     edgesOnVertex)

!        call fiona_input('mesh', 'cellsOnVertex',     cellsOnVertex)

!        call fiona_input('mesh', 'kiteAreasOnVertex', kiteAreasOnVertex)
    end if
    
    !call fiona_end_input('mesh')


    if(mpi_procs>1) then
      call mpi_find_halo_cells(cellsOnCell)
      call mpi_get_cell_destination(nCells)
      call mpi_init_graph_index
  
  
      call mpi_find_halo_edges(edgesOnCell,nEdges)
      call mpi_get_edge_destination(nEdges)

!#ifdef MPI_DEBUG
!    call mpi_test_output(TRIM("_cellsOnCell.nc"),cellsOnCell(1:nCells,1:maxEdges),maxEdges)
!#endif

    ! change global cell id to local cell id
      call mpi_adjust_cell_indexes(cellsOnCell)
      call mpi_adjust_cell_indexes(cellsOnEdge)
      call mpi_adjust_edge_indexes(edgesOnCell)
      
      call reallocate(bottomDepth,mpi_num_ncells_halo)
      call mpi_init_arrays(nCells,nEdges)
      deallocate(mpi_cellsOnCell_all)
      deallocate(mpi_edgesOnCell_all)

    else
      mpi_num_ncells_halo=nCells
    end if
    
    ! Derived quantities
    allocate(nCellsOnVertex(nVertices))
    allocate(areaEdge(nEdges))
    allocate(fCell(nCells))
    allocate(fEdge(nEdges))
    allocate(fVertex(nVertices))
    allocate(nSignEdge(maxEdges,nCells))

    nCellsOnVertex = vertexDegree
    areaEdge       = dvEdge(:) * dcEdge(:)
    fCell          = 2.0_rk * omega * sin(latCell(:))
    fEdge          = 2.0_rk * omega * sin(latEdge(:))
    fVertex        = 2.0_rk * omega * sin(latVertex(:))
    tanLat         = tan(latCell(:))
  
    nSignEdge = 0
    do iCell = 1, nCells
      do i = 1, nEdgesOnCell(iCell)
        if (iCell == cellsOnEdge(1,edgesOnCell(i,iCell))) nSignEdge(i,iCell) =  1
        if (iCell == cellsOnEdge(2,edgesOnCell(i,iCell))) nSignEdge(i,iCell) = -1
      end do
    end do

    ! Scale mesh parameters.
    xCell             = xCell             * radius
    yCell             = yCell             * radius
    zCell             = zCell             * radius
    xEdge             = xEdge             * radius
    yEdge             = yEdge             * radius
    zEdge             = zEdge             * radius
    xVertex           = xVertex           * radius
    yVertex           = yVertex           * radius
    zVertex           = zVertex           * radius
    dvEdge            = dvEdge            * radius
    dcEdge            = dcEdge            * radius
    areaCell          = areaCell          * radius**2
    areaTriangle      = areaTriangle      * radius**2
    areaEdge          = areaEdge          * radius**2
!    kiteAreasOnVertex = kiteAreasOnVertex * radius**2

    angleEdgeCos = cos(angleEdge)
    angleEdgeSin = sin(angleEdge)

  end subroutine mesh_init

  subroutine mesh_final()

    if (allocated(nEdgesOnCell))             deallocate(nEdgesOnCell)
    if (allocated(nEdgesOnEdge))             deallocate(nEdgesOnEdge)
    if (allocated(latCell))                  deallocate(latCell)
    if (allocated(lonCell))                  deallocate(lonCell)
    if (allocated(xCell))                    deallocate(xCell)
    if (allocated(yCell))                    deallocate(yCell)
    if (allocated(zCell))                    deallocate(zCell)
    if (allocated(latEdge))                  deallocate(latEdge)
    if (allocated(lonEdge))                  deallocate(lonEdge)
    if (allocated(xEdge))                    deallocate(xEdge)
    if (allocated(yEdge))                    deallocate(yEdge)
    if (allocated(zEdge))                    deallocate(zEdge)
    if (allocated(latVertex))                deallocate(latVertex)
    if (allocated(lonVertex))                deallocate(lonVertex)
    if (allocated(xVertex))                  deallocate(xVertex)
    if (allocated(yVertex))                  deallocate(yVertex)
    if (allocated(zVertex))                  deallocate(zVertex)
    if (allocated(dvEdge))                   deallocate(dvEdge)
    if (allocated(dcEdge))                   deallocate(dcEdge)
    if (allocated(areaCell))                 deallocate(areaCell)
    if (allocated(areaTriangle))             deallocate(areaTriangle)
!    if (allocated(kiteAreasOnVertex))        deallocate(kiteAreasOnVertex)
    if (allocated(angleEdge))                deallocate(angleEdge)
    if (allocated(indexToCellID))            deallocate(indexToCellID)
    if (allocated(indexToEdgeID))            deallocate(indexToEdgeID)
    if (allocated(indexToVertexID))          deallocate(indexToVertexID)
    if (allocated(cellsOnCell))              deallocate(cellsOnCell)
    if (allocated(cellsOnEdge))              deallocate(cellsOnEdge)
!    if (allocated(cellsOnVertex))            deallocate(cellsOnVertex)
    if (allocated(edgesOnCell))              deallocate(edgesOnCell)
!    if (allocated(edgesOnEdge))              deallocate(edgesOnEdge)
!    if (allocated(edgesOnVertex))            deallocate(edgesOnVertex)
    if (allocated(verticesOnCell))           deallocate(verticesOnCell)
    if (allocated(verticesOnEdge))           deallocate(verticesOnEdge)
!    if (allocated(weightsOnEdge))            deallocate(weightsOnEdge)
    if (allocated(nCellsOnVertex))           deallocate(nCellsOnVertex)
    if (allocated(areaEdge))                 deallocate(areaEdge)
    if (allocated(fCell))                    deallocate(fCell)
    if (allocated(fEdge))                    deallocate(fEdge)
    if (allocated(fVertex))                  deallocate(fVertex)
    if (allocated(nSignEdge))                deallocate(nSignEdge)
    if (allocated(edgeNormalVectors))        deallocate(edgeNormalVectors)
    if (allocated(cellTangentPlane))         deallocate(cellTangentPlane)
    if (allocated(localVerticalUnitVectors)) deallocate(localVerticalUnitVectors)
    
  end subroutine mesh_final

end module mesh_mod
