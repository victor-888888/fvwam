module mask_mod
  use params_mod
  use mesh_mod
  use fiona_mod
  use log_mod

  implicit none
  
  public mask_edge

  !size expands to mpi_num_ncells_halo
  integer, allocatable :: maskCell(:)         ! mask of cell, water = 1, land = 0
  integer, allocatable :: maskEdge(:)         ! mask of edge，water = 1, land = 0                                                  
  integer, allocatable :: maskBdy(:)          ! mask of edge，water = 1, land = 0                                                  


  contains

  subroutine mask_init()
    implicit none
    integer:: iCell, eoc, j, coc
   
    allocate(maskCell(mpi_num_ncells_halo))
    allocate(maskBdy(nCells))
    allocate(maskEdge(nEdges))

    if(mpi_procs>1) then
      if(mpi_rank==0) then
        allocate(mpi_integer_buf_nCells(mpi_total_nCells))
  
        call fiona_open_dataset('mask', file_path=mesh_file_path)
        call fiona_input('mask', 'maskCell', mpi_integer_buf_nCells)
        call mpi_input_exchange_root(mpi_integer_buf_nCells,maskCell,nCells,'oncell')
        deallocate(mpi_integer_buf_nCells)
      else
        call mpi_input_exchange_nonroot(mpi_integer_buf_nCells,maskCell,nCells,'oncell')
      end if
    else
      call fiona_open_dataset('mask', file_path=mesh_file_path)
      call fiona_input('mask', 'maskCell', maskCell)
    end if

    ! Find a bad cell and edge as land
    
    if(mpi_procs>1) then
      call mpi_data_exchange_oncell(maskCell(1:mpi_num_ncells_halo))
    end if

    call  mask_edge()   

    maskBdy(:) = maskCell(1:nCells)
    do iCell = 1, nCells
       do j = 1, nEdgesOnCell(iCell)
          eoc = edgesOnCell(j, iCell)
          if (eoc .gt. 0) then
             if (maskEdge(eoc) .eq. 0) then
                maskBdy(iCell) = 0
                exit
             end if
          end if
       end do
    end do

  end subroutine mask_init

  subroutine mask_final()

    if (allocated(maskCell))                 deallocate(maskCell)
    if (allocated(maskEdge))                 deallocate(maskEdge)

  end subroutine mask_final

  subroutine mask_edge()
    ! This subroutine must be called after the mask_open_bdy_edge subroutine is called
    ! maskEdge = 1, water
    ! maskEdge = 0, land, open bdy edge
    
    ! cell_1 and cell_2 may use cell in communication shadow, so we exchange makCell before using it
    implicit none
    integer iCell, iEdge
    integer cell_1, cell_2
    do iEdge = 1, nEdges
         cell_1 = cellsOnEdge(1,iEdge)
         cell_2 = cellsOnEdge(2,iEdge)

         if (cell_1 .gt. 0 .and. cell_2 .gt. 0 )then
             if (maskCell(cell_1) .eq. 1 .and. maskCell(cell_2) .eq. 1) then
                 maskEdge(iEdge) = 1
             else
                 maskEdge(iEdge) = 0
             end if
         else
              maskEdge(iEdge) = 0
         end if
     end do

  end subroutine mask_edge

end module mask_mod

