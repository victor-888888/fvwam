module restart_mod

  use params_mod
  use mesh_mod
  use time_mod
  use fiona_mod
  use log_mod
  use string
  use state_mod
  use mod_mpi_interfaces

  implicit none

  private

  public restart_output
  public restart_input

contains

  subroutine restart_output(wave)
    implicit none
    type(wave_type),  intent(in) :: wave
    if(mpi_rank==0) then

        call fiona_create_dataset('r0', desc=case_name,file_prefix='restart_wave')

        call fiona_add_att('r0', 'source',         'wave model')
        call fiona_add_att('r0', 'author',         'NMEFC')
        call fiona_add_att('r0', 'on_a_sphere',    'YES')
        call fiona_add_att('r0', 'sphere_radius',  radius)
        ! Dimensions
        !update to the global size, not the local size
        call fiona_add_dim('r0', 'nCells',         size=mpi_total_nCells)
        call fiona_add_dim('r0', 'nEdges',         size=mpi_total_nEdges)
        call fiona_add_dim('r0', 'nFre',           size=nFre)
        call fiona_add_dim('r0', 'nDir',           size=nDir)

        call fiona_add_dim('r0', 'TWO',            size=2)
        call fiona_add_dim('r0', 'maxEdges',       size=maxEdges)
        call fiona_add_dim('r0', 'maxEdges2',      size=maxEdges2)
        call fiona_add_dim('r0', 'Time',           add_var=.true.)
        ! Mesh parameters
        call fiona_add_var('r0', 'lonCell',        long_name='Longitude on the cell',                       units='radian', dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('r0', 'latCell',        long_name='Latitude on the cell',                        units='radian', dim_names=['nCells      '],                 data_type='real(8)')
        call fiona_add_var('r0', 'lonEdge',        long_name='Longitude on the edge',                       units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
        call fiona_add_var('r0', 'latEdge',        long_name='Latitude on the edge',                        units='radian', dim_names=['nEdges      '],                 data_type='real(8)')
        ! Dynamical variables
        call fiona_add_var('r0', 'N',              long_name='Wave action on the cell(actually use Directional Wave Spectrum)',                     units='m^2 s',      dim_names=['nCells','nDir','nFre'],            data_type='real(8)')
    
        call fiona_start_output('r0', time_elapsed_seconds(), new_file=.True.)
        
        call fiona_open_dataset('mesh', file_path=mesh_file_path)
        allocate(mpi_real_buf_nCells(mpi_total_nCells))
        call fiona_input('mesh', 'lonCell',           mpi_real_buf_nCells)
        call fiona_output('r0', 'lonCell',    mpi_real_buf_nCells)

        call fiona_input('mesh', 'latCell',           mpi_real_buf_nCells)
        call fiona_output('r0', 'latCell',         mpi_real_buf_nCells)


        allocate(mpi_real_buf_nEdges(mpi_total_nEdges))
        call fiona_input('mesh', 'lonEdge',           mpi_real_buf_nEdges)
        call fiona_output('r0', 'lonEdge',         mpi_real_buf_nEdges)

        call fiona_input('mesh', 'latEdge',           mpi_real_buf_nEdges)
        call fiona_output('r0', 'latEdge',         mpi_real_buf_nEdges)
        
        deallocate(mpi_real_buf_nCells)
        deallocate(mpi_real_buf_nEdges)
        
    end if
    
    call mpi_output_oncell('r0',  'N', wave%N(1:nCells, :, :))
   
    if (mpi_rank==0) then
      call fiona_end_output('r0')
      call log_notice('Restart output.')
    end if

  end subroutine restart_output
  
  subroutine restart_input(wave)
    implicit none
    type(wave_type),  intent(inout) :: wave
    real(real_kind),allocatable,dimension(:,:,:) :: N_all

    if(mpi_rank==0) then
        allocate(N_all(mpi_total_nCells,nDir,nFre)) 
        call fiona_open_dataset('r1', file_path='restart_wave.r0.nc')
        call fiona_input('r1', 'N',  N_all)
    end if
    call mpi_input_oncell(N_all, wave%N(1:nCells,:,:))
    
    wave%N(0,:,:) = 0.0_rk
    
    call mpi_data_exchange_oncell(wave%N(1:mpi_num_ncells_halo,1:nDir,1:nFre))

    if(mpi_rank==0) then
        deallocate(N_all) 
    end if
  end subroutine restart_input

end module restart_mod
