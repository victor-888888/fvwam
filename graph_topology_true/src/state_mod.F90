module state_mod

  use const_mod
  use mesh_mod
  use fre_dir_mod

  implicit none

  private

  public state_init
  public wave_type
  public wave
  public tend_geog
  public tend_theta

  type wave_type
     real(real_kind), allocatable :: N(:,:,:) ! Wave action, N(nCells, nDir, nFre)
  end type wave_type
     
  type(wave_type), allocatable :: wave(:)

  real(real_kind), allocatable :: tend_geog(:,:,:)   ! tendency of geog(nCells, nDir, nFre)
  real(real_kind), allocatable :: tend_theta(:,:,:)  ! tendency of geog(nCells, nDir, nFre)
 

  contains

  subroutine state_init()

    integer i

    if (.not. allocated(wave)) then
      allocate(wave(-1:2))
      do i = lbound(wave, 1), ubound(wave, 1)
        call allocate_wave(wave(i))
      end do
    end if

    if (.not. allocated(tend_geog   )) allocate(tend_geog(nCells, nDir, nFre ))
    if (.not. allocated(tend_theta  )) allocate(tend_theta(nCells, nDir, nFre))


  end subroutine state_init

  subroutine allocate_wave(wave)

    type(wave_type), intent(inout) :: wave

    if (.not. allocated(wave%N  )) allocate(wave%N(0:mpi_num_ncells_halo, nDir, nFre))

  end subroutine allocate_wave

end module state_mod
