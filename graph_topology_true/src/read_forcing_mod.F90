module  read_forcing_mod
    use const_mod
    use params_mod
    use netcdf

    implicit none
    include 'netcdf.inc' 

    contains

    subroutine read_met_forcing_var(nCells, step, u10, v10)
        implicit none
        integer, intent(in) :: nCells
        integer, intent(in) :: step ! time_step
        real(real_kind), intent(out) :: u10(:), v10(:)

        integer ncid, status
        integer u10_id    ,&
                v10_id 

        integer time_dimid, i, time_num, time_num_sum, step2
        character(max_file_path_len) :: wind_file_path
        character(30) :: dimname
        character(4)  :: filenum
        
        if (.not. multiple_nwp_forcing_file) then
            call check_nc(NF90_OPEN(trim(adjustl(nwp_forcing_file_path)),NF90_NoWRITE,ncid)  )
            ! get var ID
            call check_nc(NF90_INQ_VARID(ncid, 'u10', u10_id))
            call check_nc(NF90_INQ_VARID(ncid, 'v10', v10_id))

            call check_nc(NF90_get_var(ncid, u10_id , u10 ,start=(/1,step/) ,count=(/nCells,1/)))
            call check_nc(NF90_get_var(ncid, v10_id , v10 ,start=(/1,step/) ,count=(/nCells,1/)))
            call check_nc(NF90_CLOSE(ncid))
        else
            if (nwp_forcing_file_num .gt. 9999) then
                print*, 'Error: nwp_forcing_file_num should less than  9999!'
                stop
            end if

            time_num_sum = 0
            do i = 1, nwp_forcing_file_num
                if (i .lt. 10) then
                    write(filenum(1:1), "(I1)") 0
                    write(filenum(2:2), "(I1)") 0
                    write(filenum(3:3), "(I1)") 0
                    write(filenum(4:4), "(I1)") i
                else if (i .lt. 100) then
                    write(filenum(1:1),   "(I1)") 0
                    write(filenum(2:2),   "(I1)") 0
                    write(filenum(3:4), "(I2)") i
                else if (i .lt. 1000) then
                    write(filenum(1:1),   "(I1)") 0
                    write(filenum(2:4), "(I3)") i
                else 
                    write(filenum(1:4), "(I4)") i
                end if

                wind_file_path = trim(adjustl(nwp_forcing_file_prefix)) // filenum // '.nc'
                call check_nc(NF90_OPEN(trim(adjustl(wind_file_path)),NF90_NoWRITE,ncid))
                call check_nc(NF90_INQ_DIMID(ncid, 'Time'      ,time_dimid ))
                ! get Time dimension num
                call check_nc(NF90_Inquire_Dimension(ncid, time_dimid,  dimname, time_num)) 
                time_num_sum = time_num_sum + time_num

                

                if (time_num_sum .ge. step) then
                    !print*, 'Reading wind file from: '//wind_file_path
                    step2 = step - (time_num_sum - time_num)
                    !print*, i, step, step2, time_num_sum, time_num, wind_file_path
                    ! get var ID
                    call check_nc(NF90_INQ_VARID(ncid, 'u10', u10_id))
                    call check_nc(NF90_INQ_VARID(ncid, 'v10', v10_id))
                    call check_nc(NF90_get_var(ncid, u10_id , u10 ,start=(/1,step2/) ,count=(/nCells,1/)))
                    call check_nc(NF90_get_var(ncid, v10_id , v10 ,start=(/1,step2/) ,count=(/nCells,1/)))
                    exit
                else
                    if (i .eq. nwp_forcing_file_num)then
                        print*, 'Error: The sum of the times for all files is less than the required, please check the wind files!'
                        stop
                    end if
                end if
            end do
            call check_nc(NF90_CLOSE(ncid))
        end if

    end subroutine read_met_forcing_var
    
    subroutine read_tc_key_parameters(time_size, lonTC, latTC, rmwTC, presTC, utTC, vtTC)
        integer,         intent(in)  :: time_size
        real(real_kind), intent(out) :: lonTC(:), latTC(:)  ! radian
        real(real_kind), intent(out) :: rmwTC(:)            ! meter
        real(real_kind), intent(out) :: presTC(:)           ! Pa
        real(real_kind), intent(out) :: utTC(:), vtTC(:)    ! m/s

        integer :: time_size_data
        integer ncid, status
        character(50) :: dim_name
        integer time_dimid  ,&
                lonTC_id    ,&
                latTC_id    ,&
                rmwTC_id    ,&
                presTC_id   ,&
                utTC_id     ,&
                vtTC_id     
 
        call check_nc( NF90_OPEN(trim(adjustl(tc_forcing_file_path)),NF90_NoWRITE,ncid)  )
        ! get time size
        status=NF90_INQ_DIMID(ncid, 'Time', time_dimid)
        status=NF90_Inquire_Dimension(ncid, time_dimid, dim_name, time_size_data)
        !if (time_size_data /= time_size) then
        if (time_size_data .lt. time_size) then
            print*,"Error: Time size of TC forcing field does not match the input file ", trim(adjustl(tc_forcing_file_path))
            print*,"TC_time_size=",time_size, " Input_file_size=",time_size_data
            stop
        else if(time_size_data .gt. time_size) then
            print*,"Warning: Time size of TC forcing field less than the input file ", trim(adjustl(tc_forcing_file_path))
            print*,"TC_time_size=",time_size, " Input_file_size=",time_size_data
        end if
        ! get var ID
        status=NF90_INQ_VARID(ncid, 'lonTC'  ,lonTC_id )
        status=NF90_INQ_VARID(ncid, 'latTC'  ,latTC_id )
        status=NF90_INQ_VARID(ncid, 'rmwTC'  ,rmwTC_id )
        status=NF90_INQ_VARID(ncid, 'presTC' ,presTC_id)
        status=NF90_INQ_VARID(ncid, 'utTC'   ,utTC_id  )
        status=NF90_INQ_VARID(ncid, 'vtTC'   ,vtTC_id  )
        ! get var
        status=NF90_get_var(ncid, lonTC_id   ,lonTC  ,start=(/1/), count=(/time_size/))
        status=NF90_get_var(ncid, latTC_id   ,latTC  ,start=(/1/), count=(/time_size/))
        status=NF90_get_var(ncid, rmwTC_id   ,rmwTC  ,start=(/1/), count=(/time_size/))
        status=NF90_get_var(ncid, presTC_id  ,presTC ,start=(/1/), count=(/time_size/))
        status=NF90_get_var(ncid, utTC_id    ,utTC   ,start=(/1/), count=(/time_size/))
        status=NF90_get_var(ncid, vtTC_id    ,vtTC   ,start=(/1/), count=(/time_size/))
        call check_nc(NF90_CLOSE(ncid))
    end subroutine read_tc_key_parameters

    ! Check the status after operating the necdf file
    subroutine check_nc(status)
        implicit none
        integer,intent(in) :: status
        if (status /= nf90_noerr) then
            print*, trim( nf90_strerror(status) )
            stop
        end if
    end subroutine check_nc

end module read_forcing_mod
