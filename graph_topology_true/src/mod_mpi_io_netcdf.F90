module mod_mpi_io_netcdf
   use netcdf
   use const_mod
   use mod_mpi_test_variables
   implicit none

   interface netcdf_read
      module procedure netcdf_read_dimension
      module procedure netcdf_read_array_1d_wp
      module procedure netcdf_read_array_2d_wp
      module procedure netcdf_read_array_3d_wp
   end interface


   interface netcdf_write
      module procedure netcdf_write_array_0d_dp
      module procedure netcdf_write_array_1d_wp
      module procedure netcdf_write_array_2d_wp
      module procedure netcdf_write_array_2d_int
      module procedure netcdf_write_array_3d_wp
   end interface

contains

   subroutine netcdf_read_dimension(filename, dimname, dimlen)
!==============================================================================
      implicit none
      integer :: ncid, dimid
      integer, intent(inout) :: dimlen
      character(*), intent(in) :: filename, dimname

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_dimid(ncid, dimname, dimid), trim(filename), varname=trim(dimname))
      call netcdf_check(nf90_inquire_dimension(ncid, dimid, len=dimlen), trim(filename), varname=trim(dimname))
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_dimension

!==============================================================================
   subroutine netcdf_read_array_1d_wp(filename, varname, var, d1, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1
      real(real_kind), intent(inout) :: var(d1)
      integer(i1), allocatable, dimension(:) :: var_i1
      integer(i2), allocatable, dimension(:) :: var_i2
      integer(i4), allocatable, dimension(:) :: var_i4
      real(real_kind_4), allocatable, dimension(:) :: var_sp
      real(real_kind_8), allocatable, dimension(:) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1))
         var = real(var_i1, real_kind)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1))
         var = real(var_i2, real_kind)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1))
         var = real(var_i4,real_kind)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = real(var_sp, real_kind)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/1, ts/), count=(/d1, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = real(var_dp, real_kind)
         deallocate (var_dp)
      case default
         write (*, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (*, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_1d_wp

!==============================================================================
   subroutine netcdf_read_array_2d_wp(filename, varname, var, d1, d2, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1, d2
      real(real_kind), intent(inout) :: var(d1, d2)
      integer(i1), allocatable, dimension(:, :) :: var_i1
      integer(i2), allocatable, dimension(:, :) :: var_i2
      integer(i4), allocatable, dimension(:, :) :: var_i4
      real(real_kind_4), allocatable, dimension(:, :) :: var_sp
      real(real_kind_8), allocatable, dimension(:, :) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1, d2))
         var = real(var_i1, real_kind)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = real(var_i2, real_kind)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1, d2))
         var = real(var_i4, real_kind)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = real(var_sp, real_kind)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1, d2))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = real(var_dp, real_kind)
         deallocate (var_dp)
      case default
         write (*, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (*, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_2d_wp

!==============================================================================
   subroutine netcdf_read_array_3d_wp(filename, varname, var, d1, d2, d3, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, xtype
      integer, intent(in), optional :: ts
      character(*), intent(in) :: filename, varname
      integer, intent(in) :: d1, d2, d3
      real(real_kind), intent(inout) :: var(d1, d2, d3)
      integer(i1), allocatable, dimension(:, :, :) :: var_i1
      integer(i2), allocatable, dimension(:, :, :) :: var_i2
      integer(i4), allocatable, dimension(:, :, :) :: var_i4
      real(real_kind_4), allocatable, dimension(:, :, :) :: var_sp
      real(real_kind_8), allocatable, dimension(:, :, :) :: var_dp

      call netcdf_check(nf90_open(trim(filename), nf90_share, ncid), trim(filename), varname=trim(filename))
      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename), varname=trim(varname))
      call netcdf_check(nf90_inquire_variable(ncid, varid, xtype=xtype), trim(filename), varname=trim(varname))
      select case (xtype)
      case (nf90_byte)
         allocate (var_i1(d1, d2, d3))
         var = real(var_i1,real_kind)
         deallocate (var_i1)
      case (nf90_short)
         allocate (var_i2(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_i2, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_i2), trim(filename), varname=trim(varname))
         end if
         var = real(var_i2, real_kind)
         deallocate (var_i2)
      case (nf90_int)
         allocate (var_i4(d1, d2, d3))
         var = real(var_i4, real_kind)
         deallocate (var_i4)
      case (nf90_float)
         allocate (var_sp(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_sp, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_sp), trim(filename), varname=trim(varname))
         end if
         var = real(var_sp, real_kind)
         deallocate (var_sp)
      case (nf90_double)
         allocate (var_dp(d1, d2, d3))
         if (present(ts)) then
            call netcdf_check(nf90_get_var(ncid, varid, var_dp, &
                                           start=(/1, 1, 1, ts/), count=(/d1, d2, d3, 1/)), trim(filename), varname=trim(varname))
         else
            call netcdf_check(nf90_get_var(ncid, varid, var_dp), trim(filename), varname=trim(varname))
         end if
         var = real(var_dp, real_kind)
         deallocate (var_dp)
      case default
         write (*, *) 'ERROR! : read '//trim(varname)//' from ' &
            //trim(filename)
         write (*, *) 'unknown nf90 type = ', xtype
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename), varname=trim(filename))
   end subroutine netcdf_read_array_3d_wp

!==============================================================================
   subroutine netcdf_check(status, filename, varname, errorFlag)
!==============================================================================
      implicit none
      integer, intent(in) :: status
      integer, intent(inout), optional :: errorFlag
      character(*), intent(in) :: filename
      character(*), intent(in), optional :: varname

      if (status /= nf90_noerr) then
         write (MPI_LOG_UNIT, *) 'ERROR! : '//trim(nf90_strerror(status))
         if (present(errorFlag)) then
            errorFlag = errorFlag + status
         else
            write (MPI_LOG_UNIT, *) filename, " var/dim/file name is ", varname
            write (MPI_LOG_UNIT, *) "*** nmefc mcom io failed ***"
            stop
         end if
      end if

   end subroutine netcdf_check

!==============================================================================
   subroutine netcdf_write_array_0d_dp(filename, varname, var, fp, ts, nc_attr1)
!==============================================================================
      implicit none
      integer :: ncid, varid, dimid_t, status, status_var, status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      real(real_kind_8), intent(in) :: var
      integer(i1) :: var_i1(1)
      integer(i2) :: var_i2(1)
      integer(i4) :: var_i4(1)
      real(real_kind_4) :: var_sp(1)
      real(real_kind_8) :: var_dp(1)
      logical :: file_exist
      type(nc_attr) :: nc_attr1

      inquire (file=trim(filename), exist=file_exist)
      status = nf90_open(trim(filename), nf90_write, ncid)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               write (*, *) trim(varname)//' not exist in '//trim(filename)//', we create it now'
               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  write (*, *) 'time not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if
               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/dimid_t/), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/dimid_t/), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/dimid_t/), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/dimid_t/), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/dimid_t/), varid), trim(filename))
               case default
                  write (*, *) 'ERROR! : define '//trim(varname)//' for ' &
                     //trim(filename)
                  write (*, *) 'unknown nf90 type = '//fp
                  write (*, *) 'support nf90 type is :'
                  write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
                  write (*, *) 'i4 for nf90_int, sp for nf90_float'
                  write (*, *) 'dp for nf90_double'
                  stop
               end select
               call netcdf_check(nf90_put_att(ncid, varid, "units", trim(nc_attr1%var_units)), trim(filename))
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            write (*, *) trim(filename)//' open '//trim(filename)//' failed!'
            stop
         end if
      else
         write (*, *) 'Create output file '//trim(filename)
         call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename), varname=trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                     time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/dimid_t/), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/dimid_t/), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/dimid_t/), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/dimid_t/), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/dimid_t/), varid), trim(filename))
         case default
            write (*, *) 'ERROR! : define '//trim(varname)//' for ' &
               //trim(filename)
            write (*, *) 'unknown nf90 type = '//fp
            write (*, *) 'support nf90 type is :'
            write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (*, *) 'i4 for nf90_int, sp for nf90_float'
            write (*, *) 'dp for nf90_double'
            stop
         end select
         call netcdf_check(nf90_put_att(ncid, varid, "units", trim(nc_attr1%var_units)), trim(filename))
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      ! call netcdf_check(nf90_inquire_variable(ncid, varid), trim(filename))
      select case (trim(fp))
      case ('i1')
         var_i1(1) = int(var, 1)
      case ('i2')
         var_i2(1) = int(var, 2)
      case ('i4')
         var_i4(1) = int(var, 4)
      case ('sp')
         var_sp(1) = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, &
                                        start=(/ts/), count=(/1/)), trim(filename))
      case ('dp')
         var_dp(1) = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, &
                                        start=(/ts/), count=(/1/)), trim(filename))
      case default
         write (MPI_LOG_UNIT, *) 'ERROR! : write '//trim(varname)//' for ' &
            //trim(filename)
         write (*, *) 'unknown nf90 type = '//fp
         write (*, *) 'support nf90 type is :'
         write (*, *) 'i1 for nf90_byte, i2 for nf90_short'
         write (*, *) 'i4 for nf90_int, sp for nf90_float'
         write (*, *) 'dp for nf90_double'
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine netcdf_write_array_0d_dp

!==============================================================================
   subroutine netcdf_write_array_1d_wp(filename, varname, var, fp, d1, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, dimid_d1, dimid_t, status, status_var, status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      character(lc) :: dimname1
      integer, intent(in) :: d1
      real(real_kind), intent(in) :: var(d1)
      integer(i1), allocatable, dimension(:) :: var_i1
      integer(i2), allocatable, dimension(:) :: var_i2
      integer(i4), allocatable, dimension(:) :: var_i4
      real(real_kind_4), allocatable, dimension(:) :: var_sp
      real(real_kind_8), allocatable, dimension(:) :: var_dp
      logical :: file_exist

      dimname1 = 'nCells'

      inquire (file = trim(filename), exist = file_exist)
      status = nf90_open(trim(filename), nf90_write, ncid)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               write (MPI_LOG_UNIT, *) trim(varname)//' not exist in '//trim(filename)//', we create it now'
               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, trim(dimname1), dimid_d1)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) trim(dimname1)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) 'time not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if
               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/dimid_d1, dimid_t/), varid), trim(filename))
               case default
                  write (MPI_LOG_UNIT, *) 'ERROR! : define '//trim(varname)//' for ' &
                     //trim(filename)
                  write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
                  write (MPI_LOG_UNIT, *) 'support nf90 type is :'
                  write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
                  write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
                  write (MPI_LOG_UNIT, *) 'dp for nf90_double'
                  stop
               end select
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            write (MPI_LOG_UNIT, *) trim(filename)//' open '//trim(filename)//' failed!'
            stop
         end if
      else
         write (MPI_LOG_UNIT, *) 'Create output file '//trim(filename)
         call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename), varname=trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                     time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/dimid_d1, dimid_t/), varid), trim(filename))
         case default
            write (MPI_LOG_UNIT, *) 'ERROR! : define '//trim(varname)//' for ' &
               //trim(filename)
            write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
            write (MPI_LOG_UNIT, *) 'support nf90 type is :'
            write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
            write (MPI_LOG_UNIT, *) 'dp for nf90_double'
            stop
         end select
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      ! call netcdf_check(nf90_inquire_variable(ncid, varid), trim(filename))
      select case (trim(fp))
      case ('i1')
         allocate (var_i1(d1))
         var_i1 = int(var, 1)
         deallocate (var_i1)
      case ('i2')
         allocate (var_i2(d1))
         var_i2 = int(var, 2)
         deallocate (var_i2)
      case ('i4')
         allocate (var_i4(d1))
         var_i4 = int(var, 4)
         deallocate (var_i4)
      case ('sp')
         allocate (var_sp(d1))
         var_sp = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, &
                                        start=(/1, ts/), count=(/d1, 1/)), trim(filename))
         deallocate (var_sp)
      case ('dp')
         allocate (var_dp(d1))
         var_dp = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, &
                                        start=(/1, ts/), count=(/d1, 1/)), trim(filename))
         deallocate (var_dp)
      case default
         write (MPI_LOG_UNIT, *) 'ERROR! : write '//trim(varname)//' for ' &
            //trim(filename)
         write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
         write (MPI_LOG_UNIT, *) 'support nf90 type is :'
         write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
         write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
         write (MPI_LOG_UNIT, *) 'dp for nf90_double'
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine netcdf_write_array_1d_wp

!==============================================================================
   subroutine netcdf_write_array_2d_wp(filename, varname, var, fp, d1, d2, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, dimid_d1, dimid_d2, dimid_t, status, status_var, status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      character(lc) :: dimname1, dimname2
      integer, intent(in) :: d1, d2
      real(real_kind), intent(in) :: var(d1, d2)
      integer(i1), allocatable, dimension(:, :) :: var_i1
      integer(i2), allocatable, dimension(:, :) :: var_i2
      integer(i4), allocatable, dimension(:, :) :: var_i4
      real(real_kind_4), allocatable, dimension(:, :) :: var_sp
      real(real_kind_8), allocatable, dimension(:, :) :: var_dp
      logical :: file_exist

      dimname1 = 'nCells'

      dimname2 = 'nDir'

      inquire (file=trim(filename), exist=file_exist)
      status = nf90_open(trim(filename), nf90_write, ncid)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               write (MPI_LOG_UNIT, *) trim(varname)//' not exist in '//trim(filename)//', we create it now'

               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, trim(dimname1), dimid_d1)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) trim(dimname1)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, trim(dimname2), dimid_d2)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) trim(dimname2)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) 'time not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if

               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case default
                  write (MPI_LOG_UNIT, *) 'ERROR! : define '//trim(varname)//' for ' &
                     //trim(filename)
                  write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
                  write (MPI_LOG_UNIT, *) 'support nf90 type is :'
                  write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
                  write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
                  write (MPI_LOG_UNIT, *) 'dp for nf90_double'
                  stop
               end select
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            write (MPI_LOG_UNIT, *) trim(filename)//' open '//trim(filename)//' failed!'
            stop
         end if
      else
         write (MPI_LOG_UNIT, *) 'Create output file '//trim(filename)
         call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename), varname=trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                     time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case default
            write (MPI_LOG_UNIT, *) 'ERROR! : define '//trim(varname)//' for ' &
               //trim(filename)
            write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
            write (MPI_LOG_UNIT, *) 'support nf90 type is :'
            write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
            write (MPI_LOG_UNIT, *) 'dp for nf90_double'
            stop
         end select
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      ! call netcdf_check(nf90_inquire_variable(ncid, varid), trim(filename))
      select case (trim(fp))
      case ('i1')
         allocate (var_i1(d1, d2))
         var_i1 = int(var, 1)
         deallocate (var_i1)
      case ('i2')
         allocate (var_i2(d1, d2))
         var_i2 = int(var, 2)
         deallocate (var_i2)
      case ('i4')
         allocate (var_i4(d1, d2))
         var_i4 = int(var, 4)
         deallocate (var_i4)
      case ('sp')
         allocate (var_sp(d1, d2))
         var_sp = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, &
                                        start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename))
         deallocate (var_sp)
      case ('dp')
         allocate (var_dp(d1, d2))
         var_dp = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, &
                                        start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename))
         deallocate (var_dp)
      case default
         write (MPI_LOG_UNIT, *) 'ERROR! : write '//trim(varname)//' for ' &
            //trim(filename)
         write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
         write (MPI_LOG_UNIT, *) 'support nf90 type is :'
         write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
         write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
         write (MPI_LOG_UNIT, *) 'dp for nf90_double'
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine netcdf_write_array_2d_wp

!==============================================================================
   subroutine netcdf_write_array_2d_int(filename, varname, var, fp, d1, d2, ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, dimid_d1, dimid_d2, dimid_t, status, status_var, status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      character(lc) :: dimname1, dimname2
      integer, intent(in) :: d1, d2
      integer, intent(in) :: var(d1, d2)
      integer(i1), allocatable, dimension(:, :) :: var_i1
      integer(i2), allocatable, dimension(:, :) :: var_i2
      integer(i4), allocatable, dimension(:, :) :: var_i4
      real(real_kind_4), allocatable, dimension(:, :) :: var_sp
      real(real_kind_8), allocatable, dimension(:, :) :: var_dp
      logical :: file_exist

      dimname1 = 'nCells'

      dimname2 = 'nDir'

      inquire (file=trim(filename), exist=file_exist)
      status = nf90_open(trim(filename), nf90_write, ncid)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               write (MPI_LOG_UNIT, *) trim(varname)//' not exist in '//trim(filename)//', we create it now'

               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, trim(dimname1), dimid_d1)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) trim(dimname1)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, trim(dimname2), dimid_d2)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) trim(dimname2)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) 'time not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if

               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
               case default
                  write (MPI_LOG_UNIT, *) 'ERROR! : define '//trim(varname)//' for ' &
                     //trim(filename)
                  write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
                  write (MPI_LOG_UNIT, *) 'support nf90 type is :'
                  write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
                  write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
                  write (MPI_LOG_UNIT, *) 'dp for nf90_double'
                  stop
               end select
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            write (MPI_LOG_UNIT, *) trim(filename)//' open '//trim(filename)//' failed!'
            stop
         end if
      else
         write (MPI_LOG_UNIT, *) 'Create output file '//trim(filename)
         call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename), varname=trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                     time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/dimid_d1, dimid_d2, dimid_t/), varid), trim(filename))
         case default
            write (MPI_LOG_UNIT, *) 'ERROR! : define '//trim(varname)//' for ' &
               //trim(filename)
            write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
            write (MPI_LOG_UNIT, *) 'support nf90 type is :'
            write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
            write (MPI_LOG_UNIT, *) 'dp for nf90_double'
            stop
         end select
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      ! call netcdf_check(nf90_inquire_variable(ncid, varid), trim(filename))
      select case (trim(fp))
      case ('i1')
         allocate (var_i1(d1, d2))
         var_i1 = int(var, 1)
         deallocate (var_i1)
      case ('i2')
         allocate (var_i2(d1, d2))
         var_i2 = int(var, 2)
         deallocate (var_i2)
      case ('i4')
         call netcdf_check(nf90_put_var(ncid, varid, var, &
                                        start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename))
      case ('sp')
         allocate (var_sp(d1, d2))
         var_sp = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, &
                                        start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename))
         deallocate (var_sp)
      case ('dp')
         allocate (var_dp(d1, d2))
         var_dp = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, &
                                        start=(/1, 1, ts/), count=(/d1, d2, 1/)), trim(filename))
         deallocate (var_dp)
      case default
         write (MPI_LOG_UNIT, *) 'ERROR! : write '//trim(varname)//' for ' &
            //trim(filename)
         write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
         write (MPI_LOG_UNIT, *) 'support nf90 type is :'
         write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
         write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
         write (MPI_LOG_UNIT, *) 'dp for nf90_double'
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine netcdf_write_array_2d_int

!==============================================================================
   subroutine netcdf_write_array_3d_wp(filename, varname, var, fp, d1, d2, d3,ts)
!==============================================================================
      implicit none
      integer :: ncid, varid, dimid_d1, dimid_d2, dimid_d3,dimid_t, status, status_var, status_def, status_dim
      integer, intent(in) :: ts
      character(*), intent(in) :: filename, varname, fp
      character(lc) :: time, date, zone, timestamp
      character(lc) :: dimname1, dimname2,dimname3
      integer, intent(in) :: d1, d2,d3
      real(real_kind), intent(in) :: var(d1, d2,d3)
      integer(i1), allocatable, dimension(:, :,:) :: var_i1
      integer(i2), allocatable, dimension(:, :,:) :: var_i2
      integer(i4), allocatable, dimension(:, :,:) :: var_i4
      real(real_kind_4), allocatable, dimension(:, :,:) :: var_sp
      real(real_kind_8), allocatable, dimension(:, :,:) :: var_dp
      logical :: file_exist

      dimname1 = 'nCells'

      dimname2 = 'nDir'

      dimname3 = 'nFre'

      inquire (file=trim(filename), exist=file_exist)
      status = nf90_open(trim(filename), nf90_write, ncid)

      if (file_exist) then
         if (status == nf90_noerr) then
            status_var = nf90_inq_varid(ncid, trim(varname), varid)
            if (status_var == nf90_noerr) then
               continue
            else
               write (MPI_LOG_UNIT, *) trim(varname)//' not exist in '//trim(filename)//', we create it now'

               status_def = nf90_redef(ncid)
               status_dim = nf90_inq_dimid(ncid, trim(dimname1), dimid_d1)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) trim(dimname1)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, trim(dimname2), dimid_d2)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) trim(dimname2)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, trim(dimname3), dimid_d3)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) trim(dimname3)//' not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, trim(dimname3), d3, dimid_d3), trim(filename))
               end if
               status_dim = nf90_inq_dimid(ncid, 'time', dimid_t)
               if (status_dim .ne. nf90_noerr) then
                  write (MPI_LOG_UNIT, *) 'time not exist in '//trim(filename)//', we create it now'
                  call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
               end if

               select case (trim(fp))
               case ('i1')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                                 (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
               case ('i2')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                                 (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
               case ('i4')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                                 (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
               case ('sp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                                 (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
               case ('dp')
                  call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                                 (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
               case default
                  write (MPI_LOG_UNIT, *) 'ERROR! : define '//trim(varname)//' for ' &
                     //trim(filename)
                  write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
                  write (MPI_LOG_UNIT, *) 'support nf90 type is :'
                  write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
                  write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
                  write (MPI_LOG_UNIT, *) 'dp for nf90_double'
                  stop
               end select
               call netcdf_check(nf90_enddef(ncid), trim(filename))
            end if
         else
            write (MPI_LOG_UNIT, *) trim(filename)//' open '//trim(filename)//' failed!'
            stop
         end if
      else
         write (MPI_LOG_UNIT, *) 'Create output file '//trim(filename)
         call netcdf_check(nf90_create(trim(filename), nf90_netcdf4, ncid), trim(filename), varname=trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname1), d1, dimid_d1), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname2), d2, dimid_d2), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, trim(dimname3), d3, dimid_d3), trim(filename))
         call netcdf_check(nf90_def_dim(ncid, 'time', nf90_unlimited, dimid_t), trim(filename))
         call netcdf_check(nf90_put_att(ncid, nf90_global, "Creater", &
                                  "National Marine Envirnomental Forecasting Center Mass Conservation Ocean Model"), trim(filename))
         call date_and_time(date=date, time=time, zone=zone)
         timestamp = date(7:8)//"/"//date(5:6)//"/"//date(1:4)//" "// &
                     time(1:2)//":"//time(3:4)//":"//time(5:6)//" "//zone
         call netcdf_check(nf90_put_att(ncid, nf90_global, "TimeStamp", timestamp), trim(filename))
         select case (trim(fp))
         case ('i1')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_byte, &
                                           (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
         case ('i2')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_short, &
                                           (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
         case ('i4')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_int, &
                                           (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
         case ('sp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_float, &
                                           (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
         case ('dp')
            call netcdf_check(nf90_def_var(ncid, trim(varname), nf90_double, &
                                           (/dimid_d1, dimid_d2, dimid_d3,dimid_t/), varid), trim(filename))
         case default
            write (MPI_LOG_UNIT, *) 'ERROR! : define '//trim(varname)//' for ' &
               //trim(filename)
            write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
            write (MPI_LOG_UNIT, *) 'support nf90 type is :'
            write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
            write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
            write (MPI_LOG_UNIT, *) 'dp for nf90_double'
            stop
         end select
         call netcdf_check(nf90_enddef(ncid), trim(filename))
      end if

      call netcdf_check(nf90_inq_varid(ncid, trim(varname), varid), trim(filename))
      ! call netcdf_check(nf90_inquire_variable(ncid, varid), trim(filename))
      select case (trim(fp))
      case ('i1')
         allocate (var_i1(d1, d2,d3))
         var_i1 = int(var, 1)
         deallocate (var_i1)
      case ('i2')
         allocate (var_i2(d1, d2,d3))
         var_i2 = int(var, 2)
         deallocate (var_i2)
      case ('i4')
         allocate (var_i4(d1, d2,d3))
         var_i4 = int(var, 4)
         deallocate (var_i4)
      case ('sp')
         allocate (var_sp(d1, d2,d3))
         var_sp = real(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_sp, &
                                        start=(/1, 1, 1,ts/), count=(/d1, d2, d3,1/)), trim(filename))
         deallocate (var_sp)
      case ('dp')
         allocate (var_dp(d1, d2,d3))
         var_dp = dble(var)
         call netcdf_check(nf90_put_var(ncid, varid, var_dp, &
                                        start=(/1, 1, 1,ts/), count=(/d1, d2, d3,1/)), trim(filename))
         deallocate (var_dp)
      case default
         write (MPI_LOG_UNIT, *) 'ERROR! : write '//trim(varname)//' for ' &
            //trim(filename)
         write (MPI_LOG_UNIT, *) 'unknown nf90 type = '//fp
         write (MPI_LOG_UNIT, *) 'support nf90 type is :'
         write (MPI_LOG_UNIT, *) 'i1 for nf90_byte, i2 for nf90_short'
         write (MPI_LOG_UNIT, *) 'i4 for nf90_int, sp for nf90_float'
         write (MPI_LOG_UNIT, *) 'dp for nf90_double'
         stop
      end select
      call netcdf_check(nf90_close(ncid), trim(filename))
   end subroutine

end module 
