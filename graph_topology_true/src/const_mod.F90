module const_mod

  implicit none

  integer, parameter :: real_kind_4 = 4
  integer, parameter :: real_kind_8 = 8
  integer, parameter :: real_kind   = real_kind_4
  integer, parameter :: rk          = real_kind
  integer, parameter :: int_double   = 8
  integer, parameter :: max_file_path_len = 256
  integer, parameter :: max_name_len = 30

  real(real_kind), parameter :: radius = 6371220_rk
  real(real_kind), parameter :: g      = 9.80616_rk
  real(real_kind), parameter :: omega  = 0.00007292_rk
  real(real_kind), parameter :: pi     = 4.0_rk * atan(1.0)
  real(real_kind), parameter :: pi2    = pi * 2.0_rk
  real(real_kind), parameter :: pi05   = pi * 0.5_rk
  real(real_kind), parameter :: rad    = pi / 180.0_rk
  real(real_kind), parameter :: deg    = 180.0_rk / pi
  real(real_kind), parameter :: eps    = epsilon(1.0_rk)
  real(real_kind), parameter :: missing_value = -1e20

  real(real_kind), parameter :: rho_sea  = 1025.0_rk   ! sea water density, kg*m-3
  real(real_kind), parameter :: rho_air  = 1.293_rk    ! air density, kg*m-3
  
  real(real_kind), parameter :: Pe       = 101000.0_rk ! Ambient pressure of Typhoon, Pa
  real(real_kind), parameter :: inflowA  = 20.0_rk*rad ! Typhoon inflow angle, radian
  
end module const_mod
