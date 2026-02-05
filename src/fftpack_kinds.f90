module fftpack_kinds
   use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64
   implicit none(type, external)
   private

   public :: sp, dp, xdp, qp, lk

   !> Single precision real numbers.
   integer, parameter :: sp = selected_real_kind(6)

   !> Double precision real numbers.
   integer, parameter :: dp = selected_real_kind(15)

   !> Extended double precision real numbers.
   integer, parameter :: xdp = -1

   !> Quadruple precision real numbers.
   integer, parameter :: qp = -1

   !> Default logical kind parameter.
   integer, parameter :: lk = kind(.true.)
end module fftpack_kinds
