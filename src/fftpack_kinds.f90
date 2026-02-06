module fftpack_kinds
   implicit none(type, external)
   private
   public :: sp, dp, xdp, qp, lk

   !> Single precision real numbers.
   integer, parameter :: sp = selected_real_kind(6)

   !> Double precision real numbers.
   integer, parameter :: dp = selected_real_kind(8)

   !> Extended double precision real numbers
   integer, parameter :: xdp = -1

   !> Quadruple precision real numbers
   integer, parameter :: qp = -1

   !> Default logical kind parameter
   integer, parameter :: lk = kind(.true.)
end module fftpack_kinds
