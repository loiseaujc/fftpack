module fftpack_kind
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none(type, external)
   public
   integer, parameter :: rk = real64
end module fftpack_kind
