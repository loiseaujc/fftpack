submodule(fftpack) fftpack_fftshift

contains

   !> Shifts zero-frequency component to center of spectrum for `complex` type.
   pure module function fftshift_crk(x) result(result)
      complex(dp), intent(in) :: x(:)
      complex(dp), dimension(size(x)) :: result

      result = cshift(x, shift=-floor(0.5_dp*size(x)))

   end function fftshift_crk

   !> Shifts zero-frequency component to center of spectrum for `real` type.
   pure module function fftshift_rrk(x) result(result)
      real(dp), intent(in) :: x(:)
      real(dp), dimension(size(x)) :: result

      result = cshift(x, shift=-floor(0.5_dp*size(x)))

   end function fftshift_rrk

end submodule fftpack_fftshift
