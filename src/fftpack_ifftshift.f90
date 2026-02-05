submodule(fftpack) fftpack_ifftshift

contains

   !> Shifts zero-frequency component to beginning of spectrum for `complex` type.
   pure module function ifftshift_crk(x) result(result)
      complex(dp), intent(in) :: x(:)
      complex(dp), dimension(size(x)) :: result

      result = cshift(x, shift=-ceiling(0.5_dp*size(x)))

   end function ifftshift_crk

   !> Shifts zero-frequency component to beginning of spectrum for `real` type.
   pure module function ifftshift_rrk(x) result(result)
      real(dp), intent(in) :: x(:)
      real(dp), dimension(size(x)) :: result

      result = cshift(x, shift=-ceiling(0.5_dp*size(x)))

   end function ifftshift_rrk

end submodule fftpack_ifftshift
