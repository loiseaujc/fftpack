      subroutine dfftf(n, r, Wsave)
         use fftpack_kind, only: rk
         implicit none
         integer, intent(in) :: n
         real(rk), intent(inout) :: r(*)
         real(rk), intent(in) :: wsave(*)
         if (n == 1) return
         call rfftf1(n, r, Wsave, Wsave(n + 1), Wsave(2*n + 1))
      end subroutine dfftf
