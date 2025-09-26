      subroutine dfftb(n, r, Wsave)
         use fftpack_kind, only: rk
         implicit none
         integer, intent(in) :: n
         real(rk), intent(inout) :: r(*)
         real(rk), intent(in) :: wsave(*)
         if (n == 1) return
         call rfftb1(n, r, wsave, wsave(n + 1), wsave(2*n + 1))
      end subroutine dfftb
