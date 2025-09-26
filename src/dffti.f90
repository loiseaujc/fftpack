      subroutine dffti(n, Wsave)
         use fftpack_kind, only: rk
         implicit none
         integer, intent(in) :: n
         real(rk), intent(out) :: wsave(*)
         if (n == 1) return
         call rffti1(n, wsave(n + 1), wsave(2*n + 1))
      end subroutine dffti
