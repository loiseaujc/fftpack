      subroutine dsint(n, x, Wsave)
         use fftpack_kind, only: rk
         implicit none
         integer, intent(in) :: n
         real(rk), intent(inout) :: x(*)
         real(rk), intent(in) :: wsave(*)
         integer :: iw1, iw2, iw3, np1
         np1 = n + 1
         iw1 = n/2 + 1
         iw2 = iw1 + np1
         iw3 = iw2 + np1
         call sint1(n, x, wsave, wsave(iw1), wsave(iw2), wsave(iw3))
      end subroutine dsint
