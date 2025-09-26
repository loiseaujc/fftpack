subroutine radb2(ido, l1, cc, ch, wa1)
   use fftpack_kind, only: rk
   implicit none(type, external)
   integer, intent(in) :: ido, l1
   real(rk), intent(in) :: wa1(*), cc(ido, 2, l1)
   real(rk), intent(inout) :: ch(ido, l1, 2)

   real(rk) :: ti2, tr2
   integer :: i, ic, idp2, k
   do concurrent(k=1:l1)
      ch(1, k, 1) = cc(1, 1, k) + cc(ido, 2, k)
      ch(1, k, 2) = cc(1, 1, k) - cc(ido, 2, k)
   end do
   if (ido < 2) return
   if (ido /= 2) then
      idp2 = ido + 2
      do concurrent(i=3:ido:2, k=1:l1)
         ic = idp2 - i
         ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(ic - 1, 2, k)
         tr2 = cc(i - 1, 1, k) - cc(ic - 1, 2, k)
         ch(i, k, 1) = cc(i, 1, k) - cc(ic, 2, k)
         ti2 = cc(i, 1, k) + cc(ic, 2, k)
         ch(i - 1, k, 2) = wa1(i - 2)*tr2 - wa1(i - 1)*ti2
         ch(i, k, 2) = wa1(i - 2)*ti2 + wa1(i - 1)*tr2
      end do
      if (mod(ido, 2) == 1) return
   end if
   do concurrent(k=1:l1)
      ch(ido, k, 1) = cc(ido, 1, k) + cc(ido, 1, k)
      ch(ido, k, 2) = -(cc(1, 2, k) + cc(1, 2, k))
   end do
end subroutine radb2
