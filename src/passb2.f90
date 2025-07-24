subroutine passb2(ido, l1, cc, ch, wa1)
   use fftpack_kind
   implicit none
   real(rk) :: cc, ch, ti2, tr2, wa1
   integer :: i, ido, k, l1
   dimension cc(ido, 2, l1), ch(ido, l1, 2), wa1(*)
   if (ido > 2) then
      do concurrent(i=2:ido:2, k=1:l1)
         ch(i - 1, k, 1) = cc(i - 1, 1, k) + cc(i - 1, 2, k)
         tr2 = cc(i - 1, 1, k) - cc(i - 1, 2, k)
         ch(i, k, 1) = cc(i, 1, k) + cc(i, 2, k)
         ti2 = cc(i, 1, k) - cc(i, 2, k)
         ch(i, k, 2) = wa1(i - 1)*ti2 + wa1(i)*tr2
         ch(i - 1, k, 2) = wa1(i - 1)*tr2 - wa1(i)*ti2
      end do
   else
      do concurrent(k=1:l1)
         ch(1, k, 1) = cc(1, 1, k) + cc(1, 2, k)
         ch(1, k, 2) = cc(1, 1, k) - cc(1, 2, k)
         ch(2, k, 1) = cc(2, 1, k) + cc(2, 2, k)
         ch(2, k, 2) = cc(2, 1, k) - cc(2, 2, k)
      end do
   end if
end subroutine passb2
