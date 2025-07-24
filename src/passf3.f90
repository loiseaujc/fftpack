subroutine passf3(ido, l1, cc, ch, wa1, wa2)
   use fftpack_kind
   implicit none
   real(rk) :: cc, ch, ci2, ci3, cr2, cr3, di2, di3, &
   & dr2, dr3, ti2, tr2, wa1, wa2
   integer :: i, ido, k, l1
   dimension cc(ido, 3, l1), ch(ido, l1, 3), wa1(*), wa2(*)
   real(rk), parameter :: taur = -0.5_rk
   real(rk), parameter :: taui = -sqrt(3.0_rk)/2.0_rk
   if (ido /= 2) then
      do concurrent(i=2:ido:2, k=1:l1)
         tr2 = cc(i - 1, 2, k) + cc(i - 1, 3, k)
         cr2 = cc(i - 1, 1, k) + taur*tr2
         ch(i - 1, k, 1) = cc(i - 1, 1, k) + tr2
         ti2 = cc(i, 2, k) + cc(i, 3, k)
         ci2 = cc(i, 1, k) + taur*ti2
         ch(i, k, 1) = cc(i, 1, k) + ti2
         cr3 = taui*(cc(i - 1, 2, k) - cc(i - 1, 3, k))
         ci3 = taui*(cc(i, 2, k) - cc(i, 3, k))
         dr2 = cr2 - ci3
         dr3 = cr2 + ci3
         di2 = ci2 + cr3
         di3 = ci2 - cr3
         ch(i, k, 2) = wa1(i - 1)*di2 - wa1(i)*dr2
         ch(i - 1, k, 2) = wa1(i - 1)*dr2 + wa1(i)*di2
         ch(i, k, 3) = wa2(i - 1)*di3 - wa2(i)*dr3
         ch(i - 1, k, 3) = wa2(i - 1)*dr3 + wa2(i)*di3
      end do
   else
      do concurrent(k=1:l1)
         tr2 = cc(1, 2, k) + cc(1, 3, k)
         cr2 = cc(1, 1, k) + taur*tr2
         ch(1, k, 1) = cc(1, 1, k) + tr2
         ti2 = cc(2, 2, k) + cc(2, 3, k)
         ci2 = cc(2, 1, k) + taur*ti2
         ch(2, k, 1) = cc(2, 1, k) + ti2
         cr3 = taui*(cc(1, 2, k) - cc(1, 3, k))
         ci3 = taui*(cc(2, 2, k) - cc(2, 3, k))
         ch(1, k, 2) = cr2 - ci3
         ch(1, k, 3) = cr2 + ci3
         ch(2, k, 2) = ci2 + cr3
         ch(2, k, 3) = ci2 - cr3
      end do
   end if
end subroutine passf3
