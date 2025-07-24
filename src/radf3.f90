subroutine radf3(ido, l1, cc, ch, wa1, wa2)
   use fftpack_kind
   implicit none
   real(rk) :: cc, ch, ci2, cr2, di2, di3, dr2, dr3, &
               ti2, ti3, tr2, tr3, wa1, wa2
   integer :: i, ic, ido, idp2, k, l1
   dimension ch(ido, 3, l1), cc(ido, l1, 3), wa1(*), wa2(*)
   real(rk), parameter :: taur = -0.5_rk
   ! note: original comment said this was -sqrt(3)/2 but value was 0.86602540378443864676d0
   real(rk), parameter :: taui = sqrt(3.0_rk)/2.0_rk
   do concurrent(k=1:l1)
      cr2 = cc(1, k, 2) + cc(1, k, 3)
      ch(1, 1, k) = cc(1, k, 1) + cr2
      ch(1, 3, k) = taui*(cc(1, k, 3) - cc(1, k, 2))
      ch(ido, 2, k) = cc(1, k, 1) + taur*cr2
   end do
   if (ido == 1) return
   idp2 = ido + 2
   do concurrent(i=3:ido:2, k=1:l1)
      ic = idp2 - i
      dr2 = wa1(i - 2)*cc(i - 1, k, 2) + wa1(i - 1)*cc(i, k, 2)
      di2 = wa1(i - 2)*cc(i, k, 2) - wa1(i - 1)*cc(i - 1, k, 2)
      dr3 = wa2(i - 2)*cc(i - 1, k, 3) + wa2(i - 1)*cc(i, k, 3)
      di3 = wa2(i - 2)*cc(i, k, 3) - wa2(i - 1)*cc(i - 1, k, 3)
      cr2 = dr2 + dr3
      ci2 = di2 + di3
      ch(i - 1, 1, k) = cc(i - 1, k, 1) + cr2
      ch(i, 1, k) = cc(i, k, 1) + ci2
      tr2 = cc(i - 1, k, 1) + taur*cr2
      ti2 = cc(i, k, 1) + taur*ci2
      tr3 = taui*(di2 - di3)
      ti3 = taui*(dr3 - dr2)
      ch(i - 1, 3, k) = tr2 + tr3
      ch(ic - 1, 2, k) = tr2 - tr3
      ch(i, 3, k) = ti2 + ti3
      ch(ic, 2, k) = ti3 - ti2
   end do
end subroutine radf3
