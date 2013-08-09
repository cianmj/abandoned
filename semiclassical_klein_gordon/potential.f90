 !
 function potent(size, x)
   use parameters
   use common
   implicit none
   real(rk)    :: size, x, potent
   !
   potent = 0.5d0 * kk *(x)**2.d0
   !
 end function potent
 !
 !
 !
 subroutine spotent(v,x)
   use parameters
   use common
   implicit none
   real(rk)    :: x,v
   !
   integer(ik) :: i
   real(rk)    :: l_unit
   !
   v = depth
   !
   l_unit = w_length + b_length
   do i=0,nwell-1
      if ((x>(i*l_unit+w_ends)).and. &
           & (x<(real(i+1)*l_unit+w_ends-b_length))) v = 0._rk
   end do
   !
   if ( (x<=w_ends).or.(x>=(s_length-w_ends)) ) v = depth
   if ((x<=0.).or.(x>=s_length)) v = 10.d0*depth
   !
 end subroutine spotent

