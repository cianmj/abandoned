 !
 function potent(x)
   use parameters
   use common
   implicit none
   real(rk)    :: x,potent
   !
   integer(ik) :: i
   real(rk)    :: l_unit
   !
   potent = depth
   !
   l_unit = w_length + b_length
   do i=0,nwell-1
      if ((x>(i*l_unit+w_ends)).and. &
           & (x<(real(i+1)*l_unit+w_ends-b_length))) potent = 0._rk
   end do
   !
   if ( (x<=w_ends).or.(x>=(s_length-w_ends)) ) potent = depth
   if ((x<=0.).or.(x>=s_length)) potent = 100.*depth
   !
 end function potent
 !
 !
! function spotent(x)
!   use parameters
!   implicit none
!   real(rk) :: x, spotent
!   integer(ik) :: i
!   !
!   spotent = 0._rk
!   do i=0,nwell-1
!      if ((x>(w_length/2. + i*l_unit)).and.(x<(w_length/2.+b_length)))
!   end do 
! end function spotent(x)
 !
