
 subroutine numerical(d_pot,energy,u0,t0,uu,tt)
   use parameters
   implicit none
   real(rk), intent(in)  :: d_pot(0:npts),energy,u0,t0
   real(rk), intent(out) :: uu(0:npts),tt(0:npts)
   !
   integer(ik) :: i
   real(rk)    :: fx(npts)
   !
   uu(0)=u0; tt(0)=t0
   fx(:) = 2._rk*mass*(d_pot(:)-energy)
   do i=0,npts-1
      u(i+1) = u(i) + dx/2._rk * (2._rk*tt(i)+2._rk*(dx/2._rk)*fx(i)*u(i))
      t(i+1) = t(i) + dx/2._rk * (fx(i)*uu(i)+fx(i+1)*uu(i+1))
   end do
   !
 end subroutine numerical
