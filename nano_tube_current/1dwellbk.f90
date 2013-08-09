!
!     [    d^2       2m        ]           2m*E
!     [ - -----  + ------ V(x) ] phi(x) = ------  phi(x)
!     [   d x^2    hbar^2      ]          hbar^2
!
!     rewrite Schroedinger equation: v(x)=2m*V(x)/hbar^2; ee=2m*E/hbar^2

!
 program dwell
   use parameters
   use common
   implicit none
   real(rk), external :: potent
   real(rk), external :: wfn
   !
   integer(ik) :: i, nconv
   real(rk)    :: xx,infwn1, wconv, sign, fac, de
   real(rk)    :: nerg,energy !,cond1,cond2,old,next,nrg,u0,t0, 
!   real(rk), dimension(0:npts)    :: uu,tt
   real(rk), dimension(0:npts) :: d_wfn, d_pot
   logical :: incr = .true.
   !
   nwell = nwell
   ! Conversion to a.u.
   depth = depth / aueV
   w_length = w_length / aul
   b_length = b_length / aul
   w_ends   = w_ends / aul
   !
   w_ends = 0.001_rk*b_length
!   w_ends = b_length
   !
   !
   ! Total well length and discretization size calculated
   s_length = nwell*(w_length+b_length)-b_length + 2.*w_ends
   dx = s_length/(npts+1.)
   !
   xx = 0._rk
   do i=0,npts
      xx = xx + dx
      d_pot(i) = potent(xx)
   end do
   !
   infwn1 = pi**2./(2._rk*mass*w_length**2.)
   nconv = int((s_length - 0.25_rk * w_ends)/dx)
   !
   nerg = 1.d4              ! spacing for energy scan
   energy = 0._rk          ! initial energy
   wconv = 0._rk;sign = 1._rk; fac = 1._rk !initilization
   loop:do i = 1,int(1.d9)
      de = fac*dble(infwn1/(nerg))
      if (incr) energy = energy + sign*de
      !
      call shootingmethod(d_pot,energy,d_wfn)
write(*,*) energy, d_wfn(nconv)
read(*,*)
!exit loop
      !
      ! Newton-Raphson method
      if (abs(d_wfn(nconv))<1.d5) then
         energy = energy - fac*dble(d_wfn(nconv)*de/(d_wfn(nconv)-wconv))
         incr = .false.
      end if
!      if (abs(d_wfn(nconv)-wconv)<1.d-12) exit loop   
      wconv = d_wfn(nconv)
      if ((incr==.false.).and.(abs(d_wfn(nconv))>=1.d5)) then
         incr = .true.
         fac = 0.5*fac
      end if
      !
      if (abs(d_wfn(nconv))<1.d10) then
         write(20,fmt='(F18.7,F18.0)') energy,d_wfn(nconv)
      end if
write(10,fmt='(I16,F20.10)') i, energy
      !
      ! Step convergence method (linear)
      !
!      if (abs(d_wfn(nconv)-wconv)<1.d-12) exit loop   
!      if (fac<1.d-6) exit loop   
      if (incr) then
         if (abs(d_wfn(nconv))<1.d4) then
            if (wconv==0.) wconv = d_wfn(nconv)
            if (wconv>0.) then
               if (d_wfn(nconv)<0) then
                  sign = -sign
                  fac = .01_rk*fac
                  wconv = d_wfn(nconv)
               end if
            else
               if (d_wfn(nconv)>0) then
                  sign = -sign
                  fac = .01_rk*fac
                  wconv = d_wfn(nconv)
               end if
            end if
         end if
      end if
      !
   end do loop
   !
   write(*,*) energy,infwn1
   !
   write(*,*) d_wfn(nconv)
   call normalize(npts,d_wfn(:))
   xx = 0.
   do i=0,npts
      xx = xx + dx
      write(30,fmt='(F10.5,F10.5,F10.5,F10.5)') xx,energy,d_pot(i),d_wfn(i)
   end do 
   !
 end program dwell
 !
!      d_wfn(i) = real(wfn(i*dx))
!   call normalize(npts,d_wfn(:))
!   do i=0,npts
!      write(10,*) i*dx,real(d_wfn(i))
!   end do 
!!
!   old = 1.d10; nrg = 1.d10 ! arbitrarly large default values
!      u0 = 1._rk
!      t0 = sqrt(2._rk*mass*(d_pot(0)-energy))
!      call numerical(d_pot,energy,u0,t0,uu,tt)
!   call normalize(npts,uu(:))
!      cond1 = sqrt(2._rk*mass*(d_pot(npts)-energy))
!      cond2 = tt(npts)/uu(npts)
!      next = abs(cond1-cond2)
!      if (next<1.d-1) then
!         if (next<old) then
!            old = next
!            nrg = energy
!         else
!            energy = nrg
!            exit loop
!         end if
!      end if
