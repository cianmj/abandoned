!
!
!     [    d^2       2m        ]           2m*E
!     [ - -----  + ------ V(x) ] phi(x) = ------  phi(x)
!     [   d x^2    hbar^2      ]          hbar^2
!
!     rewrite Schroedinger equation: v(x)=2m*V(x)/hbar^2; ee=2m*E/hbar^2
!
!
 program dwell
   use parameters
   use common
   implicit none
   real(rk), external :: potent
   !
   integer(ik) :: i
   real(rk)    :: xx, d_pot(0:npts)
   real(rk)    :: eigne(nsts)
   real(rk)    :: eigns(nsts,0:npts), eigns2(nsts,0:npts)
   !
   nwell = nwell
   ! Conversion to a.u.
   depth = depth / aueV
   w_length = w_length / aul
   b_length = b_length / aul
   w_ends   = w_ends / aul
   !
   w_ends = 2._rk*b_length
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
   call finite_difference(d_pot,eigne,eigns)
   eigns2(:,:) = eigns(:,:)
   !
   eigns(:,:) = 0.
   call analytic(d_pot,eigne,eigns)
   !
   write(*,*) 'done'
   !
   do i=0,npts
      write(30,fmt='(F12.7,F20.16,F20.16,F20.16,F20.16,F20.16,F20.16,F20.16,F20.16,F20.16)') i*dx,d_pot(i),eigns(1,i),eigns(2,i),eigns(3,i),eigns2(1,i),eigns2(2,i),eigns2(3,i),eigns(7,i),eigns(8,i), 1.0
      write(40,fmt='(F12.5,F20.16,F20.16,F20.16,F20.16,F20.16,F20.16,F20.16,F20.16)') i*dx,eigne(1),eigne(2),eigne(3),eigne(4),eigne(5),eigne(6),eigne(7),eigne(8)
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
