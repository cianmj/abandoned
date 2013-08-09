!
 program dwell
   use parameters
   use common
   implicit none
   real(rk), external :: potent
   !
   logical     :: rela
   integer(ik) :: i,j,k,loop
   real(rk)    :: xx, d_pot(0:npts), vv(nsts), xvec(0:npts), t1, t2
   real(rk)    :: eigne(nsts), eigns(nsts,0:npts)
   real(rk)    :: eigneR(nsts), eignsR(nsts,0:npts)
   real(rk),    allocatable :: eigne2(:),eigns2(:,:)
   !
   call cpu_time(t1)
   ! Total well length and discretization size calculated
   s_length = 60.d0
   dx = s_length/dble(npts)
   write(*,*) s_length, dx
   !
   ! alter this values if problem with finite difference convergence
   epswf  = 1.d-10  !7
   epsnrg = 1.d-20 !17
   !
!   xx = 0._rk
   xx = -s_length/2.d0
   do i=0,npts
      xx = xx + dx
      xvec(i) = xx
      d_pot(i) = potent(s_length,xx)
   end do
   write(*,*) 'Max height of potential:',d_pot(0)
   !
   rela=.false.
   call finite_difference(rela,d_pot,eigne,eigns)
   !
   rela=.true.
   call finite_difference(rela,d_pot,eigneR,eignsR)
   if (rela) eigneR(:) = eigneR(:) - mass*c**2.d0
   !
   do i=0,npts
      write(30,fmt='(F20.10,21(F30.20))') xvec(i),d_pot(i), &
           & (eigne(j) + 100.d0*eigns(j,i),j=1,nsts), &
           & (eigneR(j) + 100.d0*eignsR(j,i),j=1,nsts)
      write(40,fmt='(F20.10,21(F30.20))') xvec(i),d_pot(i), &
           & (eigne(j),j=1,nsts),(eigneR(j),j=1,nsts)
   end do
   !
   write(*,*) 'Expected eigenenegies:'
   do i=1,nsts
      write(*,*) sqrt(kk/mass)*(dble(i-1)+0.5d0)
   end do
   !
   write(*,*) 'Amount relativistic :'
   vv(:) = sqrt(1.d0-(mass*c**2.d0/(eigne(:)+mass*c**2.d0))**2.d0)
   write(*,*) vv(:)
   write(*,*) 'Calculation complete'
   !
   call cpu_time(t2)
   write(*,*) t2-t1
   stop 'stopped'
   !
 end program dwell
 !
