 !
 subroutine analytic(d_pot,eigne,eigns)
   use parameters
   use common
   real(rk),    intent(in)  :: d_pot(0:npts), eigne(nsts)
   real(rk),    intent(out) :: eigns(nsts,0:npts)
   integer(ik) :: i, j, k, info, nums
   real(rk)    :: xx, kn(nsts), kapn(nsts)
   real(rk)    :: AA(2*nwell+1,nsts), BB(2*nwell+1,nsts)
   real(rk)    :: MM(2,3), ipiv
   !
   AA(1,:) = -1._rk; BB(1,:) = 1._rk
   MM(:,:) = 0.
   !
   main: do i = 1,nsts
      xx = 0._rk
      kn(i)   = sqrt(2.*mass*eigne(i))
      kapn(i) = sqrt(2.*mass*(depth-eigne(i)))
      cnct: do j = 1, 2*nwell
         if (mod(j,2)==1) then
           if ((j==1).or.(j==2*nwell)) then
               xx = xx + w_ends
            else
               xx = xx + b_length
            end if
            MM(1,1) = sin(kn(i)*xx)
            MM(1,2) = cos(kn(i)*xx)
            MM(1,3) = AA(j,i)*exp(-kapn(i)*xx) + BB(j,i)*exp(kapn(i)*xx)
            MM(2,1) = kn(i)*cos(kn(i)*xx)
            MM(2,2) = -kn(i)*sin(kn(i)*xx)
            MM(2,3) = -kapn(i)*AA(j,i)*exp(-kapn(i)*xx) + kapn(i)*BB(j,i)*exp(kapn(i)*xx)
         else
            xx = xx + w_length
            MM(1,1) = exp(-kapn(i)*xx)
            MM(1,2) = exp(kapn(i)*xx)
            MM(1,3) = AA(j,i)*sin(kn(i)*xx) + BB(j,i)*cos(kn(i)*xx)
            MM(2,1) = -kapn(i)*exp(-kapn(i)*xx)
            MM(2,2) = kapn(i)*exp(kapn(i)*xx)
            MM(2,3) = kn(i)*cos(kn(i)*xx)*AA(j,i) - kn(i)*sin(kn(i)*xx)*BB(j,i)
         end if
         !
         !!! LAPACK: DGESV(n,nrhs,A,lda,ipiv,B,ldb,info) !!!
         call DGESV(2,1,MM(1:2,1:2),2,ipiv,MM(1:2,3),2,info)
         if (info/=0) then
            write(*,*) "Error solving for coefficients"
            stop
         end if
         !
         AA(j+1,i) = MM(1,3)
         BB(j+1,i) = MM(2,3)
         !
      end do cnct
      !
      xx = 0.
      nums = 1
      do j = 0, npts-1
         xx = j*dx
         if ((xx>w_ends).and.(nums==1)) nums = nums + 1
         if (mod(num,1)==1) then
            if (xx>w_ends + ((nums-1)/2.)*w_length + ((nums-1.)/2)*b_length) nums = nums + 1
         else
            if (xx>w_ends + (nums/2.)*w_length + (nums/2. - 1)*b_length) nums = nums + 1
         end if
         if (mod(num,1)==1) then
            eigns(i,j) = AA(nums,i)*exp(-kapn(i)*xx) + BB(nums,i)*exp(kapn(i)*xx)
         else
            eigns(i,j) = AA(nums,i)*sin(kn(i)*xx) + BB(nums,i)*cos(kn(i)*xx)
         end if
      end do
      !
   end do main
   !
 end subroutine analytic
 !
!
!     [    d^2       2m        ]           2m*E
!     [ - -----  + ------ V(x) ] phi(x) = ------  phi(x)
!     [   d x^2    hbar^2      ]          hbar^2
!
!     rewrite Schroedinger equation: v(x)=2m*V(x)/hbar^2; ee=2m*E/hbar^2
!
 subroutine finite_difference(d_pot,eigne,eigns)
   use parameters
   use common
!   integer(ik), intent(in)  ::
   real(rk),    intent(in)  :: d_pot(0:npts)
   real(rk),    intent(out) :: eigne(nsts), eigns(nsts,0:npts)
   integer(ik) :: i, j, nconv, nerg
   real(rk)    :: infwn1, energy, de, snrg, snrg2
   real(rk)    :: xx, wconv, wconv2, d_wfn(0:npts)
   !
   real(rk), external :: wfn
   !
   infwn1 = pi**2./(2._rk*mass*w_length**2.)
   nconv = int((s_length - 0.1_rk * w_ends)/dx)
   !
   nerg = 1.d4              ! spacing for energy scan
   energy = infwn1/1.d4     ! initial energy
   write(*,*) 'Calculating eigenstate # ;'
   loop: do j=1,nsts
      write(*,*) j
      wconv = 0._rk;  ! initilization
      energy = 1.05_rk * energy    
      !
      step:do i = 1,int(1.d9)
         !
         de = dble(infwn1/nerg)
         energy = energy + de
         !
         call shootingmethod(d_pot,energy,d_wfn)
         if (wconv==0.) wconv = d_wfn(nconv)
         !
         if ( ( (wconv>0.).and.(d_wfn(nconv)<0) ).or. &
              & ( (wconv<0.).and.(d_wfn(nconv)>0) ) ) then
            snrg2 = energy
            wconv2 = d_wfn(nconv)
            bisection: do while(.true.)
               xx = (snrg + snrg2)/2._rk
               call shootingmethod(d_pot,xx,d_wfn)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !              ! if ((abs(wconv)-abs(d_wfn(nconv)))<1.d10) then
 !              if (wconv==0) then
 !                 secant: do while(.true.)
 !                    xx=energy-(energy-snrg)*d_wfn(nconv)/(d_wfn(nconv)-wconv)
 !                    wconv = d_wfn(nconv)
 !                    call shootingmethod(d_pot,xx,d_wfn)
 !                    if ((abs(wconv)-abs(d_wfn(nconv)))<1.d-6) exit secant
 !                    if ((abs(wconv)-abs(d_wfn(nconv)))<1.d-7) exit step
 !                    if (abs(d_wfn(nconv))<1.d-6) exit step
 !                    snrg = energy
 !                    energy = xx
 !                    if (abs(d_wfn(nconv))<1.d-6) exit step
 !                 end do secant
 !              end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               if (d_wfn(nconv)*wconv<0.) then
                  snrg2 = xx
                  wconv2 = d_wfn(nconv)
               else
                  snrg = xx
                  wconv = d_wfn(nconv)
               end if
               energy = xx
               if (abs(d_wfn(nconv))<1.d-1) exit step
            end do bisection
         end if
         !
         wconv = d_wfn(nconv)
         snrg = energy
         !
      end do step
      !
      eigne(j)   = energy
      call normalize(npts,d_wfn(:))
      eigns(j,:) = d_wfn(:)
      !
   end do loop
   !
 end subroutine finite_difference
 !
 !
 !
 subroutine shootingmethod(d_pot,energy,d_wfn)
   use parameters
   use common
   real(rk), intent(in)  :: d_pot(0:npts),energy
   real(rk), intent(out) :: d_wfn(0:npts)
   integer(ik) :: i
   real(rk)    :: f(0:npts)
   !
   d_wfn(0)=0._rk;d_wfn(1)=1._rk
   f(:) = 2._rk * mass * dx**2 * (d_pot(:)-energy) + 2._rk
   do i=1,npts-1
      d_wfn(i+1) = f(i)*d_wfn(i) - d_wfn(i-1)
   end do
   !
 end subroutine shootingmethod
 !
 !
 !
 subroutine numerical(d_pot,energy,u0,t0,uu,tt)
   use parameters
   use common
   implicit none
   real(rk), intent(in)  :: d_pot(0:npts),energy,u0,t0
   real(rk), intent(out) :: uu(0:npts),tt(0:npts)
   integer(ik) :: i
   real(rk)    :: fx(0:npts)
   !
   uu(0)=u0; tt(0)=t0
   fx(:) = 2._rk*mass*(d_pot(:)-energy)
   do i=0,npts-1
      uu(i+1) = uu(i) + dx/2._rk * (2._rk*tt(i)+2._rk*(dx/2._rk)*fx(i)*uu(i))
      tt(i+1) = tt(i) + dx/2._rk * (fx(i)*uu(i)+fx(i+1)*uu(i+1))
   end do
   !
 end subroutine numerical
 !
 !
 !
 function wfn(x)
   use parameters
   use common
   implicit none
!   integer(ik), intent(in) ::
   real(rk), intent(in)    :: x
   integer(ik) :: i, num
   real(rk)    :: y, width, l_unit,yw(npts)
   real(rk) :: wfn
   !
   width = w_length/3._rk      ! gaussian function width in a.u
   l_unit = w_length+b_length
   !
!   do i=0,nwell-1
!      if ( (x>i*l_unit).and.(x<(i+1)*l_unit) ) num = i
!   end do
   !
   do i=0,nwell-1
      yw(i) = w_length/2._rk + i*l_unit
   end do
   !
   wfn = (0._rk,0._rk)
 !  do i=-num,nwell-1-num
!     y = w_length/2._rk + i*l_unit
!      wfn = wfn + exp(-(y-x)**2._rk/width**2)
!   end do
   !
   do i=0,nwell-1
      wfn = wfn + exp(-(yw(i)-x)**2._rk/width**2)
   end do
   !
 end function wfn
