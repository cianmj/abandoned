 !
 !
 !     [    d^2       2m        ]           2m*E
 !     [ - -----  + ------ V(x) ] phi(x) = ------  phi(x)
 !     [   d x^2    hbar^2      ]          hbar^2
 !
 !     rewrite Schroedinger equation: v(x)=2m*V(x)/hbar^2; ee=2m*E/hbar^2
 !
 subroutine finite_difference(rela,d_pot,eigne,eigns)
   use parameters
   use common
   implicit none
   logical, intent(in)      :: rela
   real(rk),    intent(in)  :: d_pot(0:npts)
   real(rk),    intent(out) :: eigne(nsts), eigns(nsts,0:npts)
   logical     :: decr
   integer(ik) :: i, j, nconv, wrt = 0, k, count
   integer*8   :: nerg
   real(rk)    :: infwn1, energy,de,snrg,snrg2,anrg(nsts),nrgsp,nrgjp,nrgmax
   real(rk)    :: xx, wconv, wconv2, d_wfn(0:npts), asts(nsts,0:npts)
   !
   real(rk), external :: wfn
   !
   eigne(:) = 0._rk; eigns(:,:) = 0._rk; d_wfn(:) = 1.d50
   !
   nconv = npts - 1
   !
   decr = .true.
   nrgsp = 1.d3
   nrgjp = 1.d2
   nrgmax = 1.d14
   !
   nerg = nrgsp
! initial energy for scan
   energy = 0.d0     
   if (rela) energy = mass*c**2.d0
   !
   write(*,*) 'Calculating eigenstates  #',nsts
   loop: do j=1,nsts
      count = 0.
      wconv = 0._rk;  wconv2 = 0._rk! initilization
      energy = (1.q0 + 1.q0/nrgmax) * energy    
      if (j > 1) decr = .false.
      !
      step:do i = 1,int(1.d9)
         !
         de = depth/dble(nerg)
         energy = energy + de
         !
         call shootingmethod(rela,d_pot,energy,d_wfn)
            !
  write(99,*) energy,d_wfn(nconv),nerg,decr
         !
         if (wconv==0.) wconv = d_wfn(nconv)
         if (wconv2==0.) wconv2 = d_wfn(nconv)
         !
         if ( ( (wconv>0.).and.(d_wfn(nconv)<0) ).or. &
              & ( (wconv<0.).and.(d_wfn(nconv)>0) ) ) then
            snrg2 = energy
            wconv2 = d_wfn(nconv)
            bisection: do k=1, int(1.d5)
               if (k==int(9.9999999999d4)) write(*,*) 'Not converging'
               xx = (snrg + snrg2)/2._rk
               call shootingmethod(rela,d_pot,xx,d_wfn)
               if (d_wfn(nconv)*wconv<0.d0) then
                  snrg2 = xx
                  wconv2 = d_wfn(nconv)
               else
                  snrg = xx
                  wconv = d_wfn(nconv)
               end if
               if (abs(snrg-snrg2)<epsnrg) exit step
               energy = xx
               if (abs(d_wfn(nconv))<epswf) exit step
            end do bisection
            exit step
         end if
         !
         ! spacing for energy scan
         if (decr) then
            if (abs(d_wfn(nconv))>abs(wconv)) then
               count = count + 1
               if (count<1.d2) then
                  if (nerg < nrgmax) nerg = nerg * nrgjp
!                  energy = energy - 1.d1*de
                  energy = energy - 2.*de
               else
                  de = nrgsp
                  decr = .false.
               end if
               wconv = 0._rk; wconv2 = 0._rk
               cycle step
             else
                if ((energy>depth).and.(nerg > nrgsp)) nerg = nerg * 1.d0/nrgjp
            end if
         else
            if (modulo(j,nwell)==1) then
               if (abs(d_wfn(nconv))<abs(wconv)) then
                  decr = .true.
               else
                  if (nerg > nrgsp) nerg = nerg * 1.d0/nrgjp
               end if
            end if
         end if
         !
         wconv2 = wconv
         wconv = d_wfn(nconv)
         snrg2 = snrg
         snrg = energy
         !
      end do step
      !
      if (rela) then
         write(*,*) j, energy - mass*c**2.d0
      else
         write(*,*) j, energy
      end if
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
 subroutine shootingmethod(rela,d_pot,energy,d_wfn)
   use parameters
   use common
   implicit none
   real(rk), intent(in)  :: d_pot(0:npts),energy
   logical, intent(in)   :: rela
   real(rk), intent(out) :: d_wfn(0:npts)
   integer(ik) :: i
   real(rk)    :: ff(0:npts)
   !
   ff(:) = 2._rk*mass/hbar**2._rk * dx**2._rk * (d_pot(:)-energy) + 2._rk
   !
   if (rela) ff(:) = 2.d0 - dx**2.d0/c**2.d0 * ( (energy-d_pot(:))**2.d0 - &
        & mass**2.d0 * c**4.d0 )
   !
   d_wfn(0) = 0._rk; d_wfn(1) = 1._rk
   do i=1,npts-1
      d_wfn(i+1) = ff(i)*d_wfn(i) - d_wfn(i-1)
   end do
   !
 end subroutine shootingmethod
 !
 !
