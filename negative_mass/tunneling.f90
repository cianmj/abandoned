module param
  implicit none
  ! Constants
  double precision, parameter :: c=137.d0, G=2.41d-43, m0=9.10938d-31
  double precision, parameter :: Msun = 1.9891d30/m0
  double precision, parameter :: pi = 3.14159265358979325d0
  double complex, parameter   :: Im=(0.d0,1.d0)
  !
  ! Paramters !! 2000->1000;
!  integer :: lang
  integer, parameter :: NN = 2000, Nll = 2000
  integer, parameter :: Ndisc = 3*NN
  double precision, parameter :: rdef = 100.d0
!  integer, parameter :: lend = 1000, lmax = 100
  integer, parameter :: lend = 40, lmax = 20
  integer, parameter :: ks1 = 0, ke1 = lend
  double precision, parameter :: mass = 2.d-6
  double precision, parameter :: MM   = 2.d-10*Msun
  double precision, parameter :: Es   = mass*c**2
  double precision, parameter :: En   = 1.15d0*Es
  logical, parameter :: grav=.true.
  integer, parameter :: nfunc=2  !!0-bessel; 1-sin. approx; 2-exp. wave !!
  !
  double precision, parameter :: eps = 1.d-20, rfrac=1.d99 !2
  !
!!!!!!!!!! rfrac is sensitive to discritization !!!!!!
  !
  ! Definitions
  double precision, parameter :: GM   = G*MM/c**2
  double precision, parameter :: rs   = 2.d0*GM
  double precision, parameter :: kr0  = sqrt(En**2-Es**2)/c
  double precision, parameter :: kr00  = sqrt(En**2-Es**2)/c
  !
  contains
    !
    function fr_2(lang,r)
      implicit none
      integer          :: lang
      double precision :: r      ! in rs units
      double precision :: fr_2
      !
      if (r==0.d0) then
         fr_2 = 1.d99
         return
      end if
      !
      if (grav) then
         fr_2 = En**2.d0/(mass*Es) - (1.d0-1.d0/r)* &
              & (c**2.d0+(dble(lang+0.5d0)/(mass*r*rs))**2.d0)
      else
         fr_2=En**2.d0/(mass*Es)-(c**2.d0+(dble(lang+0.5d0)/(mass*r*rs))**2.d0)
      end if
      !
     end function fr_2
     !
     function fser(aa,bb,cc,ff,r)
       implicit none
       double precision :: aa,bb,cc,ff,r
       double complex   :: fser
       !
       if (grav) then
          if ( aa > bb ) then
             fser = 2.d0*aa**(1.5d0)*ff*sqrt(dcmplx(r-aa)) / &
                  & sqrt(dcmplx((aa-bb)*(aa-cc)))
          else
             fser = -2.d0*aa**(1.5d0)*ff*sqrt(dcmplx(aa-r)) / &
                  & sqrt(dcmplx((aa-bb)*(aa+cc)))
          end if
       else
          fser = sqrt(dcmplx(2.d0*aa*(r-aa)))*ff
       end if
       !
     end function fser
     !
end module param
!
!
!
program semi_class
  use param
  implicit none
  !
  integer :: i, j, k
  double precision :: theta, r, dr, prob, dth
  double precision, allocatable :: roots(:,:), rstart(:), rend(:)
  double complex, allocatable :: wav_sum(:,:),wav_sum2(:,:)
  double precision :: besselj(0:lend,1:NN),sinfn(0:lend,1:NN)  
  !
  write(*,*) 'nfunc:',nfunc
  if ( (grav).and.(nfunc.ne.2) ) then
     write(*,*) 'mismatch: grav, nfunc',grav,nfunc
     stop
  end if
  write(*,*) 'Wave vectors :',kr0,kr00
  !
  !
  ! Finding all turning points for angular mom. l= 0 to lend
  !
  write(*,*) 'Finding roots #',lend
  allocate( roots(0:lend,3) )
  call finding_roots(grav,0,lend,roots)
!  call finding_roots(grav,ks1,lend,roots)
  write(*,*) 'Roots of ', ks1
  write(*,*) roots(ks1,:)
  !
  !
  ! Setting the starting point of integration
  !
  allocate( rstart(0:lend), rend(0:lend) )
  call points(rstart,rend)
  write(*,*) 'Start int. range:', rend(0),rend(lend)
  !
  !
  ! Integration of function (v) 
  !  & creation of Legendre polynomials (v)
  !
  allocate( wav_sum(1:NN,1:Nll), wav_sum2(1:NN,1:Nll) )
  call plane_wave(rstart, rend, roots, wav_sum, wav_sum2)
  deallocate( roots )
  !
  r=0.d0
  do i=1,NN
     r = r + rend(lend)/NN
     call sph_bessel(ks1,0,kr00*rs*r,besselj(ks1,i))
     sinfn(ks1,i)= sin(kr00*rs*r-pi*(ks1)/2.d0)/(kr00*rs*r)
  end do
  !
  ! Print data to file(s)....
  !
  write(*,*) 'Printing...'
  open(unit=20,file='wav_in.txt',status='replace',action='write')
  open(unit=21,file='wav_out.txt',status='replace',action='write')
  open(unit=22,file='wav_sum.txt',status='replace',action='write')
  open(unit=23,file='wav_sum2.txt',status='replace',action='write')
  r=0.d0;
  do i=1,NN-1
     r = r + rend(lend)/NN
     theta=0.d0;
     dth = 2.d0*pi/Nll
     do j=1,Nll
        !
        theta =  theta + dth
!if (abs(theta-pi).gt.1.d-4) cycle
        !
        prob = abs(wav_sum(i,j))
        write(22,fmt='(7F30.15)') r, theta,dble(wav_sum(i,j)),&
             & dimag(wav_sum(i,j)),dble(prob),besselj(ks1,i), besselj(ks1,i)*r
        prob = abs(wav_sum2(i,j))
        write(23,fmt='(5F30.15)') r, theta, &
             & dble(wav_sum2(i,j)),dimag(wav_sum2(i,j)),dble(prob)
        !
     end do
  end do
  close( 20 )
  close( 21 )
  close( 22 )
  close( 23 )
  !
  deallocate( rstart, rend, wav_sum, wav_sum2 )
  write(*,*) 'Program complete.'
  stop
  !
end program semi_class
!
!
!
subroutine plane_wave(rstart, rend, roots, wav_sum, wav_sum2)
  use param
  implicit none
  double precision, intent(in) :: rend(0:lend), roots(0:lend,3), rstart(0:lend)
  double complex, intent(out) :: wav_sum(1:NN,1:Nll), wav_sum2(1:NN,1:Nll)
  !
  integer :: i, j, k, i_in, i_out
  double precision :: r, dr, theta, dth
  double precision :: besselj(0:lend,1:NN),sinfn(0:lend,1:NN)  
  double complex :: exp_in, exp_out, fbar
  double complex, allocatable :: fint_in(:,:), fint_out(:,:)
  double precision, allocatable :: lpoly(:,:), bslj(:,:)
  !
  !
  wav_sum(:,:) = 0.d0; wav_sum2(:,:) = 0.d0;
  allocate( lpoly(0:lend,1:Nll), bslj(0:lend, 1:NN) )
  allocate( fint_in(0:lend,0:NN), fint_out(0:lend,0:NN) )
  ! (*,'(I8)',advance="no") i
  write(*,*) 'Integrating in/out going waves... '
  call integrate( .true., rend, roots, fint_in, rstart )
  call integrate( .false., rend, roots, fint_out )
  !
  !
  ! Generating Legendre polynomials
  write(*,*) 'Generating Legendre polynomials... '
  theta = 0.d0
  do j=1,Nll
     dth = 2.d0*pi/Nll
     theta = theta + dth
     call legendre2(lend,theta,lpoly(0:lend,j))
  enddo
  lpoly(0:lend,Nll)=0.d0;   ! Since the 'NN' component is not defined !
  !
  ! Creating Bessel function
  write(*,*) 'Creating Bessel function... '
  do k=ks1,ke1
     r=0.d0
     dr = rend(lend)/NN
     asdf:do i=1,NN
        r= r + dr
        if (r<roots(k,2)) then
           bslj(k,i)=0.d0
           cycle asdf
        endif
!        call sph_bessel(k,nfunc,kr0*rs*r,bslj(k,i))
     end do asdf
  end do
  !
  !
  ! Writing out any of the data
  call write_info(NN,dr,lpoly,bslj,fint_in)
  call write_info(NN,dr,lpoly,bslj,fint_out)
  !
  !
  ! Summing the partial waves for incoming/outgoing parts
  write(*,*) 'Summing partial waves...'
  r = 0.d0
  fsum:do i=1,NN
     r = r + dr
!     if (modulo(i,1000)==1) write(*,'(I8)',advance="no") i
     if (modulo(i,int(NN/10))==1) write(*,'(I8)') i
     plsum:do j=1,Nll
        wav_sum(i,j) = 0.d0
        wav_sum2(i,j) = 0.d0
        lsum:do k=ks1,ke1 
           i_in  = int(roots(k,1)/dr)+1
           i_out  = int(roots(k,2)/dr)+1
           if (nfunc.ne.2) then
              wav_sum(i,j)=(2.d0*dble(k)+1.d0)*bslj(k,i) &
                   & *lpoly(k,j)*dcmplx(Im**(k))
           else if (nfunc==2) then
              !
              exp_in = exp(-Im*fint_in(k,i))
              if (k.le.lmax) then
                 exp_out = exp(Im*fint_out(k,i))*(-1.d0)**dble(k)* &
                      & exp(dimag(fint_in(k,1)))
              else
                 exp_out = exp(Im*fint_out(k,i))*(-1.d0)**dble(k)
                 if (dble(fint_out(k,i))==0.d0) then
                    exp_out = 0.d0
                    fbar = exp(-Im*fint_in(k,i_out)) - &
                         & exp(Im*fint_out(k,i_out))*(-1.d0)**dble(k)
                    !
                    exp_in = fbar * & 
                         & exp(dimag(fint_in(k,i)))/exp(dimag(fint_in(k,i_out)))
                 end if
              end if
              !
 !             wav_sum(i,j)= wav_sum(i,j) + dcmplx(-Im)**dble(k)*( exp_in & 
 !                  &  - exp_out*exp(-Im*kr0*roots(k,2)) ) &
 !                  & /(2.d0*Im*kr0*rs*r)
 !             wav_sum2(i,j) = wav_sum2(i,j) + wav_sum(i,j)*r
              !
              !
              wav_sum(i,j)=wav_sum(i,j) + (2.d0*dble(k)+1.d0)*lpoly(k,j) * &
                   & ( exp_in - exp_out*exp(-Im*kr0*roots(k,2)) ) &
                   & /(2.d0*Im*kr0*rs*r)
              wav_sum2(i,j)=wav_sum2(i,j) + (2.d0*dble(k)+1.d0)*lpoly(k,j) * &
                   & ( exp_in - exp_out*exp(-Im*kr0*roots(k,2)) ) &
                   & /(2.d0*Im*kr0*rs)
              !
              !
 !             wav_sum(i,j) = wav_sum(i,j) + exp_in
 !             wav_sum2(i,j) = wav_sum2(i,j) + exp_out
              !
           end if
           ! Error checking of wavfn - any NaN values?
           if (abs(dble(wav_sum(i,j))).ge.0.d0) then
              continue
           else
              write(*,*) 'Error: Wavefunction value',i,j
              write(*,*) wav_sum(i,j)
              write(*,*) fint_in(k,i), fint_out(k,i)
              read(*,*)
           endif
        end do lsum
     end do plsum
  end do fsum
  !
  ! Match values of wavefunction at boundary
!  do i = 1, NN
!     if ((i.gt.i_in).and.(i.lt.i_out)) then
!        wav_sum(i,:) = dble(wav_sum(i,:)) * &
!             & dble(wav_sum(i_out,:))/dble(wav_sum(i_out-1,:))
!     end if
!     if (i.lt.i_in) then
!        wav_sum(i,:) = wav_sum(i,:) * &
!             & dble(wav_sum(i_in,:))/dble(wav_sum(i_in-1,:))
!     end if
!  end do
  !
  return
end subroutine plane_wave
!
!
!
subroutine integrate(poly, rend, roots, fint, rstart)
  use param
  implicit none
  logical, intent(in) :: poly
  double precision, intent(in) :: rend(0:lend), roots(0:lend,3)
  double complex, intent(out) :: fint(0:lend,0:NN)
  double precision, intent(inout), optional :: rstart(0:lend)
  !
  integer :: i, j, k, irr, jr
  double precision :: r, rr, dr, drr, r1, sign
  double complex :: fdis, fval_tp
  double complex, save :: fturn(0:lend,1:2) = 0.d0 !! value at tp
  logical :: assgn, fsave
  !
  ! Initialization ...
  assgn=.true.;
  fint(:,:) = 0.d0
  dr = rend(lend)/NN
  drr = rend(lend)/Ndisc
  !
  !
  ! Integrating for function for each partial wave k
  ang:do k=ks1,ke1
     fsave=.true.
     if (poly) then
        sign = 1.d0
        j  = int(rstart(k)/dr)-1
        r = dble(j)*dr
        jr = j
        irr = int(rstart(k)/drr)
        rr = dble(irr)*drr
     else
        sign = -1.d0
        fint(k,int(dble(fturn(k,2)))) = dble(fturn(k,1))
        fint(k,1+int(dble(fturn(k,2)))) = dble(fturn(k,1))
!write(*,*) k
!write(*,*) fturn(k,2),fturn(k,1)
!read(*,*)
        if (k < lmax) then
           j=2
           r=2.d0*dr
           irr=Ndisc/NN
           rr = dble(irr)*drr
!           j=1
!           r=dr
!           irr=0;
!           rr = 0.d0
        else
           j  = int(roots(k,2)/dr) + 2
           r = dble(j)*dr
           jr = j
           irr = (j-1)*Ndisc/NN
           rr = dble(irr)*drr
        endif
     end if
     !
     !
     inn:do i=0,Ndisc-1
        !
        if ((i.ne.0).and.(modulo(i,Ndisc/NN)==0)) then
           if (poly) then
              j=j-1
              if ( (j<0).or.(i>irr) ) then
                 cycle inn
              endif
              if (assgn) fint(k,j) = fint(k,j+1)
           else
              j=j+1
              if (assgn) fint(k,j) = fint(k,j-1)
           end if
           r = r - dr*sign
        end if
        !
        if (poly) then
           if ((k.ge.lmax).and.(j < (int(roots(k,2)/dr)+1) )) then
              if (fsave) then
                 fturn(k,1) = fint(k,j+1); fturn(k,2) = j+1;
                 fsave=.false.
              end if
           end if
           !
           if ((k.lt.lmax).and.(j == 0)) then
              if (fsave) then
                 fturn(k,1) = fint(k,j+1); fturn(k,2) = j+1;
                 fsave=.false.
              end if
           end if
           !
           if ( (j<0).or.(i>irr) ) then
              assgn = .false.
              cycle inn
           else
              assgn = .true.
           end if
        else
           if ( (j>NN-1).or.(i>(Ndisc-irr)).or.(r>rend(k)) ) then
              assgn = .false.
              cycle inn
           else
              assgn = .true.
           endif
        endif
        !
        rr = rr - drr*sign
        !
        r1 = rr + drr*sign
        fdis = rs*(-sign*drr)/2.d0 * &
             & ( mass*sqrt(dcmplx(fr_2(k,rr))) + &
             & mass*sqrt(dcmplx(fr_2(k,r1))) )
        !
        fint(k,j) = fint(k,j) + fdis
        !
     end do inn
  end do ang
  !
  do k=ks1,ke1
     do i=j,NN
        if (abs(fint(k,j)).le.1.d-12) fint(k,j) = 0.d0
     end do
  end do
  !
  fint(:,0) = 0.d0
  !
  return
end subroutine integrate
