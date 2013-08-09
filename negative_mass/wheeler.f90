module common
implicit none
  ! Constants
  double precision, parameter :: c=137.d0,G=2.41d-43,m0=9.10938d-31
  double precision, parameter :: pi = 3.14159265358979325d0
  double complex, parameter   :: Im=(0.d0,1.d0)
  !
  ! Paramters
  double precision, parameter :: m    = 2.d-6
!  double precision, parameter :: Mass = 1.d50
  double precision, parameter :: Msun = 1.9891d30/m0
  double precision, parameter :: Mass = 1.d-8*Msun
  double precision, parameter :: Es   = m*c**2
  double precision :: l    = 600.d0 !4.6d0*2.5d0
  !
  double precision, parameter  :: E    = 1.0001d0*Es !!!! 1.5 !
  double precision :: GM   = G*Mass/c**2
  double precision :: rs = 2.d0*G*Mass/c**2
end module common
!
!
!
program tunneling
  use common
  implicit none
  !
  interface
     subroutine real_Ekin(N,xval,yval,nend)
       use common
       implicit none
       integer, intent(in) :: N
       integer, intent(out) :: nend
       double precision, allocatable :: xval(:),yval(:)
     end subroutine real_Ekin
  end interface
  !
  integer :: i,j,k,kk, loop,nend,nE0,ll,lmax,inout,kr=1,xmod,seed(1)
  integer*8 :: N, N2,end,NN
  double precision :: x,y,dx,dy,xmax,ymax,rturn(2),var1,var2, &
       & init,dv1,dv2,frc1,frc2,Ekpt(2),theta,lpoly,dth,xx,yy,teps
  double complex :: pr,pr2,dr,r,mod1,mod2,prob(2),expon,wavefn
  double precision, allocatable :: xE0(:),yE0(:),rr(:),tt(:), &
       & wwr(:),wwi(:), wfsq(:), rand(:)
  double complex, allocatable :: pt(:,:),prabs(:,:),Matrix(:,:,:), &
       & msum(:,:),mvec(:),ekvec(:), Mmin(:)
  double precision, external :: modsq
  !
  write(*,*) 'Converting data from spherical to cartesian coordinates...'
  !
  !
!goto 101
!  NN = (2000-1)**2
!  NN = 999000
!  NN = 1999000
NN = 3998000
! NN = 7998000
! NN = 7996000
! NN = 15996000
  allocate ( rr(NN),tt(NN),wwr(NN),wwi(NN),wfsq(NN) )
allocate( rand(NN+1) )
!
open(unit=20,file='random.txt',status='old',action='read')
read(20,*) seed(1)
close(20)
!
CALL RANDOM_SEED
CALL RANDOM_SEED (SIZE=kr)
CALL RANDOM_SEED (PUT=seed(1:kr))
CALL RANDOM_NUMBER(rand)

  in_out:do inout=1,3
!  in_out:do inout=3,4 !!!!!!!!!!
write(*,*) '#', inout
     !
!     if (inout==1) open(11,file="wvi_std.txt",status='old',action='read')
!     if (inout==2) open(11,file="wvo_std.txt",status='old',action='read')
!     if (inout==3) open(11,file="wvs_std.txt",status='old',action='read')
     if (inout==1) open(11,file="wav_in.txt",status='old',action='read')
     if (inout==2) open(11,file="wav_out.txt",status='old',action='read')
     if (inout==3) open(11,file="wav_sum.txt",status='old',action='read')
     if (inout==4) open(11,file="wav_sum2.txt",status='old',action='read')
     do i=1,NN
!write(*,*) inout,i
!     if (inout==3) then
        read(11,*) rr(i),tt(i),wwr(i),wwi(i),wfsq(i)
!     else
!        read(11,*) rr(i),tt(i),wfsq(i)
!     end if
  end do
  close (11)
  !
!  if (inout==1) open(12,file="nwi_std.txt",status='replace',action='write')
!  if (inout==2) open(12,file="nwo_std.txt",status='replace',action='write')
!  if (inout==3) open(12,file="nws_std.txt",status='replace',action='write')
!  if (inout==3) open(15,file="absw_std.txt",status='replace',action='write')
  if (inout==1) open(12,file="nwav_in.txt",status='replace',action='write')
  if (inout==2) open(12,file="nwav_out.txt",status='replace',action='write')
  if (inout==3) open(12,file="nwav_sum.txt",status='replace',action='write')
  if (inout==3) open(15,file="absw.txt",status='replace',action='write')
  if (inout==4) open(12,file="nwav_sum2.txt",status='replace',action='write')
  if (inout==4) open(15,file="absw2.txt",status='replace',action='write')
  !
  !
  teps = 0.d0
  wrt:do i=1,NN
     !
!     if ((tt(i).le.2.*pi).and.(abs(rr(i)).ne.0.d0)) then
!        if ( (tt(i).gt.pi-teps).or.(tt(i).lt.teps) ) cycle wrt
!        xx = rr(i)*cos(tt(i))
!        yy = rr(i)*sin(tt(i))
!if ( (xx==0.d0).or.(yy==0.d0) ) cycle wrt

!if (tt(i).gt.pi) cycle wrt !!!!!!!!!
     !
!if (modulo(i,37).ne.1) cycle wrt
!xmod=int(dble(rand(i)+2.d0)*2.d0) !20.d0
!if ((modulo(i,xmod).ne.1).and.(abs(tt(i)-pi).gt.pi/8.d0)) cycle wrt
!
        xx = rr(i)*sin(tt(i))
        yy = -rr(i)*cos(tt(i))

if ( (xx.eq.0.d0).and.(yy.eq.0.d0) ) cycle wrt

!if ( yy < 0.d0 ) cycle wrt
!
!if ( yy < 0.d0 ) cycle wrt
!if ( yy > 2.d0 ) cycle wrt
!if ( yy < 94.d0 ) cycle wrt
if ( abs(yy) > 10.d0 ) cycle wrt
!
if ( abs(xx) > 10.d0 ) cycle wrt
!if ( abs(xx) > 0.2d0 ) cycle wrt
!
!        write(12,fmt='(4F30.14)') xx,yy,log10(wwr(i)+1.d0),log10(wwi(i)+1.d0)
        write(12,fmt='(4F30.14)') xx,yy,wwr(i),wwi(i)
        write(15,fmt='(4F30.14)') xx,yy,wfsq(i),log10(wfsq(i)+1.d0)
!        write(12,fmt='(4F30.14)') xx,yy,wwr(i),wwi(i) !,wfsq(i)
!     else
!        write(*,*) i,rr(i),tt(i)
!        read(*,*)
!     end if
  end do wrt
  !
  close(12)
!
end do in_out
open(unit=21,file='random.txt',status='replace',action='write')
write(21,*) int(rand(NN+1)*1.d6)
close(21)
deallocate( rand )
!
deallocate( rr, tt, wwr, wwi, wfsq )
write(*,*) '...done'
stop
101 continue
  !
!!!!!!!!!!!!!! loop over parameters - for plots !!!!!!
N2=1
 !var1 = GM
 !init = l
!frc1=.8d0;frc2=0.5d0;
 !dv1 = (frc1*var1 - frc2*var1)/N2
!dv2 = (frc1*init - frc2*init)/N2
 !
 !write(*,*) var1,init
 !
 !var1 = frc2*var1
 !do k=1,N2
 !  var1 = var1 + dv1
 !  GM = var1
!  var2 = init
!  var2 = frc2*var2
 !  write(*,*) k
  do kk=1,N2
!     var2 = var2 + dv2
!     l = var2
  write(*,*) l
!!!!!!!!!!!!!! -------------------------------- !!!!!
  !
  !
  int_out: do i=1,2
  !
  r=0.;pr=0.;pr2=0.;dr=0.;
  if (i==1) then
     r=1.d0
     pr = m*(E**2.d0/Es - (1.d0-1.d0/r)*(Es+l**2.d0/(m*(r*rs)**2.d0)))
  else
     r=rturn(1)
     pr = m*(E**2.d0/Es - (1.d0-1.d0/r)*(Es+l**2.d0/(m*(r*rs)**2.d0)))
     pr2=pr
     dr=1.d0/1000.d0
  incr:do while(.true.)
        r = r + dr
        pr = m*(E**2.d0/Es - (1.d0-1.d0/r)*(Es+l**2.d0/(m*(r*rs)**2.d0)))
        if (modsq(pr)>modsq(pr2)) then
           pr2=pr
           cycle incr
        else
           exit incr
        end if
     enddo incr
  endif
  !
  ! Finding the turning point
  !
  pr2 = pr
  dr = 1.d0/10000.d0
  loop = 0
  turnpt:do while(.true.)
!     if (real(r) > 10.d0) then
!        write(*,*) "No turning point before r = 10 rs"
!        exit int_out
!     end if
     loop = loop + 1
        r = r + dr
        pr = m*(E**2.d0/Es - (1.d0-1.d0/r)*(Es+l**2.d0/(m*(r*rs)**2.d0)))
     if (modsq(pr)<modsq(pr2)) then
        pr2=pr
        cycle turnpt
     else
        r = r - 2.d0*dr
        pr = m*(E**2.d0/Es - (1.d0-1.d0/r)*(Es+l**2.d0/(m*(r*rs)**2.d0)))
        pr2 = pr
        dr = dr/2.d0
     end if
     if (sqrt(modsq(pr))<1.d-20) exit turnpt
     if (loop>1000000) then
        write(*,*) 'Cannnot converge to find turning point'
        write(*,*) 'Error: ',sqrt(modsq(pr))
        write(*,*) 'Value: ',r
        exit turnpt
     endif
  enddo turnpt
  ! 
  rturn(i) = dble(r)
  write(*,*) 'Turning point at ',rturn(i)*rs, 'schwarz_r ', rturn(i)
  !
enddo int_out
  !
!write(14,fmt='(4F24.16)') var1,var2,modsq(pr),rturn
!goto 100
  !
  !
! call energy_integral(rturn)
!stop
  !
  ! Discretize complex space
!  N = 10
!  write(*,*) '# points ',N
!  allocate (pt(0:N,-N:N))
!  xmax = rturn(1)
!  ymax = rturn(1) ! (largest imaginary value in grid) !
!  dx = xmax/N
!  dy = ymax/N
!  do i=0,N
!     do j=0,N
!        pt(i,j)  = cmplx(i*dx,j*dy)
!        pt(i,-j) = cmplx(i*dx,-j*dy)
!     end do
!  end do
  !
  !
  ! Find region where Ekin=0 !
!  nE0 = 1000
!  call real_Ekin(nE0,xE0,yE0,nend)
  !
  !
!100 continue
  ! Calculate tunneling integral along a straight line path to the
  ! points where E_kin=0.
  !
!  Ekpt(1) = rturn(1)
!  Ekpt(2) = 0.0d0
!rturn(2) = 23.30852381126683d0
!  if (rturn(1)>Ekpt(1)) Ekpt(1) = rturn(1)

!  write(*,*) 'Tunneling from', rturn(2), 'to',Ekpt(1),Ekpt(2)
  !
write(*,*) l
!lmax=1000
!l=100.d0
!theta = pi/10.d0
!call legendre(lmax,theta,lpoly)
!stop
  !
l = 305.d0 !!!!!!!!!!!!!!!!!!!!!!
  !
open(unit=20,file='wavtemp.txt',status='replace',action='write')
  lmax = int(l)
  NN = 50
  r = 0.d0
  dr = rturn(2)/NN !* 2.d0 ! for outgoing wave
  do i=0,NN-1
write(*,*) 'i=',i
     r = r + dr
     theta= 0.d0
     dth = pi/NN
     do j=0,NN-1
!write(*,*) 'j=',j
        theta = theta + dth
        wavefn = 0.d0
        do ll=0,lmax
           l=dble(ll)
           call line_integral(r,rturn,expon)
!           rturn(:)=0.d0                      ! for outgoing wave
!           call line_integral(r,rturn,expon)  !        "
!write(*,*) r,l
!write(*,*) expon
!read(*,*)
           call legendre(lmax,theta,lpoly)
           wavefn = wavefn + (2.d0*dble(l)+1.d0)*lpoly*expon
        enddo
!        write(20,fmt='(4F30.8)') real(r),theta,wavefn
        write(20,fmt='(4F30.15)') real(r),theta,wavefn
     end do
  end do
close(20)
  !
!  write(13,fmt='(4F20.10)') var1,var2,prob(1)
!  write(*,fmt='(2F20.10)') l,prob(1)
!enddo

exit

enddo
write(*,*) 'stop'
stop
!!!!!!!!!!!!------------------------------------------!!!!!!!!!!!!!

  !
  !
  ! Evaluate the momentum on the grid
  allocate(prabs(0:N,-N:N))
  do i=0,N
     do j=-N,N
        r = pt(i,j)
        mod1=modsq(E*r/(r-GM))
        mod2=modsq(l*c/r)
        prabs(i,j) = (1.d0/c)*sqrt(mod1 - mod2 - Es**2)
!write(10,fmt='(2I6,2F20.15)') i,j,prabs(i,j)
     enddo
  enddo
!stop
  !
  ! 'Adjust values of prabs'
  prabs(0,0) = cmplx(0.d0,0.d0)
!  do i=0,N
!     do j=-N,N
!        if(real(prabs(i,j)).ne.0.) prabs(i,j) = cmplx(0.d0,0.d0)
!     enddo
!  enddo
  !
  !
  ! Construct the tunneling matrices (trapezoidal rule)
  allocate( Matrix(0:N-1,-N:N,-N:N), Mmin(0:N-1) )
write(*,*) dx
read(*,*)
  do k=0,N-1
     do i=-N,N
        do j=-N,N
           Matrix(k,i,j) = dx/2.d0 * (prabs(k,i)+prabs(k+1,j))
write(11,fmt='(3I5,2F20.16)') k,i,j,Matrix(k,i,j)
        enddo
!        call normalize(Matrix(k,i,:))
     enddo
  enddo
stop
  !
  ! Subtract off smallest element in column from others
  ! S_ = S - S_min
  Mmin(:) = 1.d99
  do k=0,N-1
     do i=-N,N
        do j=-N,N
           if (real(Matrix(k,i,j))==0.) then
              if (abs(Matrix(k,i,j))<abs(Mmin(k))) Mmin(k) = Matrix(k,i,j)
           endif
        enddo
! write(*,*) Mmin(k)
     enddo
     do i=-N,N
        do j=-N,N
           if (real(Matrix(k,i,j))==0.) Matrix(k,i,j) = Matrix(k,i,j) - Mmin(k)
        enddo
     enddo
  enddo
  !
  ! Combine matrices to points where E_kin=0
  ! ... to multiply the matrices we must exponentiate its elements
  do k=0,N-1
     Mmin(k) = exp(2.d0*Mmin(k))
 write(14,fmt='(2I5,2F20.16)') k,i,Mmin(k)
     do i=-N,N
        do j=-N,N
           if (real(Matrix(k,i,j))==0.) then
              Matrix(k,i,j) = exp(-2.d0*abs(Matrix(k,i,j)))
           else
              Matrix(k,i,j) = exp(-2.d0*Im*abs(Matrix(k,i,j)))
           endif
 write(12,fmt='(3I5,2F20.16)') k,i,j,Matrix(k,i,j)
        enddo
     enddo
  enddo
stop
  !
  ! - constant region - !
  allocate ( msum(-N:N,-N:N), mvec(-N:N) )
  msum(:,:) = matmul(Matrix(end+1,:,:),Matrix(end+2,:,:))
!!!!  
  do i=end+2,N-2
     msum(:,:) = matmul(msum(:,:),Matrix(i+1,:,:))
     write(*,*) msum(10,10)
read(*,*)
  end do

!do i=-N,N
!do j=-N,N
!write(11,fmt='(2I5,3F20.10)') i,j,dble(msum(i,j)),dimag(msum(i,j))
!enddo
!enddo
!stop

  mvec(:) = 0.d0
  do i=-N,N
     mvec(:) = mvec(:) + msum(:,i)*Matrix(N-1,i,0)
 end do
write(11,fmt='(2F20.12)') mvec(:)
  !
  ! - E_kin=0 region - !
  allocate ( ekvec(0:end) )
  ekvec(:) = 0.d0
  do i=end,0
     do k=-N,N
!        ekvec(i) = Matrix(i,eky(i),k)*mvec(k) ---eky=y-value where Ekin=0 !
     enddo
     mvec(:) = matmul(mvec(:),Matrix(i,:,:))
  enddo
  !
  write(*,*) ekvec(:)
  !
  write(*,*) 'end'
  !
end program tunneling
!
!
!subroutine normalize(N,MM)
!  implicit none
!  integer, intent(in) :: N
!  double complex, intent(inout) :: MM(-N,N)
!  integer :: i
!  double precision :: sum, norm
  !
!  sum = 0.d0
!  do i=-N,N
!     sum = sum + MM(i)
!  end do
!  norm = 1.d0/sqrt(sum) ! ??? !
  !
!end subroutine normalize
!
!
function modsq(x)
  implicit none
  double complex, intent(in) :: x
  double precision :: modsq
  !
  modsq = dble(x*conjg(x))
  !
end function modsq
!
!
!
subroutine line_integral(Ekpt,rturn,wavf)
  use common
  implicit none
  double precision :: Ekpt(2),rturn(2)
  double complex, intent(out)  :: wavf
  double precision, external   :: modsq
  !
  integer                      :: i,j,N
  double precision             :: x,y,xp,yp,dx,eps,repr,impr
  double complex               :: r,tint,mod1
  double complex, allocatable  :: prcmplx(:)
  !
  N = 10000 ! >=100000 
  allocate (prcmplx(0:N))
  !
  eps=1.d-12
  ! Equation for line : y = sqrt(x*(GM-x)) !
  ! 
  xp = Ekpt(1)
!  yp = Ekpt(2)
  !
  x  = xp
  dx = (rturn(2)-xp)/N
  do i=0,N
     x = x + dx
!     y = yp*(x-rturn)/(xp-rturn)
     r = dble(x) !+ Im*y
     !
!     mod1=sqrt(dcmplx(modsq((1.d0-1.d0/r)*(Es+l**2.d0/(m*(r*rs)**2.d0)))))
!     prcmplx(i) = sqrt(m)*sqrt(dcmplx(E**2.d0/Es - mod1))
!     prcmplx(i) = 1.d0/prcmplx(i)
     prcmplx(i) = (dcmplx(E**2.d0/(m*Es) - (1.d0-1.d0/r)*&
          &(c**2.d0 + l**2.d0/(m*r*rs)**2.d0)))**(-0.5d0)
     !
     repr = dble(prcmplx(i))
     impr = dimag(prcmplx(i))
     if ( abs(repr) < eps ) prcmplx(i)=dcmplx(0.d0,impr)
     if ( abs(impr) < eps ) prcmplx(i)=dcmplx(repr,0.d0)
     !
!     write(17,fmt='(I8,4F30.15)') i,dble(r),prcmplx(i),real(abs(prcmplx(i)))
     !
  enddo
  !
  ! Integration along line
  tint = 0.d0
  do i=0,N-1
     tint = tint + dx*( prcmplx(i) + &
          & prcmplx(i+1) )/2.d0
  enddo
  tint = tint * rs !! since integrated in units of rs !
!write(*,*) tint
  !
  deallocate( prcmplx )
  !
  ! Tunneling probability
  wavf = exp(-Im*Es*tint)
  !
end subroutine line_integral
!
!
!
subroutine real_Ekin(N,xval,yval,nend)
  use common
  implicit none
  integer, intent(in) :: N
  integer, intent(out) :: nend
  double precision, allocatable :: xval(:),yval(:)
  !
  integer :: i,j
  integer*8 :: loop
  double precision :: x,y,dx,dy,eps,sg
  double precision, allocatable :: err(:), err2(:)
  double complex   :: prsq,z,Ekinsq,Ekinsq2
!  double complex, allocatable :: err2(:)
  !
  !
  !
!  Ekinsq = (c**2.*l**2.*x*(x**2.-3.d0*y**2.)+rs**2.*(x**2.+y**2.)**2. &
!       & *(Es**2.*x+E**2.*(x**2.+y**2.)))/(rs**2.*(x**2+y**2.)**3.)
  !
  sg = -1.d0
  allocate( xval(N),yval(N) )
  allocate( err(N), err2(N) )
!  xval(:)=0.;yval(:)=0.; err(:)=0.
  eps = 1.d-12
  !
  !
  x=0.d0
  dx = 5.d0/N
  xloop:do i=1,N
     x= x + dx
!     y = 0.d0
     y = 4.d0
     !
!write(*,*) i,x
     !
     z = (x + Im*y)*rs
     prsq = m*(E**2.d0/Es - (1.d0-rs/z)*(Es + l**2.d0/(m*z**2.d0)))
     Ekinsq = Es**2.d0 + prsq*c**2.d0 + (l*c/z)**2.d0
     !
     dy = 1.d0/1.d5
     Ekinsq2 = Ekinsq
     loop = 0
     itr:do while(.true.)
        loop = loop + 1
        y = y + dy*sg
        z = (x + Im*y)*rs
        prsq = m*(E**2.d0/Es - (1.d0-rs/z)*(Es+l**2.d0/(m*z**2.d0)))
        Ekinsq = Es**2.d0 + prsq*c**2.d0 + (l*c/z)**2.d0
        if (loop==1000000) then
!           write(*,*) 'Ended search points where Ekin=0'
!           write(*,*) i,x
           nend = i
!           return
           exit xloop
        endif
        if (real(Ekinsq)*real(Ekinsq2)<0.d0) then
           loop = 0
           y = y - 2.d0*dy*sg
           dy = dy/2.d0
           z = (x + Im*y)*rs
           prsq = m*(E**2.d0/Es - (1.d0-rs/z)*(Es+l**2.d0/(m*z**2.d0)))
           Ekinsq2 = Es**2.d0 + prsq*c**2.d0 + (l*c/z)**2.d0        
           yl:do while(.true.)
              loop = loop + 1
              y = y + dy*sg
              z = (x + Im*y)*rs
              prsq = m*(E**2.d0/Es - (1.d0-rs/z)*(Es+l**2.d0/(m*z**2.d0)))
              Ekinsq = Es**2.d0 + prsq*c**2.d0 + (l*c/z)**2.d0           
              if (Abs(real(Ekinsq)-real(Ekinsq2))<eps) then
                 xval(i) = x
                 yval(i) = y
                 err(i) = Abs(real(Ekinsq))
                 err2(i) = sqrt(real(Ekinsq*conjg(Ekinsq)))
!                 write(15,fmt='(I5,3F25.15)') i,x,y,Abs(real(Ekinsq))
                 exit itr
              endif
              !
!              if (modulo(loop,10000)==0) then
!                 write(*,*) loop
!                 write(*,*) x,dx,y,dy
!                 write(*,*) Ekinsq,Ekinsq2
!              endif
              !
              if (Abs(real(Ekinsq))<Abs(real(Ekinsq2))) then
                 Ekinsq2=Ekinsq
                 cycle
              else
                 y = y - 2.d0*dy*sg
                 dy = dy/2.d0
                 z = (x + Im*y)*rs
                 prsq = m*(1.d0-rs/z)*(E**2.d0/Es - &
                      & (1.d0-rs/z)*(Es+l**2.d0/(m*z**2.d0)))
                 Ekinsq = Es**2.d0 + prsq*c**2.d0 + (l*c/z)**2.d0           
                 Ekinsq2 = Ekinsq
                 cycle
              endif
           enddo yl
           !
        endif
        !
     enddo itr
     !
     write(16,fmt='(I5,4F25.15)') i,xval(i),yval(i),err(i),err2(i)
     !
  enddo xloop
  !
  deallocate( err )
  return
  !
end subroutine real_Ekin
!
!
!
subroutine energy_integral(rturn)
  use common
  implicit none
  double precision, intent(in) :: rturn(2)
  !
  integer*8 :: i,j,k,k1,itr,itr2
  double precision :: r,dr,dr2,rend,outer,inner
  double precision, allocatable :: rr(:)
  double complex, allocatable :: fn(:), eint(:),mag(:), &
       & phase(:),wavefn(:),prob(:)
  double precision, external :: modsq
  !
  itr=100000
  itr2=10000
  !
!  outer=3.d0*rturn(2)
!  inner=0.5d0*rturn(1)
  outer=10.d0
  inner=0.d0
write(*,*) inner,outer
  !
  allocate( rr(0:itr),fn(0:itr),eint(1:itr2), &
       & mag(1:itr2),phase(1:itr2),wavefn(1:itr2),prob(1:itr2) )
  !
  r = outer
  dr = (outer-inner)/itr
  do i=0,itr
     r = outer - i*dr
     rr(i) = r
     fn(i) = 1.d0/sqrt(dcmplx( E**2./(m*Es) - &
          & (1.d0-1.d0/r)*(c**2.+l**2./(m*(rs*r))**2.)) )
!  write(20,fmt='(3F20.10)') rr(i),fn(i)
  end do
  !
  !
  rend=0
  dr2 = (outer-inner)/itr2
  do i=1,itr2
     rend = outer - i*dr2
     k1=0
     inn:do k=1,itr
        if (rr(k+1).lt.rend) then
           k1=k-1
           exit inn
        endif
     end do inn
     eint(i) = 0.d0
     do j=0,k1-1
        eint(i) = eint(i) + (dr)*( fn(j) + fn(j+1) )/2.d0
     end do
     eint(i) = eint(i)*rs
!  write(21,fmt='(3F20.10)') rend,eint(i)
     mag(i) = exp(Es*dimag(eint(i)))
     phase(i) = exp(-Im*Es*dble(eint(i)))
     wavefn(i) = mag(i)*phase(i)
     prob(i) = modsq(wavefn(i))
!  write(22,fmt='(3F20.16)') rend,mag(i)
!  write(23,fmt='(3F20.16)') rend,phase(i)
!  write(24,fmt='(3F20.16)') rend,wavefn(i)
!  write(25,fmt='(3F20.16)') rend,prob(i)
  end do
  !
  stop
  !
  deallocate ( rr, fn, eint )
  return
  !
end subroutine energy_integral
!
!
subroutine legendre(lmax,theta,lpoly)
  use common
  implicit none
  integer :: i,lmax
  double precision :: theta,lpoly,lp(0:int(lmax)),lpoly2
  !
  lp(0)=1.d0
  lp(1)=cos(theta)
  lp(2)=0.5d0*(3.d0*cos(theta)**2.d0-1.d0)
  !
  if (l.ge.100) then
     lpoly = sqrt(2.d0/(pi*(dble(l)+0.5d0)*sin(theta)))* &
          & cos((dble(l)+0.5d0)*theta - pi/4.d0)
     return
  end if
  !
  if (l.ge.3) then
     do i=2,l-1
        lp(i+1) = ( (2.d0*dble(i)+1.d0)*cos(theta)*lp(i)-dble(i)*lp(i-1) ) &
             & / (dble(i)+1.d0)
     end do
  endif
  lpoly = lp(int(l))
  !
  return
  !
end subroutine legendre
!
