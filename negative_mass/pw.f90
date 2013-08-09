program pw
implicit none
integer :: i,j,k,N1,N2,lmax
double precision :: r,theta,pi,dr,dth
double complex :: Im=dcmplx(0.d0,1.d0)
double precision, allocatable :: x(:,:),y(:,:), lpoly(:,:),besselj(:,:)
double complex, allocatable :: pw0(:,:),pw1(:,:)!,pw2(:,:),pw3(:,:)
!
pi = 4.d0*atan(1.d0)
!
N1=8000
N2=8000
lmax=100
allocate(lpoly(0:lmax,N2),besselj(0:lmax,N1),x(N1,N2),y(N1,N2), &
     & pw0(N1,N2),pw1(N1,N2) ) !,pw2(N1,N2),pw3(N1,N2))
!
dr = 8000.d0/N1
dth = 2.*pi/N2
!
do k=0,lmax
   theta = 0.d0
   do j=1,N2
      theta = theta + dth
      call legendre(k,theta,lpoly(k,j))
   enddo
   r = 0.d0
   do i=1,N1
      r = r + dr
      call sph_bessel(k,0,r,besselj(k,i))
   end do
end do
!
r=0.d0
pw1(:,:) = 0.d0
do i=1,N1
   if (modulo(i,500)==0) write(*,*) i,'of',N1
   r = r + dr
   theta = 0.d0
   do j=1,N2
      theta = theta + dth
      x(i,j) = r*sin(theta)
      y(i,j) = -r*cos(theta)
      if ((y(i,j).ge.7980.d0).and.(abs(x(i,j)).le.50.d0)) then
         pw0(i,j) = exp(Im*r*cos(theta))
         do k=0,lmax
            pw1(i,j) = pw1(i,j) + Im**(k)*(2.d0*dble(k)+1.d0)* &
                 & besselj(k,i)*lpoly(k,j)
!            pw2(i,j) = pw2(i,j) + Im**(k)*(2.d0*dble(k)+1.d0)* &
!                 & (sin(r-dble(k)*pi/2.d0)/r)*lpoly(k,j)
!            pw3(i,j) = pw3(i,j) + (2.d0*dble(k)+1.d0)* &
!                 & ((exp(Im*r)-exp(-Im*r))/(2.d0*Im*r))*lpoly(k,j)
  !               & ((exp(Im*r)-(-1.d0)**(k)*exp(-Im*r))/(2.d0*Im*r))*lpoly(k,j)
         enddo
         
         write(10,fmt='(4F30.15)') x(i,j),y(i,j),dble(pw0(i,j)),abs(pw0(i,j))
         write(11,fmt='(4F30.15)') x(i,j),y(i,j),dble(pw1(i,j)),abs(pw1(i,j))
!         write(12,fmt='(4F30.15)') x(i,j),y(i,j),dble(pw2(i,j)),abs(pw2(i,j))
!         write(13,fmt='(4F30.15)') x(i,j),y(i,j),dble(pw3(i,j)),abs(pw3(i,j))
      end if
   enddo
enddo
!
!
stop
!
end program pw
!
!
subroutine legendre(ll,theta,lpoly)
  implicit none
  integer :: i
  integer, intent(in) :: ll
  double precision, intent(in) :: theta
  double precision, intent(out) :: lpoly
  double precision :: lp(0:500),ang
  double precision :: pi
  logical :: flip
  !
  pi = 4.d0*atan(1.d0)
  !
  flip = .false.
  lp(0)=1.d0
  lp(1)=cos(theta)
  lp(2)=0.5d0*(3.d0*cos(theta)**(2.d0)-1.d0)
  !
  if (ll.gt.500) then
     if ( abs(sin(theta)) < 1.d-6) then
        lpoly = 1.d0
        if ( (modulo(ll,2)==1).and.(abs(theta-pi)<1.d-6) ) lpoly = -1.d0
        return
     endif
     ang = theta
     if (ang > pi) then
        ang = ang - pi
        if (modulo(ll,2)==1) flip = .true.
     endif
     !
     lpoly = sin ( (dble(ll)+0.5d0)*ang+pi/4.d0 )/ &
          & sqrt((dble(ll)+1.d0)*pi/2.d0 * sin(ang))
     if (flip) lpoly = -1.d0*lpoly
     return
  end if
  !
  if (ll.ge.3) then
     do i=2,ll-1
        lp(i+1) = ( (2.d0*dble(i)+1.d0)*cos(theta)*lp(i)-dble(i)*lp(i-1) ) &
             & / (dble(i)+1.d0)
     end do
  endif
  lpoly = lp(int(ll))
  !
  !
  return
end subroutine legendre
!
!
subroutine sph_bessel(ll,nfunc,z,besselj)
  implicit none
  integer, intent(in) :: ll, nfunc
  double precision, intent(in)  :: z
  double precision, intent(out) :: besselj
  integer :: i
  double precision :: sb(0:10000)
  double precision  :: pi
  !
  pi = 4.d0*atan(1.d0)
  !
  if (nfunc==0) then
     sb(0) = sin(z)/z
     sb(1) = sin(z)/z**2.d0 - cos(z)/z
     do i=1,ll-1
        sb(i+1) = (2.d0*dble(i)+1.d0)*sb(i)/z - sb(i-1)
        if ( (ll>4).and.(z<(1.d-1+(ll-5)**1.16*0.3d0) ) ) then
           besselj = 0.d0
           return
        end if
     end do
     besselj = sb(ll)
     !
  else if (nfunc==1) then  ! asymptotic limit r->infty
     besselj = sin(z - ll*pi/2.d0)/z
     !
  else
     write(*,*) 'Error in subroutine sph_bessel'
     stop
  endif
  !
  return
end subroutine sph_bessel
!
