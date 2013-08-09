program pw
implicit none
integer :: i,j,k,N1,N2,lmax
double precision :: r,theta,pi,dr,dth,x,y,rend
double complex :: Im=dcmplx(0.d0,1.d0),pw0,pw1,pw2,pw3
double precision, allocatable :: lpoly(:,:),besselj(:,:)
logical :: wrt
!
pi = 4.d0*atan(1.d0)
!
wrt=.true.
rend=1000.d0
N1=rend*2.d0
N2=4000
lmax=1000
allocate(lpoly(0:lmax,N2),besselj(0:lmax,N1))
!
dr = rend/N1
dth = 2.*pi/N2
!
write(*,*) 'Generating Legendre Polyns...'
theta = 0.d0
do j=1,N2
   theta = theta + dth
   call legendre(lmax,theta,lpoly(0:lmax,j))
enddo
!
write(*,*) 'and Spherical Bessel Functions...' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
r = 0.d0
do i=1,N1
   r = r + dr
   call sph_bessel(lmax,0,r,besselj(0:lmax,i)) !!! problem with high orders !!!
end do
!
if (wrt) then
do k=0,lmax
write(20,*)
write(20,*) k
!do i=1,N1
!write(20,fmt='(I8,F30.15)') i, besselj(k,i)
do j=1,N2
write(20,fmt='(I8,F30.15)') j, lpoly(k,j)
enddo
write(20,*)
enddo
stop
end if
!
write(*,*) 'Producing plane waves'
r=0.d0
do i=1,N1
   if (modulo(i,500)==0) write(*,*) i,'of',N1
   r = r + dr
   theta = 0.d0
   do j=1,N2
      theta = theta + dth
      x = r*sin(theta)
      y = -r*cos(theta)
      if ((abs(y-(rend-40.d0)).le.20.d0).and.(abs(x).le.20.d0)) then
         pw0 = exp(Im*r*cos(theta))
         pw1=0.d0;pw2=0.d0;pw3=0.d0;
         do k=0,lmax
!            pw1 = pw1 + Im**(k)*(2.d0*dble(k)+1.d0)* &
!                 & besselj(k,i)*lpoly(k,j)
            pw2 = pw2 + Im**(k)*(2.d0*dble(k)+1.d0)* &
                 & (sin(r-dble(k)*pi/2.d0)/r)*lpoly(k,j)
            pw3 = pw3 + (2.d0*dble(k)+1.d0)* &
                 & ((exp(Im*r)-exp(-Im*r))/(2.d0*Im*r))*lpoly(k,j)
!                 & ((exp(Im*r)-(-1.d0)**(k)*exp(-Im*r))/(2.d0*Im*r))*lpoly(k,j)
! non-alternating !!!!!!!!!!
         enddo
         
         write(10,fmt='(5F30.15)') x,y,dble(pw0),dimag(pw0),abs(pw0)
!         write(11,fmt='(4F30.15)') x,y,dble(pw1),abs(pw1)
         write(12,fmt='(5F30.15)') x,y,dble(pw2),dimag(pw2),abs(pw2)
         write(13,fmt='(5F30.15)') x,y,dble(pw3),dimag(pw3),abs(pw3)
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
subroutine legendre(lmax,theta,lpoly)
  implicit none
  integer :: i
  integer, intent(in) :: lmax
  double precision, intent(in) :: theta
  double precision, intent(out) :: lpoly(0:lmax)
  double precision :: pi, ang
  logical :: flip
  !
  pi = 4.d0*atan(1.d0)
  !
  flip = .false.
  lpoly(0)=1.d0
  if (lmax<1) return
  lpoly(1)=cos(theta)
  if (lmax<2) return
  lpoly(2)=0.5d0*(3.d0*cos(theta)**(2.d0)-1.d0)
  if (lmax<3) return
  !
  loop:do i=2,lmax-1
     if ((i+1).lt.201) then
        lpoly(i+1) = ( (2.d0*dble(i)+1.d0)*cos(theta)*lpoly(i)- &
             & dble(i)*lpoly(i-1) ) / (dble(i)+1.d0)
     else
        if ( abs(sin(theta)) < 1.d-10) then
           lpoly(i+1) = 1.d0
           if ((modulo(i+1,2)==1).and.(abs(theta-pi)<1.d-10)) lpoly(i+1)=-1.d0
           cycle loop
        endif
        ang = theta
        if (ang > pi) then
           ang = ang - pi
           if (modulo(i+1,2)==1) flip = .true.
        endif
        !
        lpoly(i+1) = sin ( (dble(i+1)+0.5d0)*ang+pi/4.d0 )/ &
             & sqrt((dble(i+1)+1.d0)*pi/2.d0 * sin(ang))
        if (flip) then
           lpoly(i+1) = -1.d0*lpoly(i+1)
           flip = .false.
        end if
     end if
  end do loop
  !
  return
end subroutine legendre
!
!
subroutine sph_bessel(lmax,nfunc,z,besselj)
  implicit none
  integer, intent(in) :: lmax, nfunc
  double precision, intent(in)  :: z
  double precision, intent(out) :: besselj(0:lmax)
  integer :: i
  double precision  :: pi
  !
  pi = 4.d0*atan(1.d0)
  !
  besselj(0) = sin(z)/z
  besselj(1) = sin(z)/z**2.d0 - cos(z)/z
  do i=1,lmax-1
     if (nfunc==0) then
        besselj(i+1) = (2.d0*dble(i)+1.d0)*besselj(i)/z - besselj(i-1)
        if ( (i>4).and.(z<(1.d-1+(i-5)**1.16*0.3d0) ) ) then
           besselj(i+1) = 0.d0
        end if
     else if (nfunc==1) then  ! asymptotic limit r->infty
        besselj(i+1) = sin(z - (i+1)*pi/2.d0)/z
     end if
  end do
  !
  return
end subroutine sph_bessel
!
