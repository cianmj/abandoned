program grsolve
!
!Use fbauhube !!!!!!!! mass set to negative (sm,mass,l) !!!!
!
implicit none
!
integer, parameter :: s1 = 1,sm = -1
double complex, parameter :: c=137, G=2.41d-43
double complex, parameter :: M=1.d50, mass = -2.d-6, l=-1
!
integer :: i, j, err, err2, err3, num
double complex :: E, Es, r
double complex :: Wm, dE, ftnew(1:4),ftlag(1:4),ftqua(1:4)
double complex :: A(0:4),B(0:4),D(0:4),dv(1:4), mvv(1:4),pkin(4),zero(4)
double complex :: vtnew(4), vtlag(0:4),vtqua(0:3),vtbau(0:3),roots(1:4)
double complex, allocatable :: vt(:,:),vr(:,:)
double precision :: AA(5),BB(5),CC(5),DD(5),RR(5),RI(5),RR2(0:4),RI2(0:4),val(5)
!
num = 1000
allocate( vt(num,4),vr(num,4) )
!
Es = mass*c*c
E  = Abs(2.d0 * Es)
!
r = -9.
do i=1,num
   r = r + 10
   !
   Wm = -sm*G*M/r
   dE = sqrt(E**2. - (l*c/r)**2.)
   !
!
!   ft0 = A(4)*vt**4. + A(3)*vt**3. + A(2)*vt**2. + A(1)*vt + A(0)
!
   A(0) = l**2. * (Wm**2. - c**4.)**2.
   A(1) = 2.d0*r*E*l*(Wm*c**4.-Wm**3.)
   A(2) = r**2.*(E**2.*Wm**2.-c**4.*dE**2.) + c**2.*l**2.*(Wm**2.-2.d0*c**4.)
   A(3) = -2.d0*r*E*Wm*l*c**2.
   A(4) = r**2.*c**2.*dE**2.+l**2.*c**4.
   !
   !
   do j=0,4
      AA(j+1) = dble(A(j))
      BB(j+1) = dimag(A(j))
   enddo
   CC(:)=AA(:)
   BB(:)=BB(:)
   !
   B(:) = A(:)
   D(:) = A(:)
   !
!!!!!! The solutions to both the 'newton' and 'tclague' must be inverted !!!!
   call test_newton(AA,BB,err,RR,RI)
   do j=1,4
      vtnew(j) = RR(j) + dcmplx(0.d0,1.d0)*RI(j)
   end do
   call tclague(B,err2,vtlag)
   !
!   call quartic(D,vtqua(0:3)) ! Works for reasonable coefficients !
   !
!   err3 = bauhub(0,1,4,CC,DD,RR2,RI2,val)
!   do j=0,4
!      vtbau(j) = RR2(j) + dcmplx(0.d0,1.d0)*RI2(j)
!   end do
   !
   !
   roots(1:4) = 1.d0/vtnew(1:4)
   dv(1:4) = sqrt(c**2. - roots(1:4)**2.)
   vr(i,1:4) = c * ( 1 - roots(1:4)**2./c**2. - (mass**2.*(Wm**2.- &
        & c**2.*dv(1:4)**2.)**2.)/(E*Wm - s1*c*(dE**2.*dv(1:4)**2. &
        & + (l*Wm/r)**2.)**(1./2.))**2.)**(1./2.)
   !
   mvv(1:4) = mass * (1. - (roots(1:4)**2.+vr(i,1:4)**2.)/c**2.)**(-1./2.)
   !
   write(10,"(9F20.5)") dble(r),(dble(roots(j)), &
        & dimag(vt(i,j)),j=1,4)
   write(11,"(9F20.5)") dble(r),(dble(vr(i,j)), &
        & dimag(vr(i,j)),j=1,4)
   write(12,"(9F20.10)") dble(r),(dble(mvv(j)), &
        & dimag(mvv(j)),j=1,4)
   !
   pkin(1:4) = mvv(1:4)*sqrt( vr(i,1:4)**2. + roots(1:4)**2. )
   zero(1:4) = (pkin(1:4)*c)**2. + Es**2. - (E + G*M*mvv(1:4)/r)**2.
   write(13,"(9F30.10)") dble(r),(dble(zero(j)), &
        & dimag(zero(j)),j=1,4)
   !
   !
   !
   roots(1:4) = 1.d0/vtlag(1:4)
   dv(1:4) = sqrt(c**2. - roots(1:4)**2.)
   vr(i,1:4) = c * ( 1 - roots(1:4)**2./c**2. - (mass**2.*(Wm**2.- &
        & c**2.*dv(1:4)**2.)**2.)/(E*Wm - s1*c*(dE**2.*dv(1:4)**2. &
        & + (l*Wm/r)**2.)**(1./2.))**2.)**(1./2.)
   !
   mvv(1:4) = mass * (1. - (roots(1:4)**2.+vr(i,1:4)**2.)/c**2.)**(-1./2.)
   !
   write(20,"(9F20.5)") dble(r),(dble(roots(j)), &
        & dimag(roots(j)),j=1,4)
   write(21,"(9F20.5)") dble(r),(dble(vr(i,j)), &
        & dimag(vr(i,j)),j=1,4)
   write(22,"(9F30.10)") dble(r),(dble(mvv(j)), &
        & dimag(mvv(j)),j=1,4)
   !
   !
   pkin(1:4) = mvv(1:4)*sqrt( vr(i,1:4)**2. + roots(1:4)**2. )
   zero(1:4) = (pkin(1:4)*c)**2. + Es**2. - (E + G*M*mvv(1:4)/r)**2.
   write(23,"(9F20.10)") dble(r),(dble(zero(j)), &
        & dimag(zero(j)),j=1,4)
   !
   !
   ! Plotting from exact solution set
   
   !
end do
!
!
end program grsolve
