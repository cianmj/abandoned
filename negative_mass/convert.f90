program convert
  implicit none
  integer   :: i,inout,kr=1,xmod,seed(1)
  integer*8 :: NN
  double precision :: dth,xx,yy,teps
!  double complex :: pr,pr2,dr,r,mod1,mod2,prob(2),expon,wavefn
  double precision, allocatable :: xE0(:),yE0(:),rr(:),tt(:), &
       & wwr(:),wwi(:), wfsq(:), rand(:)
  double precision, parameter :: pi = 3.14159265358979325d0
!  double complex, allocatable :: pt(:,:),prabs(:,:),Matrix(:,:,:), &
!       & msum(:,:),mvec(:),ekvec(:), Mmin(:)
!  double precision, external :: modsq
  !
  write(*,*) 'Converting data from spherical to cartesian coordinates...'
  !
!  NN = (2000-1)**2
!  NN = 999000
!  NN = 1999000
NN = 3998000
! NN = 7998000
! NN = 7996000
! NN = 15996000
!
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
  !
  in_out:do inout=1,2
     !
     write(*,*) '#', inout
     if (inout==1) open(11,file="wav_sum.txt",status='old',action='read')
     if (inout==2) open(11,file="wav_sum2.txt",status='old',action='read')
     do i=1,NN
        read(11,*) rr(i),tt(i),wwr(i),wwi(i),wfsq(i)
     end do
     close (11)
     !
     if (inout==1) open(12,file="nwav_sum.txt",status='replace',action='write')
     if (inout==1) open(15,file="absw.txt",status='replace',action='write')
     if (inout==2) open(12,file="nwav_sum2.txt",status='replace',action='write')
     if (inout==2) open(15,file="absw2.txt",status='replace',action='write')
     !
     teps = 0.d0
     wrt:do i=1,NN
        !
!        if (modulo(i,37).ne.1) cycle wrt
 !       xmod=int(dble(rand(i)+2.d0)*4.d0)
 !       if ((modulo(i,xmod).ne.1).and.(abs(tt(i)-pi).gt.pi/8.d0)) cycle wrt
        !
        xx = rr(i)*sin(tt(i))
        yy = -rr(i)*cos(tt(i))
        !
        if ( (xx.eq.0.d0).and.(yy.eq.0.d0) ) cycle wrt
        !
        !if ( yy > 2.d0 ) cycle wrt
        !if ( yy < 94.d0 ) cycle wrt
        if ( abs(yy) > .5d0 ) cycle wrt
        !
        if ( abs(xx) > .5d0 ) cycle wrt
        !if ( abs(xx) > 0.2d0 ) cycle wrt
        !
       !   write(12,fmt='(4F30.14)') xx,yy,log10(wwr(i)+1.d0),log10(wwi(i)+1.d0)
        write(12,fmt='(4F30.14)') xx,yy,wwr(i),wwi(i)
        write(15,fmt='(4F30.14)') xx,yy,wfsq(i),log10(wfsq(i)+1.d0)
        !        write(12,fmt='(4F30.14)') xx,yy,wwr(i),wwi(i) !,wfsq(i)
     end do wrt
     !
     close(12)
     !
  end do in_out
  open(unit=21,file='random.txt',status='replace',action='write')
  write(21,*) int(rand(NN+1)*1.d6)
  close(21)
  deallocate( rand )
  deallocate( rr, tt, wwr, wwi, wfsq )
  write(*,*) '...done'
  stop
  !
end program convert
