program random
implicit none
integer :: i,seed(1),k=1
integer*8 :: N=10000000
double precision :: pi, avg
double precision,allocatable :: rand(:)
double complex :: Im=(0.d0,1.d0), sum, num
!
allocate( rand(N+1) )
pi = 4.d0*atan(1.d0)
!
open(unit=10,file='random.txt',status='old',action='read')
read(10,*) seed(1)
close(10)
!
CALL RANDOM_SEED
CALL RANDOM_SEED (SIZE=k)
CALL RANDOM_SEED (PUT=seed(1:k))
CALL RANDOM_NUMBER(rand)
!
sum = 0.d0
avg = 0.d0
do i=1,N
   num = exp(dcmplx(2.d0*pi*Im*rand(i)))
!   write(11,*) dble(num),dimag(num)
   sum = sum + num
   avg = avg + dble(rand(i))
end do
sum = sum/N
avg = avg/N
!
open(unit=10,file='random.txt',status='replace',action='write')
write(10,*) int(rand(N+1)*1.d9)
close(10)
!
write(*,*) sum,avg
!
end program random
