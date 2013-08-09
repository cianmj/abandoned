 program dfs_main
   use dfs_parameters
   use dfs_s_matrix
   implicit none
   !
   integer(ik)              :: dim,dim_dfs            ! dimension of the S-matrix in question
   integer(ik)              :: i,j,k
   integer(ik)              :: loop = 1
   real(rk), allocatable    :: r_matrix (:,:,:)    ! r-matrix that will be inputted into main program
   complex(rk), allocatable :: s_matrix (:,:,:)    ! contruct this unitary matrix from the input
   complex(rk), allocatable :: coeff    (:,:)      ! Coefficients for Eigenvectors
   complex(rk), allocatable :: initial_coeff(:,:)
   complex(rk), allocatable :: full_coeff(:,:), vsum(:)
   !
   intrinsic :: transpose, conjg
   external  :: unitary
   !
   dim = initial_states + final_states
   allocate ( r_matrix(1:dim,1:dim,1:transitions),s_matrix(1:dim,1:dim,1:transitions) )
   !
   ! Input Reactance matrix (r_matrix) into the program.
   !
   open (10,file='r_matrix.10' ,status='old',action='read')
   do i=1,6
     read(10,*) r_matrix(i,1,1),r_matrix(i,2,1),r_matrix(i,3,1), & 
               & r_matrix(i,4,1),r_matrix(i,5,1),r_matrix(i,6,1)
   end do
   close(10)
   !
   ! Construct the unitary S-matrix from the hermitian reaction matrix (r_matrix)
   !
   do i=1,transitions
     call unitary(dim,r_matrix(:,:,i),s_matrix(:,:,i))
   end do
   !
   !
   ! Calculate the dfs for the given S-matrix, and construct the appropriate coefficients.
   !
   call dfs_calculation(loop,dim,s_matrix,dim_dfs,coeff)
   !
   allocate( initial_coeff(1:dim,1:dim_dfs), full_coeff(1:dim,1:dim_dfs), vsum(dim_dfs) )
   !
   full_coeff(1:initial_states,:) = coeff(:,:); full_coeff(initial_states+1:dim,:) = 0._rk
   !
   call normalize(dim,dim_dfs,full_coeff)
   !
   initial_coeff(:,:) = full_coeff(:,:)
   !
   write(*,*) "Initial coefficients : "
   write(*,*) initial_coeff(:,:)
   !
   ! Act on the initial state of the system using the S-matrix, iterate.
   !
   do while(.true.) 
     loop = loop + 1
     !
     write(*,*) "before"; write(*,*) initial_coeff(:,:)
     vsum(:) = dcmplx(0._rk,0._rk)
     do i=1,dim_dfs
       do k=1,transitions
!         vsum(:) = matmul(s_matrix(:,:,k),initial_coeff(:,i))
!         full_coeff(:,i) = vsum(:)
!         full_coeff(:,i) = full_coeff(:,i) + matmul(s_matrix(:,:,k),initial_coeff(:,i))
         full_coeff(:,i) = matmul(s_matrix(:,:,k),initial_coeff(:,i))
!         vsum(:) = dcmplx(0._rk,0._rk)
       end do
     end do
     !
     call normalize(dim,dim_dfs,full_coeff(:,:))
     !
     write(*,*) "after"; write(*,*) full_coeff(:,:)
     read(*,*)
     !
!     if (modulo(loop,20).eq.0) then
!       write(*,*) 'loop : ',loop
!       full_coeff(1:initial_states,:) = coeff(:,:)
!       initial_coeff(:,:) = full_coeff(:,:)
!     end if
     !
     initial_coeff(:,:) = full_coeff(:,:)     
     !
   end do
   !
 end program dfs_main
!
!
! Takes in a Hermitian matrix A, and computes the associated unitary matrix B.
!
 subroutine unitary(n,A,B)
   use dfs_parameters
   implicit none
   !
   integer(ik), intent(in)     :: n
   real(rk),    intent(in)     :: A(n,n)
   complex(rk), intent(out)    :: B(n,n)
   !
   integer(ik)              :: i,j
   complex(ik)              :: identity(n,n)
   complex(rk)              :: mat_1(n,n),mat_2(n,n)
   integer(ik)              :: lwork, info
   integer(ik)              :: ipiv(n)
   complex(rk), allocatable :: work(:)
   real(rk)                 :: u_check(n,n)
   !
   intrinsic   :: dcmplx, matmul
   external    :: zgetrf, zgetri
   !
   lwork = n*n
   allocate ( work(lwork) )
   !
   identity(:,:) = 0._rk
   do i=1,n
    identity(i,i) = dcmplx(1._rk,0._rk)
   end do
   !
   mat_1(:,:) = identity(:,:) + dcmplx(0._rk,A(:,:))
   mat_2(:,:) = identity(:,:) - dcmplx(0._rk,A(:,:))
   !
   call zgetrf(n,n,mat_2,n,ipiv,info)
   if (info.ne.0) write(*,*) "Problem factorizing matrix. info=",info
   !
   call zgetri(n,mat_2,n,ipiv,work,lwork,info)
   if (info.ne.0) write(*,*) "Problem finding unitary matrix. info=",info
   !
   B(:,:) = matmul( mat_2(:,:),mat_1(:,:) )
   !
   call approx(n,B)
   !
   ! Check if the B matrix is truely Unitary
   !
   u_check(:,:) = matmul( transpose(conjg(B(:,:))),B(:,:) )
   do i=1,n
     if (real(u_check(i,i)).ne.1._rk) then
       write(*,*) "Matrix is not unitary.", i
       stop
     end if
     do j=1,n
       if (i.ne.j) then
         if (u_check(i,j).gt.epsilon) then 
           write(*,*) "Matrix not unitary"
           stop
         end if
       end if
     end do
   end do
   !
 end subroutine unitary
!
!
!
 subroutine approx(n,A)
   use dfs_parameters
   implicit none
   integer(ik), intent(in)     :: n
   complex(rk), intent(inout)  :: A(n,n)
   integer(ik) :: i,j
   real(rk)    :: rt, ct
   !
   intrinsic :: dble, dimag, abs, dcmplx
   !
   do i=1,n
     do j=1,n
       rt = dble(A(i,j)); ct = dimag(A(i,j))
       if(abs(rt).lt.epsilon) rt = 0._rk
       if(abs(ct).lt.epsilon) ct = 0._rk
       A(i,j) = dcmplx(rt,ct)
     end do
   end do
   !
 end subroutine approx
!
!
! Normalization subroutine, used to check if the coefficients remain normalized.
!
 subroutine normalize(n,m,c)
   use dfs_parameters
   implicit none
   integer(ik), intent(in)    :: n,m
   complex(rk), intent(inout) :: c(n,m)
   integer(ik)                :: i,j
   real(rk)                   :: csum
   !
   do i=1,m
     csum = 0._rk
     do j=1,n
       csum = csum + conjg(c(j,i))*c(j,i)
     end do
     if (real(csum).ne.1.0) then
       write(*,*) "coefficients have lost their normalization : ", csum
       c(:,i) = c(:,i) / sqrt(csum)
     end if
   end do
   !
 end subroutine normalize
!
