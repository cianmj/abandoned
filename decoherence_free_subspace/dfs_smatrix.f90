 module dfs_s_matrix
   use dfs_parameters
   implicit none
   private
   public dfs_calculation, eigensystem, linslv
   !
   contains
     !
     ! Given an S-matrix, this subroutine will output the coefficients such that the initial degenerate states are free
     ! from decoherence due to the interaction.
     !
     subroutine dfs_calculation(loop,dim,s_matrix,dim_dfs,coeff)
       use dfs_parameters
       implicit none
       !
       integer,     intent(in)  :: loop
       integer,     intent(in)  :: dim                                      ! Dimension of the S-matrix
       complex(rk), intent(in)  :: s_matrix(:,:,:)                          ! S-matrix that will be inputted into main program
       integer(ik)              :: i,j,k,l                                  ! summation indices
       integer(ik)              :: dim_dfs,dim_final,tdim                   ! dimiension of the DFS and final space
       complex(rk)              :: sigma(1:initial_states,1:initial_states) ! sigma matrix used for DFS calculation
       complex(rk), allocatable, intent(out) :: coeff(:,:)                  ! degenerate eigenvectors
       complex(rk), allocatable :: dfs_coeff(:,:)                           ! coefficients for DFS part
       complex(rk),    allocatable :: eval(:,:)
       real(rk)                 :: lambda(1:initial_states)                 ! Eigenvalues of sigma matrix
       complex(rk)                 :: e_states(1:initial_states,1:initial_states) ! Eigenvectors of sigma
       !
       integer(ik)              :: method = 0                               ! Determines which procedure we use to solve
                                                                            ! the sigma system of equations
       complex(rk), allocatable :: A(:,:),B(:)                              ! matrices for lin.alg. subroutine
       complex(rk)              :: work(1000)
       integer(ik)              :: info
       real(rk)                 :: rl,img
       !
       complex(rk), intrinsic   :: matmul
       intrinsic :: transpose, conjg, dimag, dble, abs, dcmplx
       external  :: dgels
       !
       dim_final = transitions*final_states
       dim_dfs   = initial_states - dim_final
       !
       if (loop.eq.1) allocate( coeff(1:initial_states,1:dim_dfs) )
       allocate (dfs_coeff(1:dim_final,1:dim_dfs),A(1:dim_final,1:dim_final), &
                 & B(1:dim_final),eval(1:initial_states,1:dim_dfs) )
       !
       ! Construction of the sigma matrix used to find the coefficients for the chaperon states.
       ! Requires a summation over the # of degenerate final states and the # of transitions.
       !
       !
       tdim = initial_states+1
       sigma(:,:) = 0._rk
       do l=1,transitions
         sigma(:,:) = sigma(:,:) + &
         &  matmul(transpose(conjg(s_matrix(tdim:dim,1:initial_states,l))),s_matrix(tdim:dim,1:initial_states,l))
       end do
       !
       ! Finding the eigenvalues of the sigma matrix.
       !
       call eigensystem(initial_states,sigma,lambda,e_states)
       !
       !
       coeff(:,:) = dcmplx(0._rk,0._rk)
       do i=1,dim_dfs
         coeff(:,i) = e_states(:,dim_final+i)
       end do
       !
       !
       ! Solving the linear equation (sigma-lamba)*coeff=0 for the dfs part of the coefficients
       ! This subroutine will solve the equation Ax=B, where A is the part of the S-matrix
       ! corresponding to the dfs part of the coefficients and B is related to the part
       ! of the S matrix corresponding to unity in the coefficient.
       !
      if ( (method.eq.1).OR.(method.eq.2)) then
         do i=1,initial_states
           sigma(i,i) = sigma(i,i) - lambda(i)
         end do
       end if
       !
       ! This first method solves the over-determined linear system of equations (sigma-lambda)*c = 0,
       ! using a LAPACK subroutine implementing least square analysis.
       !
        if (method.eq.1) then
         tdim=dim_dfs+1
         A(:,:) = sigma(:,tdim:initial_states)
         do i=1,dim_dfs
           B(:)   = -sigma(:,i)
           call dgels('N',initial_states,dim_final,1,A,initial_states &
                      &  ,B,initial_states,work,1000,info)
           dfs_coeff(:,i) = B(:)
           if (info.ne.0) write(*,*) 'info : ',info
         end do
       !
       ! This second method reduces the over-determined system to a select few equations such that it can be solved using
       ! standard matrix row reduction.  Then these answers are verified in the remaining equations.
       !
       else if (method.eq.2) then
         tdim = dim_dfs+1
         A(:,:) = sigma(1:dim_final,tdim:initial_states)
         do i=1,dim_dfs
           B(:)   = -sigma(1:dim_final,i)
           call linslv(dim_final,A,B)
           dfs_coeff(:,i) = B(:)
         end do
       end if
       !
       ! Here we construct the degenerate eigenvectors, which include the dfs part.
       !
       if ( (method.eq.1).OR.(method.eq.2)) then
         coeff(:,:) = dcmplx(0._rk, 0._rk)
         do i=1,dim_dfs
           coeff(i,i)= dcmplx(1._rk, 0._rk)
           coeff(tdim:initial_states,i) = dfs_coeff(:,i)
           do j=1,initial_states
             rl = dble(coeff(j,i)); img = dimag(coeff(j,i))
             if (abs(rl).lt.epsilon) then
               coeff(j,i) = dcmplx(0._rk, img)
               rl = 0._rk
             end if
             if (abs(img).lt.epsilon) coeff(j,i) = dcmplx(rl, 0._rk)
           end do
         end do
       end if
       !
       !
       ! Check the solutions to the DFS subspace coefficients in the extra equations of sigma
       !
       eval(:,:) = dcmplx(0._rk,0._rk)
       do j=1,dim_dfs
         eval(:,j) = matmul (sigma(:,:),coeff(:,j))
         do i=1,initial_states
           if ( (abs(dble(eval(i,j))).gt.epsilon).OR.(abs(dimag(eval(i,j))).gt.epsilon) ) then
             write(*,*) "Inaccurate solution to linear equation"
             write(*,*) eval(i,j)
           end if
         end do
       end do
       !
     end subroutine dfs_calculation
     !
     !
     ! Subroutine will compute the eigenvalues and eigenvectors for an n-by-n matrix A.
     !
     subroutine eigensystem(n,A,lambda,evec)
       use dfs_parameters
       implicit none
       !
       integer(ik), intent(in)       :: n
       complex(rk),    intent(inout) :: A(n,n)
       real(rk),    intent(out)      :: lambda(n)
       complex(rk),    intent(out)   :: evec(n,n)
       integer(ik)                   :: i
       integer(ik)                   :: ldvr,ldvl,info,lwork
       complex(rk)                   :: w(n), vl(n,n)
       real(rk), allocatable         :: rwork(:)
       complex(rk), allocatable      :: work(:)
       real(rk)                      :: l_sort, wi(n)
       complex(rk)                   :: e_sort(1:n)
       complex(rk)                   :: A_store(n,n)
       !
       A_store(:,:) = A(:,:)
       ldvr=n; ldvl=1
       lwork = 4*n
       !
       allocate( work(lwork),rwork(2*n) )
       !
       call zgeev ('N','V',n,A,n,w,vl,ldvl,evec,ldvr,work,lwork,rwork,info)
       !
       if (info.ne.0) then
         write(*,*) "Error occurred in eigensystem, info: ",info
         stop
       end if
       !
       lambda(:) = dble(w(:))
       wi(:)     = dimag(w(:))
       !
       do i=1,n
         if (abs(wi(i)).le.epsilon) then
           wi(i) = 0._rk
         else
           if (abs(w(i)).ne.0._rk) write(*,*) "Imaginary eigenvalue! :", wi(i)
         end if
         if (abs(lambda(i)).le.epsilon) then
           lambda(i) = 0._rk
         else
           if (abs(lambda(i)).le.1E-6) write(*,*) "Near zero eigenvalue :", lambda(i)
         end if
       end do
       !
       ! Ordering the eigenvalues ( lambda ) from greatest to smallest, and their corresponding eigenvectors.
       !
       sort: do i=1,n-1
         if (lambda(i).lt.lambda(i+1)) then
           if (lambda(i+1).eq.0_rk) cycle sort
           l_sort      = lambda(i)   ; e_sort(:)   = evec(:,i)
           lambda(i)   = lambda(i+1) ; evec(:,i)   = evec(:,i+1)
           lambda(i+1) = l_sort      ; evec(:,i+1) = e_sort(:)
         end if
       end do sort
       !
       A(:,:) = A_store(:,:)
       !
     end subroutine eigensystem
     !
     !
     ! This subroutine will solve the system of 'n' linear complex equations of the form Ax=B using
     ! a LAPACK subroutine.
     !
     subroutine linslv(n,A,B)
       use dfs_parameters
       implicit none
       !
       integer(ik), intent(in)    :: n
       complex(rk), intent(in)    :: A(1:n,1:n)
       complex(rk), intent(inout) :: B(1,n)
       integer(ik)                :: nrhs,lda,ldb,info
       integer(ik)                :: ipiv(1:n)
       !
       nrhs = 1;
       lda = max(1,n)
       ldb = max(1,n)
       !
       call zgesv(n,nrhs,A,lda,ipiv,B,ldb,info)
       !
       if (info.ne.0) then
         write(*,*) "Problem solving equations, info =",info
         stop
       end if
       !
     end subroutine linslv
     !
     !
     !
 end module dfs_s_matrix
