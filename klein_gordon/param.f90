 !
 module parameters
   implicit none
   !
   integer, parameter :: ik = kind(1)
   integer, parameter :: rk = kind(1d0)
   !
   integer(ik), parameter :: npts = 4000 ! #points in discretization
   integer(ik), parameter :: nsts = 10  ! num. of eigenstates calculated
   real(rk), parameter    :: mass = 0.004d0 ! effective mass of e-
   real(rk), parameter    :: kk = 2.d0/5.d0   ! harmonic constant
   !
   real(rk), parameter :: aut  = 2.418884326505d-17  ! a.u. time (s)
   real(rk), parameter :: aum  = 9.109382616d-31    ! electron mass (kg)
   real(rk), parameter :: aul  = 5.29172210818d-11  ! Bohr radius   (m)
   real(rk), parameter :: aueV = 2.7211d1           ! eV to Hartre
   !
   real(rk), parameter :: hbar = 1._rk
   real(rk), parameter :: c = 137._rk
   !
   real(rk), parameter    :: pi = 3.14159265358979_rk
   complex(rk), parameter :: Im = (0._rk,1._rk)      ! imaginary number i
   !
   integer(ik), parameter :: nwell = 2     ! Number of quantum wells
   integer(ik), parameter :: ne    = 1     ! Number of particles (electrons)
   real(rk) :: depth    = 1.d0 !.5  ! Well depth (eV)
   real(rk) :: w_length = 1.d0  ! lenght of well (m)
   real(rk) :: b_length = 1.d0  ! lenght of barrier (m)
   real(rk) :: w_ends   = 1.d0 !.1 - ! ends of well ( 2 x b_length)
   !
   integer(ik), parameter  :: npul = 2
   !
 end module parameters
 !
 module common
   use parameters
   implicit none
   !
   logical :: perturb
   integer(rk) :: nsts2
   real(rk) :: dx       ! discretization of space, dx (in a.u.)
   real(rk) :: s_length ! total length of the quantum well system
   real(rk) :: time     !
   real(rk) :: targ, epsln,epswf,epsnrg
   real(rk),allocatable :: dm2(:,:)      ! two-photon dipole matrix elements
   logical :: wrnorm = .false.
   logical :: wrnorm2 = .false.
   !
   type field_spec
      real(rk) :: power(npul),amp(npul),tau(npul), &
           & freq(npul),phase(npul),chirp(npul)
   end type field_spec
   type(field_spec) :: f
   !
   contains
     !
     subroutine normalize(n,c)
       implicit none
       integer(ik), intent(in)    :: n
       real(rk), intent(inout)    :: c(0:n)
       integer(ik)                :: j
       real(rk)                   :: csum
       csum = 0._rk
       do j=0,n-1
          csum = csum + dx*( c(j)*c(j) + c(j+1)*c(j+1) ) / 2._rk
       end do
       if (real(csum).ne.1.) then
          if ((real(csum)==0.).and.(wrnorm==.false.)) then
             write(*,*) 'Not normalized: all coeff=0'
             wrnorm = .true.
             return
          end if
          c(:) = c(:) / sqrt(csum)
       end if
       return
     end subroutine normalize
     !
     subroutine normalizec(n,c)
       implicit none
       integer(ik), intent(in)    :: n
       complex(rk), intent(inout)    :: c(0:n)
       integer(ik)                :: j
       real(rk)                   :: csum
       csum = 0._rk
       do j=0,n-1
          csum = csum + dx*( c(j)*conjg(c(j)) + c(j+1)*conjg(c(j+1)) ) / 2._rk
       end do
       if (real(csum).ne.1.) then
          if ((real(csum)==0.).and.(wrnorm==.false.)) then
             write(*,*) 'Not normalized: all coeff=0'
             wrnorm = .true.
             return
          end if
          c(:) = c(:) / sqrt(csum)
       end if
       return
     end subroutine normalizec
     !
     subroutine cnorm(n,c)
       implicit none
       integer, intent(in)        :: n
       complex(rk), intent(inout) :: c(n)
       integer(ik)                :: j
       real(rk)                   :: csum
       !
       csum = 0._rk
       do j=1,n
          csum = csum + conjg(c(j))*c(j)
       end do
       !
       if (real(csum).ne.1.) then
          if ((real(csum)==0.).and.(wrnorm2==.false.)) then
             write(*,*) 'Not normalized: all coeff=0'
             wrnorm2 = .true.
             return
          end if
          c(:) = c(:) / sqrt(csum)
       end if
       !
       return
     end subroutine cnorm
     !
 end module common
