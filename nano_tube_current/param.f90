 !
 module parameters
   implicit none
   !
   integer, parameter :: ik = kind(1)
   integer, parameter :: rk = kind(1d0)
   !
   real(rk), parameter :: aut  = 2.418884326505d-17  ! a.u. time (s)
   real(rk), parameter :: aum  = 9.109382616d-31    ! electron mass (kg)
   real(rk), parameter :: aul  = 5.29172210818d-11  ! Bohr radius   (m)
   real(rk), parameter :: aueV = 2.7210d1           ! eV to Hartre
   !
   real(rk), parameter    :: pi = 3.14159265358979_rk
   complex(rk), parameter :: Im = (0.,1.)           ! imaginary number i
   !
   integer(rk)            :: nwell = 4     ! Number of quantum wells
   integer(rk), parameter :: ne    = 1     ! Number of particles (electrons)
   real(rk), parameter    :: mass = 1._rk  ! mass of electron
   real(rk) :: depth    = 0.1_rk  ! Well depth (eV)
   real(rk) :: w_length = 10.d-9 ! lenght of well (m)
   real(rk) :: b_length = 2.d-9 ! lenght of barrier (m)
   real(rk) :: w_ends   = 5.d-9 ! ends of well
   !
   integer(ik), parameter :: npts = 10000 ! num. of points in discretization
   integer(ik), parameter :: nsts = 3   ! num. of eigenstates calculated
   !
 end module parameters
 !
 module common
   use parameters
   implicit none
   !
   real(rk) :: dx       ! discretization of space, dx (in a.u.)
   real(rk) :: s_length ! total length of the quantum well system
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
       do j=0,n
          csum = csum + c(j)*c(j)
!          csum = csum + conjg(c(j))*c(j)
       end do
       if (real(csum).ne.1.) then
          c(:) = c(:) / sqrt(csum)
       end if
       return
     end subroutine normalize
     !
 end module common
