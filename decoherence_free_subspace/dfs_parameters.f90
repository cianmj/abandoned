 module dfs_parameters
   integer, parameter :: ik = kind(1)
   integer, parameter :: rk = kind(1d0)
   !
   !
   integer(ik) :: initial_states = 4       ! number of degenerate states in our initial subspace
   integer(ik) :: final_states   = 2       ! number of energetically degenerate final states
   integer(ik) :: transitions    = 1       ! number of distinct transitions from initial to final subspace
   !
   real(rk)           :: epsilon = 1E-15   ! Accuracy limit (if value<epsilon then value= 0._rk)
   !
   ! Verify whether we are capable of constucting a DFS of dimension
   ! ( initial_states - (final_states*transitions) ) >= 0 
   !
 end module dfs_parameters
