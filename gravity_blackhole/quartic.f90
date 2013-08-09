!
! General Solution to solve a quartic polynomial of the form:
! AA(4)x^4 + AA(3)x^3 + AA(2)x^2 + AA(1)x + AA(0)
! 
! Solution is known as Ferrari's Method.
! http://en.wikipedia.org/wiki/Quartic_equation
!
subroutine quartic(AA,RR)
  implicit none
  !
  double complex, intent(in)  :: AA(0:4)
  double complex, intent(out) :: RR(4)
  double precision, parameter :: eps = 1.d-12
  !
  double complex :: A,B,C,D,E
  double complex :: alpha, beta, gamma, P, Q, R, U, y, W
  !
  A=AA(4); B=AA(3); C=AA(2); D=AA(1); E=AA(0)
  !
  alpha = -3.d0*B**2./(8.d0*A**2.) + C/A
  beta  = B**3./(8.d0*A**3.) - B*C/(2.d0*A**2) + D/A
  gamma = -3.d0*B**4./(256.d0*A**4.) + C*B**2./(16.d0*A**3.) &
       & - B*D/(4.d0*A**2.) + E/A
  !
  if (abs(beta)<eps) then
     RR(1) = - B/(4.d0*A) + sqrt((-alpha + sqrt(alpha**2.-4.d0*gamma))/2.d0)
     RR(2) = - B/(4.d0*A) - sqrt((-alpha + sqrt(alpha**2.-4.d0*gamma))/2.d0)
     RR(3) = - B/(4.d0*A) + sqrt((-alpha - sqrt(alpha**2.-4.d0*gamma))/2.d0)
     RR(4) = - B/(4.d0*A) - sqrt((-alpha - sqrt(alpha**2.-4.d0*gamma))/2.d0)
     return
  end if
  !
  P = -alpha**2./12.d0 - gamma
  Q = -alpha**3./108.d0 + alpha*gamma/3.d0 - beta**2./8.d0
  R = Q/2.d0 + sqrt(Q**2./4.d0 + P**3./27.d0)
  !
  U = R**(1.d0/3.d0)
  !
  if (abs(P)<eps) then
     y = -5.d0*alpha/6.d0 - Q**(1.d0/3.d0)
  else
     y = -5.d0*alpha/6.d0 + P/(3.d0*U) - U
  end if
  !
  W = sqrt(alpha + 2.d0*y)
  !
  RR(1) = -B/(4.d0*A) + ( + W + sqrt(-(3.d0*alpha+2.d0*y + 2.d0*beta/W)) )/2.d0
  RR(2) = -B/(4.d0*A) + ( + W - sqrt(-(3.d0*alpha+2.d0*y + 2.d0*beta/W)) )/2.d0
  RR(3) = -B/(4.d0*A) + ( - W + sqrt(-(3.d0*alpha+2.d0*y - 2.d0*beta/W)) )/2.d0
  RR(4) = -B/(4.d0*A) + ( - W - sqrt(-(3.d0*alpha+2.d0*y - 2.d0*beta/W)) )/2.d0
  !
  return
  !
end subroutine quartic
