!****************************************************
!*   Find all roots of a complex polynomial using   *
!*   Newton's iterative formulation.                *
!****************************************************   
SUBROUTINE TEST_NEWTON(A,B,IERR,RR,RI)
parameter(N=4)

real*8  A(N+1), B(N+1), RR(N+1), RI(N+1)
real*8  WORK(4*N+5)
integer IERR

  call NEWTON (N,A,B,RR,RI,IERR,WORK)

  return

END SUBROUTINE TEST_NEWTON


      SUBROUTINE NEWTON (N,A,B,X,Y,IER,W)
      IMPLICIT REAL *8 (A-H,O-Z)
!
!     CALCULATE ALL THE ROOTS OF A COMPLEX POLYNOMIAL USING
!     NEWTON'S ITERATIVE FORMULATION

      DIMENSION A(*),B(*),X(*),Y(*),W(*)
      N1 = N+2
      N2 = N+N+3
      N3 = N+N+N+4
      CALL CNWTON (A,B,X,Y,N,IER,W(1),W(N1),W(N2),W(N3))
      RETURN
      END
!
      SUBROUTINE CNWTON (A,B,X,Y,N,IER,G,H,T,F)
      IMPLICIT REAL *8 (A-H,O-Z)
      DIMENSION A(*),B(*),X(*),Y(*),G(*),H(*),T(*),F(*)
      DATA IT1 /25/,IT2 /100/,TOL /1.D-2/,EPS /0.D0/
      IER = 0
      IF (EPS.NE.0.D0) GO TO 7
      EPS = 1.D0
    5 EPS = 0.5D0*EPS
      EPS1 = EPS+1.D0
      IF (EPS1.GT.1.D0) GO TO 5
    7 L = N
      X0 = 1.D0
      Y0 = 1.D0
      G(1) = A(1)
      T(1) = A(1)
      H(1) = B(1)
      F(1) = B(1)
   10 IF (L.EQ.1) GO TO 55
      M = L+1
      I = N-L+1
      K1 = 0
      K2 = 0
      X1 = X0
      Y1 = Y0
   15 DO 20 K = 2,M
      G(K) = A(K)+X1*G(K-1)-Y1*H(K-1)
   20 H(K) = B(K)+Y1*G(K-1)+X1*H(K-1)
      DO 25 K = 2,L
      T(K) = G(K)+X1*T(K-1)-Y1*F(K-1)
   25 F(K) = H(K)+Y1*T(K-1)+X1*F(K-1)
      GD = T(L)*T(L)+F(L)*F(L)
      GX = G(M)*T(L)+H(M)*F(L)
      GY = H(M)*T(L)-G(M)*F(L)
      DX = -GX/GD
      DY = -GY/GD
      X2 = X1+DX
      Y2 = Y1+DY
      R1 = (DABS(DX)+DABS(DY))/(DABS(X2)+DABS(Y2))
      K2 = K2+1
      X1 = X2
      Y1 = Y2
      IF (R1.LT.TOL) GO TO 30
      IF (K2.LE.IT2) GO TO 15
      IER = I
      GO TO 60
   30 K1 = K1+1
      IF (K1.LE.IT1) GO TO 15
      X(I) = X2
      Y(I) = Y2
      IF (DABS(X2).LE.EPS) X2 = EPS
      R2 = DABS(Y2/X2)
      IF (R2.LT.0.1D0) GO TO 40
      IF (R2*EPS.GT.1.D0) X(I) = 0.D0
      X0 = X2
      Y0 = Y2
      GO TO 45
   40 IF (R2.LE.EPS) Y(I) = 0.D0
      X0 = X2+Y2
      Y0 = X2-Y2
   45 DO 50 K = 1,L
      A(K) = G(K)
   50 B(K) = H(K)
      L = L-1
      GO TO 10
   55 R3 = G(1)*G(1)+H(1)*H(1)
      X(N) = -(G(1)*G(2)+H(1)*H(2))/R3
      Y(N) =  (H(1)*G(2)-G(1)*H(2))/R3
      IF (DABS(X(N)).LE.EPS) X(N) = EPS
      R2 = DABS(Y(N)/X(N))
      IF (R2*EPS.GT.1.D0) X(N) = 0.D0
      IF (R2.LE.EPS)      Y(N) = 0.D0
   60 CONTINUE
      DO L=1,N
      P=X(L)
      Q=Y(L)
      R=P*P+Q*Q
      IF(L.EQ.1) GO TO 70
      DO 65 I=L,2,-1
      CABS=X(I-1)**2+Y(I-1)**2
      IF(R.GE.CABS) GO TO 75
      X(I)=X(I-1)
      Y(I)=Y(I-1)
   65 CONTINUE
   70 I=1
   75 X(I)=P
      Y(I)=Q
      ENDDO
      END

!end of file tnewton.f90
