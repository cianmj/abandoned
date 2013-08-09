!****************************************************
!*   Find all roots of a complex polynomial using   *
!*   Laguerre formulation in complex domain.        *
!****************************************************
      SUBROUTINE TCLAGUE(A,imp,RAC)
      integer, parameter :: N=4

      COMPLEX*16 A(0:N), RAC(0:N), WORK(0:3*N)
      real*8 EPS,EPS2
      integer IMP, ITMAX

      ITMAX=20                  !Maximum number of iterations
      EPS=1.d-8; EPS2=1.d-6     !Minimum, maximum relative error
      RAC(0)=CMPLX(0.,0.)       !Approximate value of 1st root

      call CLAGUE(N,A,ITMAX,EPS,EPS2,IMP,RAC,WORK)

      RETURN
      END SUBROUTINE TCLAGUE

      SUBROUTINE CLAGUE(N,A,ITMAX,EPS,EPS2,IMP,X,W)

      COMPLEX*16 A(0:N),X(0:N),W(0:*)
      REAL*8 EPS,EPS2
      N1=N+1
      N2=N1+N
      CALL       LGUERR(N,A,ITMAX,EPS,EPS2,X,IMP,W(0),W(N1),W(N2))
      RETURN
      END

      SUBROUTINE LGUERR(K,A,ITMAX,EPS,EPS2,X,IMP,B,C,D)
      COMPLEX*16 A(0:K),B(0:K),C(0:K),D(0:K),X(0:K) &
      ,XK,XR,F,FP,FS,H,DEN,SQ,D2
      REAL*8 EPS,EPS2,TEST
      IMP=0
      B(0)=A(0)
      C(0)=B(0)
      D(0)=C(0)
      IF(K.EQ.1) THEN
      X(1)=-A(1)/A(0)
      RETURN
      ENDIF
      N=K

      XK=X(0)
    1 IK=0
      EPS1=EPS
      IT=0

    2 DO I=1,N
      B(I)=A(I)+XK*B(I-1)
      IF(I.LE.N-1) C(I)=B(I)+XK*C(I-1)
      IF(I.LE.N-2) D(I)=C(I)+XK*D(I-1)
      ENDDO

      F=B(N)
      FP=C(N-1)
      FS=D(N-2)
      H=(FLOAT(N-1)*FP)**2-FLOAT(N*(N-1))*F*FS
      IF(ABS(H).EQ.0.D0) THEN
      DO I=N,1,-1
      X(I)=-A(1)/(FLOAT(N)*A(0))
      ENDDO
      RETURN
      ENDIF
      IF(FP.EQ.0.D0) THEN
      DEN=SQRT(H)
      ELSE
      SQ=SQRT(H)
      DEN=FP+SQ
      D2=FP-SQ
      IF(ABS(DEN).LT.ABS(D2)) DEN=D2
      ENDIF
      IF(DEN.EQ.0.D0) THEN
      XK=XK+(0.1D0,0.1D0)
      GO TO 2
      ENDIF
      IK=IK+1
      XR=XK-FLOAT(N)*F/DEN
      TEST=ABS(XR-XK)
      XK=XR
      IF(TEST.LT.EPS1*ABS(XR)) GO TO 3
!     WRITE(*,'(1X,I4,5E15.6)') IK,H,DEN,XK,TEST
      IF(IK.LE.ITMAX)  GO TO 2
      EPS1=EPS1*10.D0
      IK=0
!     WRITE(*,'(1X,A,E15.6)') 'EPS1=',EPS1
      IF(EPS1.LT.EPS2) GO TO 2
      IMP=1
      RETURN

    3 X(N)=XR
      N=N-1
      IF(N.EQ.1) GO TO 4
      IF(N.LE.0) RETURN
      DO I=0,N
      A(I)=B(I)
      ENDDO
      GO TO 1

    4 X(N)=-B(1)/B(0)
      RETURN
      END

! End of file tclague.f90
