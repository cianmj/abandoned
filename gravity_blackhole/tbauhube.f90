!***************************************************************************
!*               Test program for Bauhuber's method                        *
!* ----------------------------------------------------------------------- *
!*  This program uses Bauhuber's Method to find all real or complex        *
!*  roots of a polynomial of degree n:                                     *
!***************************************************************************
Program TBauhube

Use FBauhube

IMPLICIT REAL*8 A-H,O-Z

Dimension ar(0:NMAX), ai(0:NMAX), rootr(0:NMAX), rooti(0:NMAX), val(0:NMAX)
Integer rc, skala 


  print *,' '
  print *,'  ----------------------------------------------------------'
  print *,'  Complex and Real Roots of a Polynomial (Bauhuber''s method)'
  print *,'  ----------------------------------------------------------'

  n = 4      !order of polynomial

  !define ar vector
  ar(0) =  0.000089248d0
  ar(1) = -0.04368d0
  ar(2) =  2.9948d0
  ar(3) = -6.798d0
  ar(4) =  1.d0

  !ai vector is null
  do i=0, n 
    ai(i) = 0.d0
  end do

  print *,'  Polynomial coefficients:'
  do i = 0, n
    write(*,10)  i, ar(i), ai(i)
  end do

  do j = 0, 1
    skala = j
    rc = bauhub (0, skala, n, ar, ai, rootr, rooti, val)
    if (rc == 0) then
      if (skala == 0) then
        print *,'  Roots (without scaling)'
      else
        print *,'  Roots (with scaling)'
      end if

      print *,'    No   Real part    Imaginary part   Function value'

      do i = 0, n-1
        write(*,20) i, rootr(i), rooti(i), val(i)
      end do
    else
      print *,'  *** Error in bauhube, rc=', rc
    end if
  end do

  print *,'  ----------------------------------------------------------'
  print *,' '
  stop

10 Format('   a(',I2,') = ', E13.6,'  ',E13.6)
20 Format('   ',I4,'  ',E13.6,'  ',E13.6,'  ',E13.6)

END

!end of file tbauhube.f90
